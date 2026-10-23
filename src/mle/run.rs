//! Serial RUN1 orchestration and failure-safe publication.

use std::collections::HashSet;
use std::fmt::Write as _;
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Error, ErrorKind, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::config::MleConfig;
use crate::data::ParsedCountTable;
use crate::error::MageckError;

use super::design_matrix::DesignMatrix;

pub const DEFAULT_GENES_VARMODELING: usize = 0;
pub const DEFAULT_MLE_ALPHA: f64 = 0.01;
pub const MAX_TEMP_CREATE_ATTEMPTS: usize = 32;

static TMP_COUNTER: AtomicU64 = AtomicU64::new(0);

#[cfg(test)]
std::thread_local! {
    // Isolated to the calling unit-test thread; compiled out of ordinary CLI builds.
    static RECONCILIATION_PUBLICATION_FAULT: std::cell::Cell<u8> = const { std::cell::Cell::new(0) };
}

trait Filesystem {
    fn create_new(&self, path: &Path) -> std::io::Result<File>;
    fn rename(&self, from: &Path, to: &Path) -> std::io::Result<()>;
    fn remove_file(&self, path: &Path) -> std::io::Result<()>;
    fn write_temp(&self, file: File, data: &str) -> std::io::Result<()> {
        let mut writer = BufWriter::new(file);
        writer.write_all(data.as_bytes())?;
        writer.flush()
    }
}

struct RealFilesystem;

impl Filesystem for RealFilesystem {
    fn create_new(&self, path: &Path) -> std::io::Result<File> {
        OpenOptions::new().write(true).create_new(true).open(path)
    }
    fn rename(&self, from: &Path, to: &Path) -> std::io::Result<()> {
        fs::rename(from, to)
    }
    fn remove_file(&self, path: &Path) -> std::io::Result<()> {
        fs::remove_file(path)
    }
}

#[derive(Default)]
struct PublicationTracker {
    gene_tmp: Option<PathBuf>,
    sgrna_tmp: Option<PathBuf>,
    gene_target_published: bool,
    sgrna_target_published: bool,
}

fn temp_path(target: &Path) -> Result<PathBuf, MageckError> {
    let name = target.file_name().ok_or_else(|| {
        MageckError::Config(format!(
            "Output target has no file name: {}",
            target.display()
        ))
    })?;
    Ok(target.with_file_name(format!(
        "{}.tmp.{}.{}",
        name.to_string_lossy(),
        std::process::id(),
        TMP_COUNTER.fetch_add(1, Ordering::Relaxed)
    )))
}

fn create_temp<F: Filesystem>(fs: &F, target: &Path) -> Result<(PathBuf, File), MageckError> {
    for _ in 0..MAX_TEMP_CREATE_ATTEMPTS {
        let candidate = temp_path(target)?;
        match fs.create_new(&candidate) {
            Ok(file) => return Ok((candidate, file)),
            Err(error) if error.kind() == ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(MageckError::Io(error)),
        }
    }
    Err(MageckError::Io(Error::new(ErrorKind::AlreadyExists, format!(
        "Failed to create temporary output file after {} attempts for {}; no final output was published",
        MAX_TEMP_CREATE_ATTEMPTS, target.display()
    ))))
}

fn cleanup<F: Filesystem>(
    fs: &F,
    tracker: &PublicationTracker,
    gene_target: &Path,
    primary: Error,
) -> MageckError {
    let mut failures = Vec::new();
    for path in [
        tracker.gene_tmp.as_deref(),
        tracker.sgrna_tmp.as_deref(),
        tracker.gene_target_published.then_some(gene_target),
    ]
    .into_iter()
    .flatten()
    {
        if let Err(error) = fs.remove_file(path) {
            failures.push((path, error));
        }
    }
    if failures.is_empty() {
        return MageckError::Io(primary);
    }
    let mut message = format!("Primary publication error: {primary}");
    for (path, error) in failures {
        message.push_str(&format!(
            "; cleanup failed for '{}': {error}",
            path.display()
        ));
    }
    MageckError::Io(Error::new(primary.kind(), message))
}

fn publish_with<F: Filesystem>(
    fs: &F,
    prefix: &str,
    gene: &str,
    sgrna: &str,
) -> Result<(), MageckError> {
    let gene_target = PathBuf::from(format!("{prefix}.gene_summary.txt"));
    let sgrna_target = PathBuf::from(format!("{prefix}.sgrna_summary.txt"));
    if fs::symlink_metadata(&gene_target).is_ok() {
        return Err(MageckError::Config(format!(
            "Output file already exists: {}",
            gene_target.display()
        )));
    }
    if fs::symlink_metadata(&sgrna_target).is_ok() {
        return Err(MageckError::Config(format!(
            "Output file already exists: {}",
            sgrna_target.display()
        )));
    }
    if let Some(parent) = gene_target.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut tracker = PublicationTracker::default();
    let (gene_tmp, gene_file) = create_temp(fs, &gene_target)?;
    tracker.gene_tmp = Some(gene_tmp);
    let (sgrna_tmp, sgrna_file) = match create_temp(fs, &sgrna_target) {
        Ok(value) => value,
        Err(MageckError::Io(error)) => return Err(cleanup(fs, &tracker, &gene_target, error)),
        Err(error) => return Err(error),
    };
    tracker.sgrna_tmp = Some(sgrna_tmp);
    let write = |file: File, data: &str| -> std::io::Result<()> {
        let mut writer = BufWriter::new(file);
        writer.write_all(data.as_bytes())?;
        writer.flush()
    };
    if let Err(error) = write(gene_file, gene).and_then(|_| write(sgrna_file, sgrna)) {
        return Err(cleanup(fs, &tracker, &gene_target, error));
    }
    let gene_tmp = tracker.gene_tmp.as_ref().expect("tracked gene temporary");
    if let Err(error) = fs.rename(gene_tmp, &gene_target) {
        return Err(cleanup(fs, &tracker, &gene_target, error));
    }
    tracker.gene_tmp = None;
    tracker.gene_target_published = true;
    let sgrna_tmp = tracker.sgrna_tmp.as_ref().expect("tracked sgrna temporary");
    if let Err(error) = fs.rename(sgrna_tmp, &sgrna_target) {
        return Err(cleanup(fs, &tracker, &gene_target, error));
    }
    tracker.sgrna_tmp = None;
    tracker.sgrna_target_published = true;
    Ok(())
}

// CNV-only extension of RUN1: all three temporaries are written/flushed before
// sequential renames. The default two-output path above remains unchanged.
fn publish_cnv_with<F: Filesystem>(
    filesystem: &F,
    prefix: &str,
    gene: &str,
    sgrna: &str,
    info: &str,
) -> Result<(), MageckError> {
    let targets = ["gene_summary.txt", "sgrna_summary.txt", "cnv_run_info.txt"]
        .map(|suffix| PathBuf::from(format!("{prefix}.{suffix}")));
    for target in &targets {
        if fs::symlink_metadata(target).is_ok() {
            return Err(MageckError::Config(format!(
                "Output file already exists: {}",
                target.display()
            )));
        }
    }
    if let Some(parent) = targets[0].parent() {
        fs::create_dir_all(parent)?;
    }
    let mut temporaries = Vec::new();
    let mut published = Vec::new();
    let result = (|| {
        let mut files = Vec::new();
        for target in &targets {
            let (path, file) = create_temp(filesystem, target)?;
            temporaries.push(Some(path));
            files.push(file);
        }
        for (file, data) in files.into_iter().zip([gene, sgrna, info]) {
            filesystem.write_temp(file, data)?;
        }
        for (temporary, target) in temporaries.iter_mut().zip(&targets) {
            filesystem.rename(temporary.as_ref().expect("tracked temporary"), target)?;
            *temporary = None;
            published.push(target);
        }
        Ok(())
    })();
    match result {
        Err(MageckError::Io(primary)) => {
            let mut message = format!("Primary publication error: {primary}");
            let mut failed = false;
            for path in temporaries
                .iter()
                .filter_map(Option::as_ref)
                .chain(published)
            {
                if let Err(error) = filesystem.remove_file(path) {
                    failed = true;
                    write!(
                        message,
                        "; cleanup failed for '{}': {error}",
                        path.display()
                    )
                    .unwrap();
                }
            }
            Err(MageckError::Io(if failed {
                Error::new(primary.kind(), message)
            } else {
                primary
            }))
        }
        other => other,
    }
}

fn publish_diagnostics(target: &Path, data: &str) -> Result<(), MageckError> {
    if fs::symlink_metadata(target).is_ok() {
        return Err(MageckError::Config(format!(
            "Output file already exists: {}",
            target.display()
        )));
    }
    if let Some(parent) = target.parent() {
        fs::create_dir_all(parent)?;
    }
    let (temporary, file) = create_temp(&RealFilesystem, target)?;
    let result = RealFilesystem
        .write_temp(file, data)
        .and_then(|_| RealFilesystem.rename(&temporary, target));
    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result.map_err(MageckError::Io)
}

fn replace_owned_diagnostics(target: &Path, data: &str) -> Result<(), MageckError> {
    let (temporary, file) = create_temp(&RealFilesystem, target)?;
    let result = RealFilesystem
        .write_temp(file, data)
        .and_then(|_| RealFilesystem.rename(&temporary, target));
    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result.map_err(MageckError::Io)
}
fn controls(config: &MleConfig) -> Result<Option<Vec<String>>, MageckError> {
    if config.control_sgrna.is_none() {
        return Ok(None);
    }
    let path = config
        .control_sgrna
        .as_ref()
        .expect("validated control normalization configuration");
    let content = fs::read_to_string(path)?;
    let mut seen = HashSet::new();
    Ok(Some(
        content
            .lines()
            .map(str::trim)
            .filter(|id| !id.is_empty())
            .filter(|id| seen.insert((*id).to_owned()))
            .map(str::to_owned)
            .collect(),
    ))
}

fn training_window(requested: usize, available: usize) -> usize {
    if requested == 0 {
        0
    } else {
        requested.saturating_add(1).min(available)
    }
}

/// Runs the complete, deterministic, serial classical-MLE pipeline.
pub fn run(config: &MleConfig) -> Result<(), MageckError> {
    let mut profile = crate::execution::Profile::new("mle");
    crate::config::validate_test_threads(config.threads)?;
    if config.seed.is_some() {
        return Err(MageckError::Config(
            "PRNG seed configuration (--seed) is not supported in serial RUN1".to_owned(),
        ));
    }
    let design = DesignMatrix::from_reader(BufReader::new(File::open(&config.design_matrix)?))?;
    let count_table =
        ParsedCountTable::from_reader(BufReader::new(File::open(&config.count_table)?))?;
    let mut bundle = super::ingest::ingest_mle_input(&count_table, &design)?;
    drop(count_table);
    profile.stage("parse-ingest");
    let control_sgrnas = controls(config)?;
    let size_factors = super::normalize::compute_mle_size_factors_for_bundle(
        &bundle,
        config.norm_method,
        control_sgrnas.as_deref(),
    )?;
    // Keep every guide in normalization and retain the original file order for
    // variance training; oversized genes are excluded only from final fitting.
    let all_genes = bundle.genes.clone();
    let mut exclusions = Vec::new();
    bundle.genes.retain(|gene| {
        let control = gene.gene_id == "Control"
            || control_sgrnas
                .as_ref()
                .is_some_and(|ids| gene.sgrna_ids.iter().all(|g| ids.contains(g)));
        let reason = if control {
            Some("CONTROL_ONLY")
        } else if gene.sgrna_ids.len() > config.max_sgrnas_per_gene {
            Some("PREDECLARED_GUIDE_LIMIT")
        } else {
            None
        };
        if let Some(reason) = reason {
            if reason == "PREDECLARED_GUIDE_LIMIT" {
                eprintln!(
                    "Skipping MLE gene '{}': {} guides exceeds --max-sgrnas-per-gene threshold {}",
                    gene.gene_id,
                    gene.sgrna_ids.len(),
                    config.max_sgrnas_per_gene
                );
            }
            exclusions.push(super::diagnostics::Record::new(
                gene,
                "selection",
                None,
                Some(reason.into()),
            ));
            false
        } else {
            true
        }
    });
    let planned: Vec<_> = bundle.genes.iter().map(|g| g.gene_id.clone()).collect();
    let diagnostics_target =
        PathBuf::from(format!("{}.mle_diagnostics.json", config.output_prefix));
    for suffix in [
        "mle_diagnostics.json",
        "gene_summary.txt",
        "sgrna_summary.txt",
    ] {
        let path = PathBuf::from(format!("{}.{suffix}", config.output_prefix));
        if fs::symlink_metadata(&path).is_ok() {
            return Err(MageckError::Config(format!(
                "Output file already exists: {}",
                path.display()
            )));
        }
    }
    if config.cnv_norm.is_some() {
        let path = PathBuf::from(format!("{}.cnv_run_info.txt", config.output_prefix));
        if fs::symlink_metadata(&path).is_ok() {
            return Err(MageckError::Config(format!(
                "Output file already exists: {}",
                path.display()
            )));
        }
    }
    let identities = super::diagnostics::object(vec![
        ("normalization", format!("{:?}", config.norm_method).into()),
        ("genes_varmodeling", config.genes_varmodeling.into()),
        ("max_sgrnas_per_gene", config.max_sgrnas_per_gene.into()),
        (
            "control_sgrna_sha256",
            config
                .control_sgrna
                .as_ref()
                .map(|p| fs::read(p).map(|v| super::diagnostics::sha256(&v).into()))
                .transpose()?
                .unwrap_or(super::diagnostics::Json::Null),
        ),
        (
            "cnv_sha256",
            config
                .cnv_norm
                .as_ref()
                .map(|p| fs::read(p).map(|v| super::diagnostics::sha256(&v).into()))
                .transpose()?
                .unwrap_or(super::diagnostics::Json::Null),
        ),
        (
            "count_table_sha256",
            super::diagnostics::sha256(&fs::read(&config.count_table)?).into(),
        ),
        (
            "design_matrix_sha256",
            super::diagnostics::sha256(&fs::read(&config.design_matrix)?).into(),
        ),
        (
            "size_factors",
            super::diagnostics::Json::Array(size_factors.iter().map(|&x| x.into()).collect()),
        ),
    ]);
    let mut serial_input = bundle.as_serial_mle_input(
        &size_factors,
        0, // Set the training window against all genes below, including skipped genes.
        DEFAULT_MLE_ALPHA,
        true,
        false,
    )?;
    // Zero keeps scalar alpha=0.01. For positive limits the oracle checks after
    // adding the gene, so 2000 trains on 2001 genes in original file order.
    serial_input.genes_varmodeling = training_window(config.genes_varmodeling, all_genes.len());
    let cnv = config
        .cnv_norm
        .as_ref()
        .map(|path| {
            super::cnv_input::prepare(
                BufReader::new(File::open(path)?),
                bundle.design.condition_labels(),
                &bundle.genes,
            )
        })
        .transpose()?;
    let preparation_info = cnv
        .as_ref()
        .map(|_| {
            super::cnv_input::run_info(
                &bundle,
                &all_genes,
                &size_factors,
                config,
                serial_input.genes_varmodeling,
            )
        })
        .transpose()?;
    profile.stage("normalize-select-cnv-prepare");
    let execution = crate::execution::Execution::workflow(
        config.threads,
        bundle.genes.len().max(serial_input.genes_varmodeling),
    )?;
    let mut serial_output =
        match super::serial::run_mle_with_execution(&serial_input, &all_genes, &execution) {
            Ok(output) => output,
            Err(error) => {
                let records: Vec<_> = bundle
                    .genes
                    .iter()
                    .map(|g| {
                        super::diagnostics::Record::new(
                            g,
                            "validation",
                            None,
                            Some(error.to_string()),
                        )
                    })
                    .collect();
                let diagnostic = super::diagnostics::sidecar(
                    &records,
                    &[],
                    &planned,
                    &exclusions,
                    bundle.design.condition_labels(),
                    identities,
                    "INVALID_INPUT",
                );
                publish_diagnostics(&diagnostics_target, &diagnostic)?;
                return Err(error);
            }
        };
    profile.stage("fit");
    if serial_output.results.is_empty() && serial_output.failure.is_none() {
        serial_output.failure = Some("EMPTY_ELIGIBLE_UNIVERSE".into());
    }
    if let Some(error) = &serial_output.failure {
        let diagnostic = super::diagnostics::sidecar(
            &serial_output.records,
            &serial_output.results,
            &planned,
            &exclusions,
            bundle.design.condition_labels(),
            identities.clone(),
            if serial_output.results.is_empty() && serial_output.records.is_empty() {
                "EMPTY_ELIGIBLE_UNIVERSE"
            } else {
                "FAILED"
            },
        );
        publish_diagnostics(&diagnostics_target, &diagnostic)?;
        return Err(MageckError::Math(format!(
            "MLE planned-target failure: {error}; diagnostics only; no primary output published"
        )));
    }
    for r in &serial_output.results {
        if r.status == super::irls::IrlsStatus::ConvergedSameModelRetry {
            eprintln!(
                "MLE gene '{}' CONVERGED_SAME_MODEL_RETRY after {} attempts",
                r.gene_id, r.attempts
            );
        }
    }

    super::correction::adjust_mle_output_fdr(&mut serial_output)?;
    let cnv_info = cnv
        .as_ref()
        .map(|prepared| prepared.correct(&mut serial_output.results))
        .transpose()?;
    profile.stage("inference-cnv");
    let mut gene =
        super::output::serialize_gene_summary_to_string(&bundle.genes, &serial_output.results)?;
    // The serializer derives condition headers from fitted results. Retain the
    // full schema even when every gene was skipped and there are no results.
    if bundle.genes.is_empty() {
        gene.pop();
        for condition in bundle.design.condition_labels() {
            write!(gene, "\t{condition}|beta\t{condition}|z\t{condition}|p-value\t{condition}|fdr\t{condition}|wald-p-value\t{condition}|wald-fdr").unwrap();
        }
        gene.push('\n');
    }
    let fitted_sgrna =
        super::output::serialize_sgrna_summary_to_string(&bundle.genes, &serial_output.results)?;
    let mut fitted_rows = fitted_sgrna.lines().skip(1);
    let mut sgrna = String::from("Gene\tsgRNA\teff\n");
    for input in &all_genes {
        for guide in &input.sgrna_ids {
            if !planned.contains(&input.gene_id) {
                // Oracle convention: skipped guides retain eff=1 (unfitted,
                // default efficiency), observed for all 1,000 Brunello NTCs.
                writeln!(sgrna, "{}\t{guide}\t1", input.gene_id).unwrap();
            } else {
                sgrna.push_str(fitted_rows.next().expect("serialized fitted guide"));
                sgrna.push('\n');
            }
        }
    }
    profile.stage("serialize");
    // Retain the accepted auto metadata transaction while keeping the mandatory
    // numerical sidecar durable before any primary scientific output.
    let mut suffixes = vec!["gene_summary.txt", "sgrna_summary.txt", "run_info.txt"];
    if config.cnv_norm.is_some() {
        suffixes.push("cnv_run_info.txt");
    }
    let stage = if config.threads.is_some_and(crate::cli::ThreadMode::is_auto) {
        if let Some(parent) = Path::new(&config.output_prefix).parent() {
            fs::create_dir_all(parent)?;
        }
        Some(crate::run::OutputStage::new(
            &config.output_prefix,
            &suffixes,
        )?)
    } else {
        None
    };
    let output_prefix = stage
        .as_ref()
        .map_or(config.output_prefix.as_str(), |s| s.prefix.as_str());
    let diagnostic = super::diagnostics::sidecar(
        &serial_output.records,
        &serial_output.results,
        &planned,
        &exclusions,
        bundle.design.condition_labels(),
        identities.clone(),
        "READY_FOR_PUBLICATION",
    );
    // The mandatory sidecar is durable before any primary scientific file.
    publish_diagnostics(&diagnostics_target, &diagnostic)?;
    #[cfg(test)]
    if RECONCILIATION_PUBLICATION_FAULT.with(|fault| fault.get() == 1) {
        assert!(stage.is_some());
        assert!(fs::read_to_string(&diagnostics_target)?.contains("READY_FOR_PUBLICATION"));
        fs::write(
            format!("{}.run_info.txt", config.output_prefix),
            "OWNER SENTINEL\n",
        )?;
    }

    let result = (|| {
        match cnv_info {
            Some(mut info) => {
                info.push_str(
                    preparation_info
                        .as_deref()
                        .expect("CNV preparation metadata"),
                );
                publish_cnv_with(&RealFilesystem, output_prefix, &gene, &sgrna, &info)
            }
            None => publish_with(&RealFilesystem, output_prefix, &gene, &sgrna),
        }?;
        if let Some(stage) = &stage {
            fs::write(
                format!("{}.run_info.txt", stage.prefix),
                format!("key\tvalue\ncommand\tmle\n{}", execution.metadata()),
            )?;
            stage.publish(&config.output_prefix, &suffixes)?;
        }
        Ok(())
    })();
    let final_status = if result.is_ok() {
        "SUCCESS"
    } else {
        "PUBLICATION_FAILED"
    };
    let final_diagnostic = super::diagnostics::sidecar(
        &serial_output.records,
        &serial_output.results,
        &planned,
        &exclusions,
        bundle.design.condition_labels(),
        identities,
        final_status,
    );
    #[cfg(test)]
    if RECONCILIATION_PUBLICATION_FAULT.with(|fault| fault.get() == 2) {
        assert!(result.is_ok());
        assert!(Path::new(&format!("{}.gene_summary.txt", config.output_prefix)).exists());
        assert!(Path::new(&format!("{}.sgrna_summary.txt", config.output_prefix)).exists());
        if stage.is_some() {
            assert!(Path::new(&format!("{}.run_info.txt", config.output_prefix)).exists());
        }
        fs::remove_file(&diagnostics_target)?;
        fs::create_dir(&diagnostics_target)?;
    }
    // Replace only the sidecar this invocation just created; never an existing owner's file.
    if let Err(error) = replace_owned_diagnostics(&diagnostics_target, &final_diagnostic) {
        if result.is_ok() {
            for suffix in [
                "gene_summary.txt",
                "sgrna_summary.txt",
                "cnv_run_info.txt",
                "run_info.txt",
            ] {
                if (suffix != "cnv_run_info.txt" || config.cnv_norm.is_some())
                    && (suffix != "run_info.txt" || stage.is_some())
                {
                    fs::remove_file(format!("{}.{suffix}", config.output_prefix))?;
                }
            }
        }
        return Err(error);
    }
    profile.stage("publish");
    profile.finish();
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

    #[test]
    fn scalar_default_and_opt_in_activate_the_expected_model() {
        assert_eq!(DEFAULT_GENES_VARMODELING, 0);
        assert_eq!(training_window(0, 2), 0);
        assert_eq!(training_window(1, 2), 2);
        assert_eq!(training_window(2000, 3000), 2001);
        assert_eq!(training_window(usize::MAX, 2), 2);
        assert_eq!(training_window(1, 0), 0);

        let counts = ParsedCountTable::from_reader(std::io::Cursor::new(
            "sgRNA\tGene\tbaseline\ttreat\nkeep1\tKEEP\t100\t80\nkeep2\tKEEP\t110\t90\nskip1\tOVERSIZED\t100\t95\nskip2\tOVERSIZED\t120\t110\nskip3\tOVERSIZED\t90\t80\n",
        )).unwrap();
        let design = DesignMatrix::from_reader(std::io::Cursor::new(
            "Samples\tbaseline\ttreatment\nbaseline\t1\t0\ntreat\t1\t1\n",
        ))
        .unwrap();
        let mut bundle = super::super::ingest::ingest_mle_input(&counts, &design).unwrap();
        let all_genes = bundle.genes.clone();
        bundle.genes.retain(|gene| gene.sgrna_ids.len() <= 2);
        for requested in [DEFAULT_GENES_VARMODELING, 1, usize::MAX] {
            let mut input = bundle
                .as_serial_mle_input(&[1., 1.], 0, DEFAULT_MLE_ALPHA, true, false)
                .unwrap();
            input.genes_varmodeling = training_window(requested, all_genes.len());
            let output = super::super::serial::run_mle_serial(&input, &all_genes).unwrap();
            assert_eq!(output.meanvar_model.is_some(), requested > 0);
            assert_eq!(output.results.len(), 1);
            assert_eq!(output.results[0].gene_id, "KEEP");
        }
    }

    struct CollisionFilesystem;
    impl Filesystem for CollisionFilesystem {
        fn create_new(&self, _: &Path) -> std::io::Result<File> {
            Err(Error::new(ErrorKind::AlreadyExists, "collision"))
        }
        fn rename(&self, _: &Path, _: &Path) -> std::io::Result<()> {
            unreachable!()
        }
        fn remove_file(&self, _: &Path) -> std::io::Result<()> {
            unreachable!()
        }
    }

    #[test]
    fn temporary_creation_exhaustion_has_the_pinned_error() {
        let target = Path::new("result.gene_summary.txt");
        let error = create_temp(&CollisionFilesystem, target).unwrap_err();
        match error {
            MageckError::Io(error) => {
                assert_eq!(error.kind(), ErrorKind::AlreadyExists);
                assert_eq!(error.to_string(), "Failed to create temporary output file after 32 attempts for result.gene_summary.txt; no final output was published");
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn cleanup_preserves_primary_error_kind_and_context() {
        struct FailingCleanup;
        impl Filesystem for FailingCleanup {
            fn create_new(&self, _: &Path) -> std::io::Result<File> {
                unreachable!()
            }
            fn rename(&self, _: &Path, _: &Path) -> std::io::Result<()> {
                unreachable!()
            }
            fn remove_file(&self, _: &Path) -> std::io::Result<()> {
                Err(Error::new(ErrorKind::PermissionDenied, "cleanup"))
            }
        }
        let tracker = PublicationTracker {
            gene_tmp: Some("gene.tmp".into()),
            sgrna_tmp: Some("sgrna.tmp".into()),
            gene_target_published: true,
            sgrna_target_published: false,
        };
        let error = cleanup(
            &FailingCleanup,
            &tracker,
            Path::new("gene.out"),
            Error::other("primary"),
        );
        match error {
            MageckError::Io(error) => {
                assert_eq!(error.kind(), ErrorKind::Other);
                assert_eq!(error.to_string(), "Primary publication error: primary; cleanup failed for 'gene.tmp': cleanup; cleanup failed for 'sgrna.tmp': cleanup; cleanup failed for 'gene.out': cleanup");
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    struct RenameFailure {
        call: AtomicUsize,
        fail_on: usize,
    }
    impl Filesystem for RenameFailure {
        fn create_new(&self, path: &Path) -> std::io::Result<File> {
            OpenOptions::new().write(true).create_new(true).open(path)
        }
        fn rename(&self, from: &Path, to: &Path) -> std::io::Result<()> {
            if self.call.fetch_add(1, AtomicOrdering::SeqCst) + 1 == self.fail_on {
                Err(Error::new(ErrorKind::PermissionDenied, "rename injected"))
            } else {
                fs::rename(from, to)
            }
        }
        fn remove_file(&self, path: &Path) -> std::io::Result<()> {
            fs::remove_file(path)
        }
    }

    #[test]
    fn first_and_second_rename_failures_leave_no_outputs() {
        for fail_on in [1, 2] {
            let prefix = std::env::temp_dir().join(format!(
                "mageckrs-run1-rename-{fail_on}-{}",
                std::process::id()
            ));
            let result = publish_with(
                &RenameFailure {
                    call: AtomicUsize::new(0),
                    fail_on,
                },
                &prefix.to_string_lossy(),
                "gene",
                "sgrna",
            );
            assert!(
                matches!(result, Err(MageckError::Io(error)) if error.kind() == ErrorKind::PermissionDenied)
            );
            assert!(!PathBuf::from(format!("{}.gene_summary.txt", prefix.display())).exists());
            assert!(!PathBuf::from(format!("{}.sgrna_summary.txt", prefix.display())).exists());
        }
    }
    struct PublicationFailure {
        create: AtomicUsize,
        write: AtomicUsize,
        rename: AtomicUsize,
        stage: &'static str,
        fail_on: usize,
        fail_cleanup: bool,
    }
    impl Filesystem for PublicationFailure {
        fn create_new(&self, path: &Path) -> std::io::Result<File> {
            if self.stage == "create"
                && self.create.fetch_add(1, AtomicOrdering::SeqCst) + 1 == self.fail_on
            {
                return Err(Error::new(ErrorKind::PermissionDenied, "create injected"));
            }
            RealFilesystem.create_new(path)
        }
        fn write_temp(&self, file: File, data: &str) -> std::io::Result<()> {
            if self.stage == "write"
                && self.write.fetch_add(1, AtomicOrdering::SeqCst) + 1 == self.fail_on
            {
                return Err(Error::new(
                    ErrorKind::PermissionDenied,
                    "write/flush injected",
                ));
            }
            RealFilesystem.write_temp(file, data)
        }
        fn rename(&self, from: &Path, to: &Path) -> std::io::Result<()> {
            if self.stage == "rename"
                && self.rename.fetch_add(1, AtomicOrdering::SeqCst) + 1 == self.fail_on
            {
                return Err(Error::new(ErrorKind::PermissionDenied, "rename injected"));
            }
            RealFilesystem.rename(from, to)
        }
        fn remove_file(&self, path: &Path) -> std::io::Result<()> {
            if self.fail_cleanup {
                return Err(Error::other("cleanup injected"));
            }
            RealFilesystem.remove_file(path)
        }
    }
    #[test]
    fn cnv_create_write_and_rename_failures_roll_back_all_three_outputs() {
        for stage in ["create", "write", "rename"] {
            for fail_on in 1..=3 {
                let dir = tempfile::tempdir().unwrap();
                let fs = PublicationFailure {
                    create: AtomicUsize::new(0),
                    write: AtomicUsize::new(0),
                    rename: AtomicUsize::new(0),
                    stage,
                    fail_on,
                    fail_cleanup: false,
                };
                let error = publish_cnv_with(
                    &fs,
                    dir.path().join("out").to_str().unwrap(),
                    "gene",
                    "sgrna",
                    "info",
                )
                .unwrap_err();
                assert!(
                    matches!(error, MageckError::Io(e) if e.kind() == ErrorKind::PermissionDenied)
                );
                assert_eq!(
                    fs::read_dir(dir.path()).unwrap().count(),
                    0,
                    "{stage} {fail_on}"
                );
            }
        }
    }
    #[test]
    fn cnv_third_rename_cleanup_errors_preserve_primary_and_all_paths() {
        let dir = tempfile::tempdir().unwrap();
        let fs = PublicationFailure {
            create: AtomicUsize::new(0),
            write: AtomicUsize::new(0),
            rename: AtomicUsize::new(0),
            stage: "rename",
            fail_on: 3,
            fail_cleanup: true,
        };
        let error = publish_cnv_with(
            &fs,
            dir.path().join("out").to_str().unwrap(),
            "gene",
            "sgrna",
            "info",
        )
        .unwrap_err();
        let MageckError::Io(error) = error else {
            panic!("wrong error");
        };
        assert_eq!(error.kind(), ErrorKind::PermissionDenied);
        let message = error.to_string();
        assert!(message.starts_with("Primary publication error: rename injected"));
        assert_eq!(message.matches("cleanup failed for").count(), 3);
        assert!(message.contains("cnv_run_info.txt.tmp."));
        assert!(message.contains("gene_summary.txt'"));
        assert!(message.contains("sgrna_summary.txt'"));
    }
    fn reconciliation_publication_config(
        root: &Path,
        mode: Option<crate::cli::ThreadMode>,
    ) -> MleConfig {
        fs::write(root.join("count.tsv"), "sgRNA\tGene\ts0\ts1\ts2\ts3\nn1\tNORMAL\t99\t119\t109\t129\nn2\tNORMAL\t149\t169\t159\t179\n").unwrap();
        fs::write(
            root.join("design.tsv"),
            "Samples\tbaseline\ttreatment\ns0\t1\t0\ns1\t1\t1\ns2\t1\t0\ns3\t1\t1\n",
        )
        .unwrap();
        MleConfig {
            count_table: root.join("count.tsv"),
            design_matrix: root.join("design.tsv"),
            output_prefix: root.join("result").display().to_string(),
            norm_method: crate::cli::NormMethod::None,
            max_sgrnas_per_gene: 100,
            genes_varmodeling: 0,
            control_sgrna: None,
            cnv_norm: None,
            threads: mode,
            seed: None,
        }
    }
    #[test]
    fn combined_auto_late_metadata_collision_preserves_sentinel_and_failed_diagnostics() {
        for maximum in [None, Some(1), Some(2), Some(8)] {
            let d = tempfile::tempdir().unwrap();
            let config = reconciliation_publication_config(
                d.path(),
                Some(crate::cli::ThreadMode::Auto { maximum }),
            );
            RECONCILIATION_PUBLICATION_FAULT.with(|fault| fault.set(1));
            let result = run(&config);
            RECONCILIATION_PUBLICATION_FAULT.with(|fault| fault.set(0));
            assert!(result.is_err());
            assert!(!d.path().join("result.gene_summary.txt").exists());
            assert!(!d.path().join("result.sgrna_summary.txt").exists());
            assert_eq!(
                fs::read_to_string(d.path().join("result.run_info.txt")).unwrap(),
                "OWNER SENTINEL\n"
            );
            let diagnostic =
                fs::read_to_string(d.path().join("result.mle_diagnostics.json")).unwrap();
            assert!(diagnostic.contains("PUBLICATION_FAILED"));
            assert!(!diagnostic.contains("READY_FOR_PUBLICATION"));
            assert!(diagnostic.contains("\"m\":0"));
            assert_eq!(fs::read_dir(d.path()).unwrap().count(), 4);
        }
    }
    #[test]
    fn combined_final_sidecar_replacement_failure_rolls_back_published_science() {
        use crate::cli::ThreadMode;
        for mode in [
            None,
            Some(ThreadMode::Fixed(1)),
            Some(ThreadMode::Fixed(2)),
            Some(ThreadMode::Fixed(4)),
            Some(ThreadMode::Fixed(8)),
            Some(ThreadMode::Auto { maximum: None }),
            Some(ThreadMode::Auto { maximum: Some(2) }),
        ] {
            let d = tempfile::tempdir().unwrap();
            let config = reconciliation_publication_config(d.path(), mode);
            RECONCILIATION_PUBLICATION_FAULT.with(|fault| fault.set(2));
            let result = run(&config);
            RECONCILIATION_PUBLICATION_FAULT.with(|fault| fault.set(0));
            assert!(result.is_err());
            for suffix in ["gene_summary.txt", "sgrna_summary.txt", "run_info.txt"] {
                assert!(!d.path().join(format!("result.{suffix}")).exists());
            }
            assert!(d.path().join("result.mle_diagnostics.json").is_dir());
            assert_eq!(fs::read_dir(d.path()).unwrap().count(), 3);
        }
    }
}
