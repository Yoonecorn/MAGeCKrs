//! Count-table comparison orchestration for the test CLI.

use std::collections::HashSet;
use std::ffi::CString;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Write};
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::cli::NormMethod;
use crate::config::TestConfig;
use crate::count::normalize::normalize;
use crate::data::{CountMatrix, ParsedCountTable};
use crate::error::MageckError;
use crate::stats;

use super::aggregation::{aggregate_lfc, aggregate_results};
use super::rra::control_pool::prepare_control_pools;
use super::rra::{
    apply_directional_bh, compute_empirical_with_execution, compute_gene_lo_values, NullModelConfig,
};
use super::serialize::write_rra_summaries;
use super::validity_mask::build_comparison_validity;
use super::{calculate_sgrna_pvalues, fit_meanvar, CountTable};

const MAX_TEMP_CREATE_ATTEMPTS: usize = 32;
static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

fn target(prefix: &str, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{prefix}.{suffix}"))
}

fn reject_existing(path: &Path) -> Result<(), MageckError> {
    match fs::symlink_metadata(path) {
        Ok(_) => Err(MageckError::Config(format!(
            "Output already exists: {}",
            path.display()
        ))),
        Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(e.into()),
    }
}

fn run_info(config: &TestConfig, genes: usize, guides: usize) -> String {
    format!("key\tvalue\nmageckrs_version\t{}\ncommand\ttest\nrra_seed\t{:#018x}\nnum_genes\t{genes}\nnum_guides\t{guides}\ncontrol_null\t{}\noutput_prefix\t{}\n",
        env!("CARGO_PKG_VERSION"), config.rra_seed, config.control_sgrna.is_some(), config.output_prefix)
}

pub(crate) fn exclusive_rename(from: &Path, to: &Path) -> io::Result<()> {
    unsafe extern "C" {
        fn renamex_np(
            from: *const std::ffi::c_char,
            to: *const std::ffi::c_char,
            flags: std::ffi::c_uint,
        ) -> std::ffi::c_int;
    }
    let path = |p: &Path| {
        CString::new(p.as_os_str().as_bytes())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "NUL in output path"))
    };
    let from = path(from)?;
    let to = path(to)?;
    // SAFETY: both NUL-terminated strings live through the call. macOS SDK
    // sys/stdio.h defines 0x00000004 as RENAME_EXCL, refusing existing targets.
    if unsafe { renamex_np(from.as_ptr(), to.as_ptr(), 0x00000004) } == 0 {
        Ok(())
    } else {
        Err(io::Error::last_os_error())
    }
}

fn write_metadata(path: &Path, contents: &str) -> Result<(), MageckError> {
    reject_existing(path)?;
    for _ in 0..MAX_TEMP_CREATE_ATTEMPTS {
        let mut name = path
            .file_name()
            .ok_or_else(|| MageckError::Config("Metadata path has no filename".into()))?
            .to_os_string();
        name.push(format!(
            ".tmp.{}.{}",
            std::process::id(),
            TEMP_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        let temporary = path.with_file_name(name);
        let file = match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
        {
            Ok(file) => file,
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e.into()),
        };
        let result = (|| -> io::Result<()> {
            let mut writer = BufWriter::new(file);
            writer.write_all(contents.as_bytes())?;
            writer.flush()?;
            drop(writer);
            exclusive_rename(&temporary, path)
        })();
        if let Err(primary) = result {
            let mut message = primary.to_string();
            if let Err(cleanup) = fs::remove_file(&temporary) {
                message.push_str(&format!(
                    "; cleanup failed for {}: {cleanup}",
                    temporary.display()
                ));
            }
            return Err(MageckError::Io(io::Error::new(primary.kind(), message)));
        }
        return Ok(());
    }
    Err(io::Error::new(
        io::ErrorKind::AlreadyExists,
        "Metadata temporary creation exhausted 32 attempts",
    )
    .into())
}

fn rollback_summaries(prefix: &str, primary: MageckError) -> MageckError {
    let mut message = primary.to_string();
    let mut failed = false;
    for suffix in ["gene_summary.txt", "sgrna_summary.txt"] {
        let path = target(prefix, suffix);
        if let Err(cleanup) = fs::remove_file(&path) {
            failed = true;
            message.push_str(&format!(
                "; cleanup failed for {}: {cleanup}",
                path.display()
            ));
        }
    }
    if failed {
        MageckError::Io(io::Error::other(message))
    } else {
        primary
    }
}

/// Single-control comparisons estimate dispersion from the selected control
/// and treatments, retaining control-only means for the actual comparison.
fn comparison_nb(
    counts: &CountMatrix,
    controls: &[usize],
    treatments: &[usize],
) -> Result<super::NbModel, MageckError> {
    if controls.len() != 1 {
        return fit_meanvar(counts, controls);
    }
    let pooled: Vec<usize> = controls.iter().chain(treatments).copied().collect();
    super::NbModel::fit_pooled_dispersion(counts, controls[0], &pooled)
}

fn validate_outputs(config: &TestConfig) -> Result<(), MageckError> {
    crate::config::validate_test_threads(config.threads)?;
    if config.output_prefix.contains(['\t', '\r', '\n', '\0']) {
        return Err(MageckError::Config(
            "Output prefix cannot contain TSV control characters or NUL".into(),
        ));
    }
    if !matches!(config.adjust_method.as_str(), "fdr" | "bh") {
        return Err(MageckError::Config(format!(
            "Unsupported --adjust-method {:?}; test supports fdr or bh",
            config.adjust_method
        )));
    }
    for suffix in ["gene_summary.txt", "sgrna_summary.txt", "run_info.txt"] {
        reject_existing(&target(&config.output_prefix, suffix))?;
    }
    Ok(())
}

pub(crate) fn run(config: &TestConfig) -> Result<(), MageckError> {
    validate_outputs(config)?;
    let raw = ParsedCountTable::from_reader(BufReader::new(File::open(&config.count_table)?))?;
    analyze(config, raw, &config.output_prefix)
}

/// Mechanical neutral-table adapter; shared analysis below retains every statistical operation.
pub(crate) fn run_from_table(
    config: &TestConfig,
    raw: ParsedCountTable,
    metadata_prefix: &str,
) -> Result<(), MageckError> {
    validate_outputs(config)?;
    analyze(config, raw, metadata_prefix)
}

fn analyze(
    config: &TestConfig,
    raw: ParsedCountTable,
    metadata_prefix: &str,
) -> Result<(), MageckError> {
    let mut profile = crate::execution::Profile::new("test");
    let table = CountTable::new(
        raw.sgrna_ids().to_vec(),
        raw.genes().to_vec(),
        raw.sample_labels().to_vec(),
        raw.counts().clone(),
        &config.treatment_id,
        &config.day0_label,
    )?;
    // File syntax skips blank/comment lines; each retained identifier is opaque.
    // Preserve duplicates and unresolved occurrences for P4F8 diagnostics.
    let controls = config
        .control_sgrna
        .as_ref()
        .map(|path| -> Result<Vec<String>, MageckError> {
            Ok(fs::read_to_string(path)?
                .lines()
                .filter(|line| !line.trim().is_empty() && !line.trim_start().starts_with('#'))
                .map(str::to_owned)
                .collect())
        })
        .transpose()?;
    let normalization_controls = controls.as_ref().map(|ids| {
        let ids: HashSet<&str> = ids.iter().map(String::as_str).collect();
        raw.sgrna_ids()
            .iter()
            .enumerate()
            .filter_map(|(i, id)| ids.contains(id.as_str()).then_some(i))
            .collect::<Vec<_>>()
    });
    let normalized = normalize(
        raw.counts(),
        config.norm_method,
        normalization_controls.as_deref(),
    )?;
    if config.norm_method == NormMethod::Control
        && normalization_controls.as_ref().map_or(true, Vec::is_empty)
    {
        return Err(MageckError::Config(
            "Control normalization requires a matched control guide".into(),
        ));
    }
    let validity = build_comparison_validity(
        &normalized,
        &table,
        config.remove_zero,
        config.remove_zero_threshold,
        "test",
    )?;
    profile.stage("prepare-normalize-mask");
    let nb = comparison_nb(
        &normalized,
        table.control_indices(),
        table.treatment_indices(),
    )?;
    profile.stage("nb");
    let guide_pvalues = calculate_sgrna_pvalues(&normalized, table.treatment_indices(), &nb)?;
    profile.stage("guide-pvalues");
    let neg_p: Vec<f64> = guide_pvalues.iter().map(|g| g.p_lower).collect();
    let pos_p: Vec<f64> = guide_pvalues.iter().map(|g| g.p_higher).collect();
    let neg_ranks = stats::assign_masked_rank_fractions(&neg_p, validity.mask())?;
    let pos_ranks = stats::assign_masked_rank_fractions(&pos_p, validity.mask())?;
    let alpha = stats::calc_directional_alpha(
        &neg_p,
        &pos_p,
        validity.mask(),
        config.gene_test_fdr_threshold,
    )
    .map_err(|e| MageckError::Math(format!("Directional alpha: {e:?}")))?;
    let observed = compute_gene_lo_values(
        raw.genes(),
        &neg_ranks,
        &pos_ranks,
        alpha,
        validity.n_valid(),
    )?;
    profile.stage("rank-alpha-observed");
    let pools = controls
        .as_ref()
        .map(|ids| {
            let ids: Vec<&str> = ids.iter().map(String::as_str).collect();
            prepare_control_pools(
                &ids,
                raw.sgrna_ids(),
                validity.mask(),
                &neg_ranks,
                &pos_ranks,
            )
        })
        .transpose()?;
    let null = if let Some(pools) = &pools {
        log::info!(
            "RRA controls: {:?}; unresolved identifiers: {:?}",
            pools.diagnostics,
            pools.unresolved_identifiers
        );
        NullModelConfig::Control {
            neg_control_percentiles: &pools.neg_percentiles,
            pos_control_percentiles: &pools.pos_percentiles,
            diagnostics: &pools.diagnostics,
        }
    } else {
        NullModelConfig::Ordinary
    };
    profile.stage("control-pools");
    let max_guides = observed
        .iter()
        .map(|g| g.num_valid_sgrnas)
        .max()
        .unwrap_or(1);
    let batch_jobs = observed
        .iter()
        .filter(|g| g.num_valid_sgrnas > 0)
        .map(|g| (65536 / g.num_valid_sgrnas).clamp(1, 64))
        .max()
        .unwrap_or(1);
    let jobs = if config.threads.is_some_and(crate::cli::ThreadMode::is_auto) {
        config.permutation.min(batch_jobs).max(observed.len())
    } else {
        config.permutation.max(observed.len())
    };
    let execution = crate::execution::Execution::stage(
        config.threads,
        jobs,
        64,
        max_guides
            .saturating_mul(32)
            .saturating_add(8 * 1024 * 1024),
        "rra.empirical",
    )?;
    let empirical = compute_empirical_with_execution(
        &observed,
        &null,
        config.permutation,
        config.rra_seed,
        validity.n_valid(),
        alpha,
        &execution,
    )?;
    profile.stage("empirical");
    let fdr = apply_directional_bh(&empirical)?;
    let lfc = aggregate_lfc(&raw, table.control_indices(), table.treatment_indices())?;
    let genes = aggregate_results(&fdr, &lfc)?;
    profile.stage("bh-lfc-aggregate");
    write_rra_summaries(
        Path::new(&config.output_prefix),
        &genes,
        &raw,
        table.control_indices(),
        table.treatment_indices(),
        &neg_ranks,
    )?;
    let metadata = run_info(
        &TestConfig {
            output_prefix: metadata_prefix.into(),
            ..config.clone()
        },
        genes.len(),
        raw.n_guides(),
    ) + &execution.metadata();
    if let Err(primary) = write_metadata(&target(&config.output_prefix, "run_info.txt"), &metadata)
    {
        return Err(rollback_summaries(&config.output_prefix, primary));
    }
    profile.stage("serialize-publish");
    profile.finish();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dropout_fixture() -> CountMatrix {
        // Four guides in one gene, with already normalized control/treatments.
        CountMatrix::from_vec(
            vec![100.0; 4]
                .into_iter()
                .chain([5.0; 4])
                .chain([3.0; 4])
                .chain([4.0; 4])
                .collect(),
            4,
            4,
        )
        .unwrap()
    }

    #[test]
    fn p4f22_single_control_pools_four_columns_and_keeps_control_mean() {
        let counts = dropout_fixture();
        let nb = comparison_nb(&counts, &[0], &[1, 2, 3]).unwrap();
        let pooled = fit_meanvar(&counts, &[0, 1, 2, 3]).unwrap();
        assert_eq!(nb.mean(), &[100.0; 4]);
        assert_eq!(pooled.mean(), &[28.0; 4]);
        // These moments identify all four columns, rather than duplicated
        // control observations or treatment-only fitting.
        assert_eq!(nb.variance(), pooled.variance());
        assert!((nb.variance()[0] - 2304.6666666666665).abs() < 1e-10);
        assert_eq!(nb.raw_dispersion(), pooled.raw_dispersion());
        assert_eq!(nb.fitted_dispersion(), pooled.fitted_dispersion());
        assert_eq!(nb.final_dispersion(), pooled.final_dispersion());
        assert_eq!(nb.size(), pooled.size());
        assert_eq!(nb.eligible_guides_count(), 4);
        assert!(nb.final_dispersion().iter().all(|&d| d > 0.1));
        assert!((nb.final_dispersion()[0] - 2.90391156462585).abs() < 1e-12);
        println!(
            "F_NB1 columns=4; F_NB2 dispersion={:.17}; control_mean={:.17}; pooled_variance={:.17}",
            nb.final_dispersion()[0],
            nb.mean()[0],
            nb.variance()[0]
        );
    }

    #[test]
    fn p4f22_single_control_dropout_pvalue_uses_control_mean() {
        let counts = dropout_fixture();
        let nb = comparison_nb(&counts, &[0], &[1, 2, 3]).unwrap();
        let values = calculate_sgrna_pvalues(&counts, &[1, 2, 3], &nb).unwrap();
        for value in &values {
            assert_eq!(value.control_mean, 100.0);
            assert_eq!(value.treatment_mean, 4.0);
            assert!(value.p_lower > 0.0 && value.p_lower < 0.05);
            assert!((value.p_lower - 0.010995619140642118).abs() < 1e-10);
            assert!((value.lfc - (5.0_f64 / 101.0).log2()).abs() < 1e-12);
        }
        println!("F_NB3 p_lower={:.17}", values[0].p_lower);
    }

    #[test]
    fn p4f22_multi_control_model_and_pvalues_match_baseline() {
        let counts = CountMatrix::from_vec(
            vec![
                100.0, 20.0, 30.0, 5.0, 80.0, 24.0, 60.0, 8.0, 4.0, 22.0, 45.0, 7.0,
            ],
            3,
            4,
        )
        .unwrap();
        let baseline = fit_meanvar(&counts, &[1, 0]).unwrap();
        let nb = comparison_nb(&counts, &[1, 0], &[2]).unwrap();
        assert_eq!(nb, baseline);
        assert_eq!(
            calculate_sgrna_pvalues(&counts, &[2], &nb).unwrap(),
            calculate_sgrna_pvalues(&counts, &[2], &baseline).unwrap()
        );
    }

    #[test]
    fn p4f22_pool_excludes_unselected_samples_and_handles_nonfirst_control() {
        let counts = dropout_fixture();
        let selected = comparison_nb(&counts, &[0], &[1, 2, 3]).unwrap();
        let reordered = CountMatrix::from_vec(
            vec![5.0; 4]
                .into_iter()
                .chain([f64::NAN; 4])
                .chain([100.0; 4])
                .chain([3.0; 4])
                .chain([4.0; 4])
                .collect(),
            5,
            4,
        )
        .unwrap();
        let nb = comparison_nb(&reordered, &[2], &[0, 3, 4]).unwrap();
        assert_eq!(nb, selected);
    }

    #[test]
    fn p4f22_pool_rejects_missing_duplicate_and_out_of_bounds_samples() {
        let counts = dropout_fixture();
        for (control, samples) in [(0, vec![1, 2, 3]), (0, vec![0, 0, 1]), (4, vec![4, 1])] {
            assert!(matches!(
                super::super::NbModel::fit_pooled_dispersion(&counts, control, &samples),
                Err(MageckError::Math(_))
            ));
        }
        assert!(comparison_nb(&counts, &[0], &[]).is_err());
    }

    #[test]
    fn single_control_nb_keeps_means_and_pools_variance() {
        let counts = CountMatrix::from_vec(vec![2.0, 10.0, 1.0, 20.0], 2, 2).unwrap();
        let nb = comparison_nb(&counts, &[0], &[1]).unwrap();
        assert_eq!(nb.mean(), &[2.0, 10.0]);
        let pooled = fit_meanvar(&counts, &[0, 1]).unwrap();
        assert_eq!(nb.variance(), pooled.variance());
        assert_eq!(nb.final_dispersion(), pooled.final_dispersion());
        let multi = comparison_nb(&counts, &[0, 1], &[]).unwrap();
        assert_eq!(multi, fit_meanvar(&counts, &[0, 1]).unwrap());
    }

    #[test]
    fn metadata_atomic_write_and_existing_output_rejection() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run_info.txt");
        write_metadata(&path, "key\tvalue\n").unwrap();
        assert_eq!(fs::read_to_string(&path).unwrap(), "key\tvalue\n");
        assert!(write_metadata(&path, "replace").is_err());
        assert_eq!(fs::read_to_string(&path).unwrap(), "key\tvalue\n");
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
    }

    #[test]
    fn exclusive_metadata_rename_failure_cleans_temporary() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out");
        fs::create_dir(&path).unwrap();
        assert!(write_metadata(&path, "new").is_err());
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
    }

    #[test]
    fn metadata_failure_rolls_back_owned_summaries() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("out").to_string_lossy().into_owned();
        for suffix in ["gene_summary.txt", "sgrna_summary.txt"] {
            fs::write(target(&prefix, suffix), "new").unwrap();
        }
        let error = rollback_summaries(&prefix, MageckError::Config("metadata error".into()));
        assert!(error.to_string().contains("metadata error"));
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
    }
}
