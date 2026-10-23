use crate::cli::{NormMethod, Trim5Mode};
use crate::config::CountConfig;
use crate::error::MageckError;
use std::path::Path;

use crate::count::control::parse_control_sgrna_file;
use crate::count::matcher::{build_index, parse_library_file};
use crate::count::normalize::normalize;
use crate::count::reader::count_grouped_fastq_with_trim;
use crate::count::writer::{write_count_outputs, CountOutputPaths};

pub(crate) fn format_unknown_control_warning(path: &Path, count: usize) -> String {
    format!(
        "Control file '{}' contained {} unknown entries.",
        path.display(),
        count
    )
}

pub(crate) fn format_unused_control_warning(path: &Path, method: NormMethod) -> String {
    format!(
        "Control file '{}' was supplied but will not influence '{}' normalization.",
        path.display(),
        method
    )
}

/// Orchestrates the `count` command pipeline from configuration to output.
pub(crate) struct PreparedCount {
    pub library: crate::count::Library,
    pub grouped: crate::count::GroupedCountResult,
    normalized: crate::data::CountMatrix,
    umi_summaries: Vec<crate::count::reader::UmiFileSummary>,
    thread_metadata: String,
}

pub(crate) fn prepare(
    config: &CountConfig,
    unmapped: Option<&mut dyn std::io::Write>,
) -> Result<PreparedCount, MageckError> {
    crate::config::validate_count_options(config.count_n, config.unmapped_to_file, config.threads)?;
    // 1. Initial sanity check for Control normalization
    if matches!(config.norm_method, NormMethod::Control) && config.control_sgrna.is_none() {
        return Err(MageckError::Config(
            "Control normalization requires a control-sgrna file".to_string(),
        ));
    }

    // 2. Parse Library
    let library = parse_library_file(&config.list_seq)?;

    // 3. Build Guide Index
    let guide_index = build_index(
        &library,
        config.test_mismatch == 1 && config.trim_5 != Trim5Mode::Auto,
    )?;

    // 4. Resolve Control Guides if supplied
    let resolved_controls = if let Some(control_path) = &config.control_sgrna {
        let resolved = parse_control_sgrna_file(control_path, library.records())?;

        if !resolved.unknown_entries.is_empty() {
            log::warn!(
                "{}",
                format_unknown_control_warning(control_path, resolved.unknown_entries.len())
            );
        }

        if !matches!(config.norm_method, NormMethod::Control) {
            log::warn!(
                "{}",
                format_unused_control_warning(control_path, config.norm_method)
            );
            None
        } else {
            Some(resolved.guide_indexes)
        }
    } else {
        None
    };

    let offsets = match config.trim_5 {
        Trim5Mode::None => "0".to_string(),
        Trim5Mode::Fixed(offset) => offset.to_string(),
        Trim5Mode::Auto => "0..=40".to_string(),
    };
    log::info!(
        "Count guide start offsets: {offsets} (trim-5 {:?})",
        config.trim_5
    );

    let execution =
        crate::execution::Execution::stage(config.threads, 1, 1, 8 * 1024 * 1024, "count.serial")?;
    let thread_metadata = execution.metadata();
    // 5. Count Grouped FASTQ; preserve the original reader when no option is enabled.
    let (grouped_result, umi_summaries) = if config.umi.is_some() || config.unmapped_to_file {
        crate::count::reader::count_grouped_with_options(config, &guide_index, unmapped)?
    } else {
        (
            count_grouped_fastq_with_trim(
                &config.sample_groups,
                &guide_index,
                config.trim_5,
                config.sgrna_len,
                config.test_mismatch == 1,
            )?,
            Vec::new(),
        )
    };

    // 6. Normalize
    let normalized_matrix = normalize(
        &grouped_result.counts,
        config.norm_method,
        resolved_controls.as_deref(),
    )?;

    Ok(PreparedCount {
        library,
        grouped: grouped_result,
        normalized: normalized_matrix,
        umi_summaries,
        thread_metadata,
    })
}

pub(crate) fn write_prepared(
    config: &CountConfig,
    prepared: &PreparedCount,
) -> Result<CountOutputPaths, MageckError> {
    let grouped_result = &prepared.grouped;
    let library = &prepared.library;
    let normalized_matrix = &prepared.normalized;
    let offsets = match config.trim_5 {
        Trim5Mode::None => "0".to_string(),
        Trim5Mode::Fixed(n) => n.to_string(),
        Trim5Mode::Auto => "0..=40".to_string(),
    };
    // 7. Write Outputs
    let sample_labels: Vec<&str> = config
        .sample_groups
        .iter()
        .map(|g| g.label.as_str())
        .collect();

    let outputs = write_count_outputs(
        &config.output_prefix,
        &grouped_result.counts,
        normalized_matrix,
        library.records(),
        &sample_labels,
        &grouped_result.summaries,
    )?;
    std::fs::write(
        format!("{}.run_info.txt", config.output_prefix),
        format!(
            "key\tvalue\ntrim_5_mode\t{:?}\ntrim_5_offsets\t{offsets}\n",
            config.trim_5
        ) + &prepared.thread_metadata,
    )?;
    if let Some(umi) = &config.umi {
        use std::io::Write;
        let mut writer = std::io::BufWriter::new(std::fs::File::create(format!(
            "{}.umi_summary.txt",
            config.output_prefix
        ))?);
        writeln!(
            writer,
            "File\tLabel\tValidMappedReads\tNewMolecules\tDuplicateMolecules"
        )?;
        for row in &prepared.umi_summaries {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                row.path.display(),
                row.label,
                row.valid_mapped_reads,
                row.new_molecules,
                row.duplicate_molecules
            )?;
        }
        writer.flush()?;
        let mut metadata = std::fs::OpenOptions::new()
            .append(true)
            .open(format!("{}.run_info.txt", config.output_prefix))?;
        writeln!(metadata, "umi_start\t{}\numi_length\t{}\numi_dedup_scope\t{}\numi_max_keys\t{}\numi_equality\tuppercase_ACGT_exact\numi_count_semantics\tmolecules\numi_qc_semantics\tphysical_reads_before_deduplication\numi_cap_is_process_memory_bound\tfalse", umi.start, umi.length,
            match umi.scope { crate::cli::UmiDedupScope::File => "file", crate::cli::UmiDedupScope::Sample => "sample" }, umi.max_keys)?;
    }
    Ok(outputs)
}

pub(crate) fn output_suffixes(config: &CountConfig) -> Vec<&'static str> {
    let mut suffixes = vec![
        "count.txt",
        "count.normalized.txt",
        "countsummary.txt",
        "run_info.txt",
    ];
    if config.umi.is_some() {
        suffixes.push("umi_summary.txt");
    }
    if config.unmapped_to_file {
        suffixes.push("unmapped.fastq");
    }
    suffixes
}

/// Ordinary COUNT retains its existing publication contract. Opt-in workflows stage all outputs.
pub fn run(config: &CountConfig) -> Result<CountOutputPaths, MageckError> {
    use std::io::Write;
    crate::config::validate_count_options(config.count_n, config.unmapped_to_file, config.threads)?;
    if let Some(umi) = &config.umi {
        umi.validate()?;
    }
    if config.umi.is_none() && !config.unmapped_to_file {
        let prepared = prepare(config, None)?;
        return write_prepared(config, &prepared);
    }
    let suffixes = output_suffixes(config);
    let stage = crate::run::OutputStage::new(&config.output_prefix, &suffixes)?;
    let mut staged = config.clone();
    staged.output_prefix = stage.prefix.clone();
    let mut writer = if config.unmapped_to_file {
        Some(std::io::BufWriter::new(std::fs::File::create(format!(
            "{}.unmapped.fastq",
            stage.prefix
        ))?))
    } else {
        None
    };
    let prepared = prepare(&staged, writer.as_mut().map(|w| w as &mut dyn Write))?;
    if let Some(w) = writer.as_mut() {
        w.flush()?;
    }
    write_prepared(&staged, &prepared)?;
    stage.publish(&config.output_prefix, &suffixes)?;
    Ok(CountOutputPaths {
        count_path: format!("{}.count.txt", config.output_prefix).into(),
        normalized_path: format!("{}.count.normalized.txt", config.output_prefix).into(),
        summary_path: format!("{}.countsummary.txt", config.output_prefix).into(),
    })
}
