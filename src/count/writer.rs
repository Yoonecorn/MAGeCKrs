//! Pure-data TSV output serialization for count pipeline outputs.
//!
//! Implements Phase 3D deliverables:
//! - `{prefix}.count.txt`: Raw integer count table.
//! - `{prefix}.count.normalized.txt`: Normalized count table formatted to 4 decimal places (`{:.4}`).
//! - `{prefix}.countsummary.txt`: Per-sample QC mapping statistics and zero-count guides.
//!
//! All files are serialized as UTF-8, tab-separated values (`\t`), terminated with Unix newlines (`\n`).

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::count::matcher::SgRnaRecord;
use crate::count::reader::PhysicalFastqSummary;
use crate::data::count_matrix::CountMatrix;
use crate::error::MageckError;

/// Paths to the three generated count output files.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CountOutputPaths {
    /// Path to `{prefix}.count.txt`.
    pub count_path: PathBuf,
    /// Path to `{prefix}.count.normalized.txt`.
    pub normalized_path: PathBuf,
    /// Path to `{prefix}.countsummary.txt`.
    pub summary_path: PathBuf,
}

const MAX_SAFE_INT_F64: f64 = 9_007_199_254_740_992.0; // 2^53

/// Validates that a raw count value is a finite, non-negative integer that fits exactly in f64.
#[inline]
fn validate_raw_count(val: f64) -> Result<u64, MageckError> {
    if !val.is_finite() || val < 0.0 {
        return Err(MageckError::Config(format!(
            "Raw count matrix contains invalid count {val}: must be non-negative and finite"
        )));
    }
    if val.fract() != 0.0 {
        return Err(MageckError::Config(format!(
            "Raw count matrix contains non-integral count {val}: counts must represent integers"
        )));
    }
    if val > MAX_SAFE_INT_F64 {
        return Err(MageckError::Config(format!(
            "Raw count matrix value {val} exceeds exact f64 integer representation (2^53)"
        )));
    }
    Ok(val as u64)
}

/// Validates that a string does not contain TSV control characters (tabs, newlines).
#[inline]
fn validate_tsv_field(field: &str, field_name: &str) -> Result<(), MageckError> {
    if field.contains('\t') || field.contains('\n') || field.contains('\r') {
        return Err(MageckError::Config(format!(
            "TSV field '{}' ({}) contains invalid control characters",
            field_name, field
        )));
    }
    Ok(())
}

/// Validates that a normalized count value is a finite, non-negative number.
#[inline]
fn validate_normalized_count(val: f64) -> Result<(), MageckError> {
    if !val.is_finite() || val < 0.0 {
        return Err(MageckError::Config(format!(
            "Normalized count matrix contains invalid count {val}: must be non-negative and finite"
        )));
    }
    Ok(())
}

/// Serializes raw guide count matrix to a tab-separated file (`{prefix}.count.txt`).
///
/// # Output Contract
/// - Header: `sgRNA\tGene\t{label1}\t{label2}...\n`
/// - Rows: One row per library record in supplied library order.
/// - Values: Non-negative integers serialized without padding (e.g. `0`, `1`, `42`).
/// - Rejects negative, NaN, infinite, or non-integral values with [`MageckError::Config`].
/// - Validates dimensional constraints before creating the file on disk.
///
/// # Default Behavior
/// Default file creation overwrites an existing file if present.
pub fn write_raw_counts<P: AsRef<Path>, S: AsRef<str>>(
    path: P,
    matrix: &CountMatrix,
    records: &[SgRnaRecord],
    sample_labels: &[S],
) -> Result<(), MageckError> {
    let n_samples = matrix.n_samples();
    let n_guides = matrix.n_guides();

    if sample_labels.len() != n_samples {
        return Err(MageckError::Config(format!(
            "Sample label count ({}) does not match matrix sample count ({n_samples})",
            sample_labels.len()
        )));
    }
    if records.len() != n_guides {
        return Err(MageckError::Config(format!(
            "Library record count ({}) does not match matrix guide count ({n_guides})",
            records.len()
        )));
    }

    // Validate all values before creating destination file to avoid partial output
    for &val in matrix.as_slice() {
        validate_raw_count(val)?;
    }
    for label in sample_labels {
        validate_tsv_field(label.as_ref(), "sample label")?;
    }
    for record in records {
        validate_tsv_field(&record.id, "sgRNA ID")?;
        validate_tsv_field(&record.gene, "Gene")?;
    }

    let file = File::create(path.as_ref()).map_err(MageckError::Io)?;
    let mut writer = BufWriter::new(file);

    // Write header: sgRNA\tGene\t{label1}\t{label2}...\n
    write!(writer, "sgRNA\tGene")?;
    for label in sample_labels {
        write!(writer, "\t{}", label.as_ref())?;
    }
    writeln!(writer)?;

    // Write rows in library guide order: {id}\t{gene}\t{count_s0}\t{count_s1}...\n
    for (g, record) in records.iter().enumerate() {
        write!(writer, "{}\t{}", record.id, record.gene)?;
        for s in 0..n_samples {
            let val = matrix.entry(s, g);
            let val_u64 = validate_raw_count(val)?;
            write!(writer, "\t{val_u64}")?;
        }
        writeln!(writer)?;
    }

    writer.flush().map_err(MageckError::Io)?;
    Ok(())
}

/// Serializes normalized guide count matrix to a tab-separated file (`{prefix}.count.normalized.txt`).
///
/// # Output Contract
/// - Header: `sgRNA\tGene\t{label1}\t{label2}...\n`
/// - Rows: One row per library record in supplied library order.
/// - Values: Formatted with exactly four decimal places (`{:.4}`), preserving full precision in matrix.
/// - Rejects negative, NaN, or infinite values with [`MageckError::Config`].
/// - Validates dimensional constraints before creating the file on disk.
///
/// # Default Behavior
/// Default file creation overwrites an existing file if present.
pub fn write_normalized_counts<P: AsRef<Path>, S: AsRef<str>>(
    path: P,
    matrix: &CountMatrix,
    records: &[SgRnaRecord],
    sample_labels: &[S],
) -> Result<(), MageckError> {
    let n_samples = matrix.n_samples();
    let n_guides = matrix.n_guides();

    if sample_labels.len() != n_samples {
        return Err(MageckError::Config(format!(
            "Sample label count ({}) does not match matrix sample count ({n_samples})",
            sample_labels.len()
        )));
    }
    if records.len() != n_guides {
        return Err(MageckError::Config(format!(
            "Library record count ({}) does not match matrix guide count ({n_guides})",
            records.len()
        )));
    }

    // Validate all values before creating destination file to avoid partial output
    for &val in matrix.as_slice() {
        validate_normalized_count(val)?;
    }
    for label in sample_labels {
        validate_tsv_field(label.as_ref(), "sample label")?;
    }
    for record in records {
        validate_tsv_field(&record.id, "sgRNA ID")?;
        validate_tsv_field(&record.gene, "Gene")?;
    }

    let file = File::create(path.as_ref()).map_err(MageckError::Io)?;
    let mut writer = BufWriter::new(file);

    // Write header: sgRNA\tGene\t{label1}\t{label2}...\n
    write!(writer, "sgRNA\tGene")?;
    for label in sample_labels {
        write!(writer, "\t{}", label.as_ref())?;
    }
    writeln!(writer)?;

    // Write rows in library guide order: {id}\t{gene}\t{norm_s0:.4}\t{norm_s1:.4}...\n
    for (g, record) in records.iter().enumerate() {
        write!(writer, "{}\t{}", record.id, record.gene)?;
        for s in 0..n_samples {
            let val = matrix.entry(s, g);
            validate_normalized_count(val)?;
            write!(writer, "\t{val:.4}")?;
        }
        writeln!(writer)?;
    }

    writer.flush().map_err(MageckError::Io)?;
    Ok(())
}

/// Serializes per-sample alignment and quality summary to `{prefix}.countsummary.txt`.
///
/// # Output Contract
/// - Exact Header: `File\tLabel\tReads\tMapped\tPercentage\tZerocounts\n`
/// - Rows: One row per sample/FASTQ, preserving supplied input order.
/// - Reads: Total raw reads from `summary.total_reads`.
/// - Mapped: `summary.mapped_reads` (`exact_matches + one_mismatch_matches`).
/// - Percentage: Deterministic implementation decision formatting mapped fraction `summary.mapped_percentage()`
///   in `[0.0, 1.0]` to four decimal places (`{:.4}`). If `Reads == 0`, evaluates to `0.0000` (finite, represents zero).
/// - Zerocounts: Number of guide rows for that sample where raw count is exactly `0.0`.
///
/// # Default Behavior
/// Default file creation overwrites an existing file if present.
pub fn write_count_summary<P: AsRef<Path>>(
    path: P,
    summaries: &[PhysicalFastqSummary],
) -> Result<(), MageckError> {
    // Validate TSV safety and invariants
    for (s, phys_summary) in summaries.iter().enumerate() {
        validate_tsv_field(&phys_summary.label, "sample label")?;

        let p_str = phys_summary.path.to_str().ok_or_else(|| {
            MageckError::Config(format!(
                "FASTQ path '{}' is not valid UTF-8",
                phys_summary.path.to_string_lossy()
            ))
        })?;
        validate_tsv_field(p_str, "FASTQ path")?;

        let summary = &phys_summary.summary;
        let sum_parts = summary.mapped_reads
            + summary.ambiguous_reads
            + summary.unmapped_reads
            + summary.short_reads;
        if summary.total_reads != sum_parts {
            return Err(MageckError::Config(format!(
                "Sample {s} summary invariant violated: total_reads ({}) != mapped ({}) + ambiguous ({}) + unmapped ({}) + short ({})",
                summary.total_reads, summary.mapped_reads, summary.ambiguous_reads, summary.unmapped_reads, summary.short_reads
            )));
        }
        if summary.mapped_reads != summary.exact_matches + summary.one_mismatch_matches {
            return Err(MageckError::Config(format!(
                "Sample {s} summary invariant violated: mapped_reads ({}) != exact_matches ({}) + one_mismatch_matches ({})",
                summary.mapped_reads, summary.exact_matches, summary.one_mismatch_matches
            )));
        }
    }

    let file = File::create(path.as_ref()).map_err(MageckError::Io)?;
    let mut writer = BufWriter::new(file);

    // Exact MAGeCK2 header
    writeln!(writer, "File\tLabel\tReads\tMapped\tPercentage\tZerocounts")?;

    for phys_summary in summaries {
        let fq_display = phys_summary.path.to_str().unwrap(); // validated above
        let label = &phys_summary.label;
        let summary = &phys_summary.summary;
        let zerocounts = phys_summary.zerocounts;

        // Deterministic Percentage formatting: mapped fraction in [0.0, 1.0] formatted with {:.4}
        let pct = summary.mapped_percentage();

        writeln!(
            writer,
            "{fq_display}\t{label}\t{}\t{}\t{pct:.4}\t{zerocounts}",
            summary.total_reads, summary.mapped_reads
        )?;
    }

    writer.flush().map_err(MageckError::Io)?;
    Ok(())
}

/// Convenience function that validates dimensions across all outputs and writes
/// `{prefix}.count.txt`, `{prefix}.count.normalized.txt`, and `{prefix}.countsummary.txt`.
///
/// # Validation
/// Enforces all dimensional and numerical invariants prior to creating any file on disk.
pub fn write_count_outputs<P: AsRef<Path>, S: AsRef<str>>(
    output_prefix: P,
    raw_matrix: &CountMatrix,
    normalized_matrix: &CountMatrix,
    records: &[SgRnaRecord],
    sample_labels: &[S],
    summaries: &[PhysicalFastqSummary],
) -> Result<CountOutputPaths, MageckError> {
    // 1. Validate raw vs normalized dimensions
    if raw_matrix.n_samples() != normalized_matrix.n_samples()
        || raw_matrix.n_guides() != normalized_matrix.n_guides()
    {
        return Err(MageckError::Config(format!(
            "Raw matrix dimensions [{}, {}] do not match normalized matrix dimensions [{}, {}]",
            raw_matrix.n_samples(),
            raw_matrix.n_guides(),
            normalized_matrix.n_samples(),
            normalized_matrix.n_guides()
        )));
    }

    let n_samples = raw_matrix.n_samples();
    let n_guides = raw_matrix.n_guides();

    // 2. Validate metadata lengths against matrix dimensions
    if sample_labels.len() != n_samples {
        return Err(MageckError::Config(format!(
            "Sample label count ({}) does not match matrix sample count ({n_samples})",
            sample_labels.len()
        )));
    }
    if records.len() != n_guides {
        return Err(MageckError::Config(format!(
            "Library record count ({}) does not match matrix guide count ({n_guides})",
            records.len()
        )));
    }

    // 3. Validate numerical values before file creation
    for &val in raw_matrix.as_slice() {
        validate_raw_count(val)?;
    }
    for &val in normalized_matrix.as_slice() {
        validate_normalized_count(val)?;
    }

    // 4. Validate summary invariants and TSV strings for summaries
    for (s, phys_summary) in summaries.iter().enumerate() {
        validate_tsv_field(&phys_summary.label, "sample label")?;

        let p_str = phys_summary.path.to_str().ok_or_else(|| {
            MageckError::Config(format!(
                "FASTQ path '{}' is not valid UTF-8",
                phys_summary.path.to_string_lossy()
            ))
        })?;
        validate_tsv_field(p_str, "FASTQ path")?;

        let summary = &phys_summary.summary;
        let sum_parts = summary.mapped_reads
            + summary.ambiguous_reads
            + summary.unmapped_reads
            + summary.short_reads;
        if summary.total_reads != sum_parts {
            return Err(MageckError::Config(format!(
                "Sample {s} summary invariant violated: total_reads ({}) != mapped ({}) + ambiguous ({}) + unmapped ({}) + short ({})",
                summary.total_reads, summary.mapped_reads, summary.ambiguous_reads, summary.unmapped_reads, summary.short_reads
            )));
        }
        if summary.mapped_reads != summary.exact_matches + summary.one_mismatch_matches {
            return Err(MageckError::Config(format!(
                "Sample {s} summary invariant violated: mapped_reads ({}) != exact_matches ({}) + one_mismatch_matches ({})",
                summary.mapped_reads, summary.exact_matches, summary.one_mismatch_matches
            )));
        }
    }

    // 5. Validate remaining TSV strings
    for label in sample_labels {
        validate_tsv_field(label.as_ref(), "sample label")?;
    }
    for record in records {
        validate_tsv_field(&record.id, "sgRNA ID")?;
        validate_tsv_field(&record.gene, "Gene")?;
    }

    // Construct destination file paths
    let mut count_os = output_prefix.as_ref().as_os_str().to_os_string();
    count_os.push(".count.txt");
    let count_path = PathBuf::from(count_os);

    let mut norm_os = output_prefix.as_ref().as_os_str().to_os_string();
    norm_os.push(".count.normalized.txt");
    let normalized_path = PathBuf::from(norm_os);

    let mut summary_os = output_prefix.as_ref().as_os_str().to_os_string();
    summary_os.push(".countsummary.txt");
    let summary_path = PathBuf::from(summary_os);

    // 6. Serialize all files
    write_raw_counts(&count_path, raw_matrix, records, sample_labels)?;
    write_normalized_counts(&normalized_path, normalized_matrix, records, sample_labels)?;
    write_count_summary(&summary_path, summaries)?;

    Ok(CountOutputPaths {
        count_path,
        normalized_path,
        summary_path,
    })
}
