//! FASTQ record streaming, read-window extraction, and guide-matching.
//!
//! Implements Phase 3B of the count pipeline:
//! - Streaming FASTQ parsing via `needletail` without full-file heap buffers.
//! - Fixed-offset read window extraction `[sgrna_start .. sgrna_start + sgrna_len)`
//!   without unnecessary string or vector allocations.
//! - Guide matching dispatching to [`GuideIndex`] (with exact precedence, optional
//!   1-mismatch matching, and silent skipping of ambiguous reads).
//! - Sequential raw count accumulation in original library guide order.
//! - Classification of every read into mutually exclusive summary counters.
use crate::cli::Trim5Mode;
use crate::error::MageckError;

use std::path::{Path, PathBuf};

use crate::config::FastqSampleGroup;

use needletail::errors::ParseErrorKind;
use needletail::parser::FastxReader;

use crate::count::matcher::GuideIndex;
use crate::data::count_matrix::{CountMatrix, CountMatrixError};
use crate::data::MatchResult;

/// Summary statistics for read processing across a single FASTQ input.
///
/// Every read processed is classified into exactly one mutually exclusive category:
/// - [`exact_matches`](Self::exact_matches) (mapped)
/// - [`one_mismatch_matches`](Self::one_mismatch_matches) (mapped, when enabled)
/// - [`ambiguous_reads`](Self::ambiguous_reads) (unmapped / silently skipped)
/// - [`unmapped_reads`](Self::unmapped_reads) (unmapped / no match in index)
/// - [`short_reads`](Self::short_reads) (unmapped / sequence shorter than window)
///
/// Mathematical invariants:
/// - `mapped_reads == exact_matches + one_mismatch_matches`
/// - `total_reads == mapped_reads + ambiguous_reads + unmapped_reads + short_reads`
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ReadSummary {
    /// Total number of valid FASTQ records processed.
    pub total_reads: usize,
    /// Number of reads successfully mapped to an sgRNA (`exact_matches + one_mismatch_matches`).
    pub mapped_reads: usize,
    /// Number of reads exactly matching an sgRNA in the index.
    pub exact_matches: usize,
    /// Number of reads uniquely matching an sgRNA with 1 mismatch (when enabled).
    pub one_mismatch_matches: usize,
    /// Number of reads matching multiple sgRNAs with 1 mismatch (silently skipped).
    pub ambiguous_reads: usize,
    /// Number of reads with no matching sgRNA in the index.
    pub unmapped_reads: usize,
    /// Number of reads shorter than the extraction window (`sgrna_start + sgrna_len`).
    pub short_reads: usize,
}

impl ReadSummary {
    /// Returns the mapped fraction of total reads in `[0.0, 1.0]`.
    ///
    /// Despite the legacy method name `mapped_percentage`, MAGeCK's `Percentage`
    /// count-summary column stores a mapped fraction on the scale `[0.0, 1.0]`,
    /// not a conventional percentage in `[0.0, 100.0]`.
    ///
    /// Returns `0.0` if `total_reads` is 0. For any valid `ReadSummary` where
    /// `mapped_reads <= total_reads`, the result is finite and in `[0.0, 1.0]`.
    #[inline]
    pub fn mapped_percentage(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            self.mapped_reads as f64 / self.total_reads as f64
        }
    }
}

const MAX_SAFE_INT_F64: f64 = 9_007_199_254_740_992.0; // 2^53

/// Summary record for a single physical FASTQ file within a sample group.
///
/// Retains the supplied physical file path, the sample group label, the independent
/// read summary statistics, and the unaggregated zero-count metric for that physical file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhysicalFastqSummary {
    /// Supplied path to the physical FASTQ file.
    pub path: PathBuf,
    /// Shared sample label for the sample group containing this file.
    pub label: String,
    /// Read classification summary statistics for this physical file.
    pub summary: ReadSummary,
    /// Number of library guides with zero counts in this physical file.
    pub zerocounts: usize,
}

/// Result of grouped sequential FASTQ processing.
///
/// Contains the aggregated raw count matrix (one row per sample group) and the
/// ordered per-physical-file summary records across all sample groups.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupedCountResult {
    /// Raw count matrix where each row corresponds to one sample group.
    pub counts: CountMatrix,
    /// Ordered per-physical-file summaries across all sample groups.
    pub summaries: Vec<PhysicalFastqSummary>,
}

/// Extracts a borrowed subslice `&read_seq[start .. start + len]` representing the guide window.
///
/// # Invariants
/// - Uses 0-based indexing for `start`.
/// - Uses standard half-open range `[start, start + len)`.
/// - Returns `None` if `start + len` overflows `usize` or exceeds `read_seq.len()` (short read).
/// - Does not allocate memory or perform UTF-8 validation.
#[inline]
pub(crate) fn extract_window(read_seq: &[u8], start: usize, len: usize) -> Option<&[u8]> {
    let end = start.checked_add(len)?;
    if read_seq.len() < end {
        None
    } else {
        Some(&read_seq[start..end])
    }
}

/// First exact match in ascending offset order; each read contributes at most once.
fn match_auto_window(seq: &[u8], index: &GuideIndex, sgrna_len: usize) -> MatchResult {
    for offset in 0..=40 {
        let Some(window) = extract_window(seq, offset, sgrna_len) else {
            break;
        };
        if let result @ MatchResult::Exact(_) = index.lookup(window) {
            return result;
        }
    }
    MatchResult::NoMatch
}

/// Core streaming record processor over an active [`FastxReader`].
fn process_fastx_records<R: FastxReader + ?Sized>(
    reader: &mut R,
    source_name: &str,
    index: &GuideIndex,
    trim_5: Trim5Mode,
    sgrna_len: usize,
    one_mismatch: bool,
) -> Result<(Vec<f64>, ReadSummary), MageckError> {
    let sgrna_start = match trim_5 {
        Trim5Mode::Fixed(offset) => offset,
        Trim5Mode::None | Trim5Mode::Auto => 0,
    };
    if sgrna_len == 0 {
        return Err(MageckError::Fastq(
            "sgRNA length must be greater than 0".to_string(),
        ));
    }
    if sgrna_start.checked_add(sgrna_len).is_none() {
        return Err(MageckError::Fastq(format!(
            "Window range overflow: start ({sgrna_start}) + len ({sgrna_len}) exceeds usize::MAX"
        )));
    }

    let n_guides = index.len();
    let mut counts = vec![0.0; n_guides];
    let mut summary = ReadSummary::default();

    while let Some(record_res) = reader.next() {
        let record = record_res.map_err(|e| {
            MageckError::Fastq(format!(
                "Error parsing FASTQ record in '{source_name}': {e}"
            ))
        })?;

        summary.total_reads = summary
            .total_reads
            .checked_add(1)
            .ok_or_else(|| MageckError::Fastq("total_reads counter overflow".to_string()))?;

        let seq = record.seq();
        let result = if trim_5 == Trim5Mode::Auto {
            match_auto_window(&seq, index, sgrna_len)
        } else {
            let window = match extract_window(&seq, sgrna_start, sgrna_len) {
                Some(window) => window,
                None => {
                    summary.short_reads = summary.short_reads.checked_add(1).ok_or_else(|| {
                        MageckError::Fastq("short_reads counter overflow".to_string())
                    })?;
                    continue;
                }
            };
            index.lookup(window)
        };

        match result {
            MatchResult::Exact(guide_idx) => {
                if let Some(c) = counts.get_mut(guide_idx) {
                    *c += 1.0;
                }
                summary.exact_matches = summary.exact_matches.checked_add(1).ok_or_else(|| {
                    MageckError::Fastq("exact_matches counter overflow".to_string())
                })?;
                summary.mapped_reads = summary.mapped_reads.checked_add(1).ok_or_else(|| {
                    MageckError::Fastq("mapped_reads counter overflow".to_string())
                })?;
            }
            MatchResult::OneMismatch(guide_idx) => {
                if one_mismatch {
                    if let Some(c) = counts.get_mut(guide_idx) {
                        *c += 1.0;
                    }
                    summary.one_mismatch_matches =
                        summary.one_mismatch_matches.checked_add(1).ok_or_else(|| {
                            MageckError::Fastq("one_mismatch_matches counter overflow".to_string())
                        })?;
                    summary.mapped_reads =
                        summary.mapped_reads.checked_add(1).ok_or_else(|| {
                            MageckError::Fastq("mapped_reads counter overflow".to_string())
                        })?;
                } else {
                    summary.unmapped_reads =
                        summary.unmapped_reads.checked_add(1).ok_or_else(|| {
                            MageckError::Fastq("unmapped_reads counter overflow".to_string())
                        })?;
                }
            }
            MatchResult::Ambiguous => {
                summary.ambiguous_reads =
                    summary.ambiguous_reads.checked_add(1).ok_or_else(|| {
                        MageckError::Fastq("ambiguous_reads counter overflow".to_string())
                    })?;
            }
            MatchResult::NoMatch => {
                summary.unmapped_reads =
                    summary.unmapped_reads.checked_add(1).ok_or_else(|| {
                        MageckError::Fastq("unmapped_reads counter overflow".to_string())
                    })?;
            }
        }
    }

    Ok((counts, summary))
}

/// Process a single FASTQ file using streaming I/O.
///
/// Supports both uncompressed (`.fastq`, `.fq`) and gzip-compressed (`.fastq.gz`, `.fq.gz`)
/// inputs through `needletail::parse_fastx_file`.
///
/// If the FASTQ file is empty, returns all-zero counts and a zeroed [`ReadSummary`].
///
/// # Arguments
/// - `path`: Path to FASTQ file.
/// - `index`: Reference guide index.
/// - `sgrna_start`: 0-based offset from start of read to guide sequence (e.g. `trim_5`).
/// - `sgrna_len`: Length of guide sequence to extract (e.g. `sgrna_len`).
/// - `one_mismatch`: Whether 1-mismatch matching is enabled.
///
/// # Returns
/// - Raw guide counts vector `Vec<f64>` indexed in original library order.
/// - [`ReadSummary`] reporting mutually exclusive read classification counts.
pub fn process_one_fastq<P: AsRef<Path>>(
    path: P,
    index: &GuideIndex,
    sgrna_start: usize,
    sgrna_len: usize,
    one_mismatch: bool,
) -> Result<(Vec<f64>, ReadSummary), MageckError> {
    process_one_fastq_with_trim(
        path,
        index,
        Trim5Mode::Fixed(sgrna_start),
        sgrna_len,
        one_mismatch,
    )
}

fn process_one_fastq_with_trim<P: AsRef<Path>>(
    path: P,
    index: &GuideIndex,
    trim_5: Trim5Mode,
    sgrna_len: usize,
    one_mismatch: bool,
) -> Result<(Vec<f64>, ReadSummary), MageckError> {
    let path_ref = path.as_ref();
    let mut reader = match needletail::parse_fastx_file(path_ref) {
        Ok(r) => r,
        Err(e) if e.kind == ParseErrorKind::EmptyFile => {
            return Ok((vec![0.0; index.len()], ReadSummary::default()));
        }
        Err(e) => {
            return Err(MageckError::Fastq(format!(
                "Failed to open FASTQ file '{}': {e}",
                path_ref.display()
            )));
        }
    };

    process_fastx_records(
        &mut *reader,
        &path_ref.display().to_string(),
        index,
        trim_5,
        sgrna_len,
        one_mismatch,
    )
}

/// Process FASTQ records streamed from an arbitrary byte reader.
///
/// Useful for testing with in-memory buffers or piping from non-file streams.
/// If the stream is empty, returns all-zero counts and a zeroed [`ReadSummary`].
#[allow(dead_code)]
pub(crate) fn process_fastq_stream<R: std::io::Read + Send + 'static>(
    reader: R,
    index: &GuideIndex,
    sgrna_start: usize,
    sgrna_len: usize,
    one_mismatch: bool,
) -> Result<(Vec<f64>, ReadSummary), MageckError> {
    let mut fastx_reader = match needletail::parse_fastx_reader(reader) {
        Ok(r) => r,
        Err(e) if e.kind == ParseErrorKind::EmptyFile => {
            return Ok((vec![0.0; index.len()], ReadSummary::default()));
        }
        Err(e) => {
            return Err(MageckError::Fastq(format!(
                "Failed to initialize FASTQ stream reader: {e}"
            )));
        }
    };

    process_fastx_records(
        &mut *fastx_reader,
        "<stream>",
        index,
        Trim5Mode::Fixed(sgrna_start),
        sgrna_len,
        one_mismatch,
    )
}

/// Assembles multiple per-sample count vectors into a contiguous [`CountMatrix`].
///
/// Each slice in `sample_counts` must have length equal to `n_guides`.
pub(crate) fn assemble_count_matrix(
    sample_counts: &[Vec<f64>],
    n_guides: usize,
) -> Result<CountMatrix, CountMatrixError> {
    let n_samples = sample_counts.len();
    for counts in sample_counts {
        if counts.len() != n_guides {
            return Err(CountMatrixError::DimensionMismatch {
                expected: n_guides,
                actual: counts.len(),
            });
        }
    }
    let total_len = n_samples
        .checked_mul(n_guides)
        .ok_or(CountMatrixError::DimensionOverflow {
            n_samples,
            n_guides,
        })?;
    let mut flat_data = Vec::with_capacity(total_len);
    for sample in sample_counts {
        flat_data.extend_from_slice(sample);
    }
    CountMatrix::from_vec(flat_data, n_samples, n_guides)
}

/// Process multiple FASTQ files sequentially, assembling results into a [`CountMatrix`].
///
/// # Arguments
/// - `fastq_paths`: Slice of paths to FASTQ files.
/// - `index`: Reference guide index.
/// - `sgrna_start`: 0-based offset from start of read to guide sequence.
/// - `sgrna_len`: Length of guide sequence to extract.
/// - `one_mismatch`: Whether 1-mismatch matching is enabled.
///
/// # Returns
/// - Sample-major [`CountMatrix`] with shape `[n_samples, n_guides]`.
/// - Vector of [`ReadSummary`] structs, one per FASTQ file in `fastq_paths` order.
pub fn count_all_fastq<P: AsRef<Path>>(
    fastq_paths: &[P],
    index: &GuideIndex,
    sgrna_start: usize,
    sgrna_len: usize,
    one_mismatch: bool,
) -> Result<(CountMatrix, Vec<ReadSummary>), MageckError> {
    let mut sample_counts = Vec::with_capacity(fastq_paths.len());
    let mut summaries = Vec::with_capacity(fastq_paths.len());

    for path in fastq_paths {
        let (counts, summary) =
            process_one_fastq(path, index, sgrna_start, sgrna_len, one_mismatch)?;
        sample_counts.push(counts);
        summaries.push(summary);
    }

    let matrix = assemble_count_matrix(&sample_counts, index.len())
        .map_err(|e| MageckError::Fastq(format!("Failed to assemble count matrix: {e}")))?;

    Ok((matrix, summaries))
}

/// Safely accumulates counts from a single physical file into the target sample slice.
///
/// Validates dimensional compatibility, ensures values remain finite and non-negative integers,
/// and checks that accumulated values do not exceed 2^53 (exact integer representation in f64).
fn accumulate_counts_into_slice(
    target_slice: &mut [f64],
    file_counts: &[f64],
) -> Result<(), MageckError> {
    if file_counts.len() != target_slice.len() {
        return Err(MageckError::Fastq(format!(
            "Count length mismatch: file produced {} counts, expected {}",
            file_counts.len(),
            target_slice.len()
        )));
    }

    // Pass 1: Validate all operands and results
    for (&dest, &src) in target_slice.iter().zip(file_counts) {
        if !src.is_finite() || src < 0.0 || src.fract() != 0.0 {
            return Err(MageckError::Fastq(format!(
                "Count value is invalid (must be a finite, non-negative integer): {src}"
            )));
        }
        if !dest.is_finite() || dest < 0.0 || dest.fract() != 0.0 {
            return Err(MageckError::Fastq(format!(
                "Existing target value is invalid (must be a finite, non-negative integer): {dest}"
            )));
        }
        if dest > MAX_SAFE_INT_F64 - src {
            return Err(MageckError::Fastq(
                "Accumulated count exceeds exact f64 integer representation (2^53)".to_string(),
            ));
        }
        let sum = dest + src;
        if !sum.is_finite() || sum < 0.0 || sum.fract() != 0.0 {
            return Err(MageckError::Fastq(format!(
                "Accumulated count is invalid (must be a finite, non-negative integer): {sum}"
            )));
        }
        if sum > MAX_SAFE_INT_F64 {
            return Err(MageckError::Fastq(format!(
                "Accumulated count {sum} exceeds exact f64 integer representation (2^53)"
            )));
        }
    }

    // Pass 2: Mutate only after complete validation succeeds
    for (dest, &src) in target_slice.iter_mut().zip(file_counts) {
        *dest += src;
    }

    Ok(())
}

/// Sequentially processes grouped FASTQ files, aggregating counts per sample group.
///
/// # Invariants
/// - Sample groups become distinct rows in the returned [`CountMatrix`].
/// - Technical replicates within each sample group have their raw counts added exactly.
/// - Each physical FASTQ produces an independent [`PhysicalFastqSummary`].
/// - Per-file `zerocounts` is computed from that file's unaggregated guide counts.
/// - Summaries preserve sample-group order and physical-file order within each group.
/// - If any file fails to parse or process, returns an error without partial success.
pub fn count_grouped_fastq(
    sample_groups: &[FastqSampleGroup],
    index: &GuideIndex,
    sgrna_start: usize,
    sgrna_len: usize,
    one_mismatch: bool,
) -> Result<GroupedCountResult, MageckError> {
    count_grouped_fastq_with_trim(
        sample_groups,
        index,
        Trim5Mode::Fixed(sgrna_start),
        sgrna_len,
        one_mismatch,
    )
}

/// Grouped counting with per-read auto or fixed offset selection.
pub(crate) fn count_grouped_fastq_with_trim(
    sample_groups: &[FastqSampleGroup],
    index: &GuideIndex,
    trim_5: Trim5Mode,
    sgrna_len: usize,
    one_mismatch: bool,
) -> Result<GroupedCountResult, MageckError> {
    let sgrna_start = match trim_5 {
        Trim5Mode::Fixed(offset) => offset,
        Trim5Mode::None | Trim5Mode::Auto => 0,
    };
    if sgrna_len == 0 {
        return Err(MageckError::Fastq(
            "sgRNA length must be greater than 0".to_string(),
        ));
    }
    if sgrna_start.checked_add(sgrna_len).is_none() {
        return Err(MageckError::Fastq(format!(
            "Window range overflow: start ({sgrna_start}) + len ({sgrna_len}) exceeds usize::MAX"
        )));
    }

    let n_groups = sample_groups.len();
    let n_guides = index.len();
    let total_elements = n_groups
        .checked_mul(n_guides)
        .ok_or_else(|| MageckError::Fastq(format!(
            "CountMatrix dimension overflow: n_groups ({n_groups}) * n_guides ({n_guides}) overflows usize"
        )))?;

    let mut matrix = CountMatrix::from_vec(vec![0.0; total_elements], n_groups, n_guides)
        .map_err(|e| MageckError::Fastq(format!("Failed to initialize CountMatrix: {e}")))?;

    let total_files: usize = sample_groups.iter().map(|g| g.fastq_paths.len()).sum();
    let mut summaries = Vec::with_capacity(total_files);

    for (group_idx, group) in sample_groups.iter().enumerate() {
        let group_slice = matrix.get_sample_slice_mut(group_idx).ok_or_else(|| {
            MageckError::Fastq(format!("Sample group index {group_idx} out of bounds"))
        })?;

        for path in &group.fastq_paths {
            let (file_counts, summary) =
                process_one_fastq_with_trim(path, index, trim_5, sgrna_len, one_mismatch)?;

            accumulate_counts_into_slice(group_slice, &file_counts)?;

            let zerocounts = file_counts.iter().filter(|&&c| c == 0.0).count();

            summaries.push(PhysicalFastqSummary {
                path: path.clone(),
                label: group.label.clone(),
                summary,
                zerocounts,
            });
        }
    }

    Ok(GroupedCountResult {
        counts: matrix,
        summaries,
    })
}

/// UMI QC records retain physical-file attribution independently of read QC.
#[derive(Debug, Clone)]
pub(crate) struct UmiFileSummary {
    pub path: PathBuf,
    pub label: String,
    pub valid_mapped_reads: usize,
    pub new_molecules: usize,
    pub duplicate_molecules: usize,
}

fn increment(counter: &mut usize) -> Result<(), MageckError> {
    *counter = counter
        .checked_add(1)
        .ok_or_else(|| MageckError::Fastq("Read counter overflow".into()))?;
    Ok(())
}

/// Opt-in sequential processor. Ordinary COUNT continues through the original reader.
pub(crate) fn count_grouped_with_options(
    config: &crate::config::CountConfig,
    index: &GuideIndex,
    mut unmapped: Option<&mut dyn std::io::Write>,
) -> Result<(GroupedCountResult, Vec<UmiFileSummary>), MageckError> {
    use crate::cli::UmiDedupScope;
    use crate::data::umi::UmiDeduper;
    if let Some(umi) = &config.umi {
        umi.validate()?;
    }
    let guide_start = match config.trim_5 {
        Trim5Mode::Fixed(n) => n,
        _ => 0,
    };
    if config.sgrna_len == 0 || guide_start.checked_add(config.sgrna_len).is_none() {
        return Err(MageckError::Fastq(
            "Invalid guide length or coordinate addition overflow".into(),
        ));
    }
    let mut sample_counts = Vec::new();
    let mut summaries = Vec::new();
    let mut umi_summaries = Vec::new();
    for group in &config.sample_groups {
        let mut counts = vec![0.0; index.len()];
        let mut deduper = config
            .umi
            .as_ref()
            .map(|u| UmiDeduper::new(index.len(), u.length));
        for path in &group.fastq_paths {
            if config
                .umi
                .as_ref()
                .is_some_and(|u| u.scope == UmiDedupScope::File)
            {
                deduper = config
                    .umi
                    .as_ref()
                    .map(|u| UmiDeduper::new(index.len(), u.length));
            }
            let mut summary = ReadSummary::default();
            let mut seen_guides = vec![false; index.len()];
            let mut molecules = UmiFileSummary {
                path: path.clone(),
                label: group.label.clone(),
                valid_mapped_reads: 0,
                new_molecules: 0,
                duplicate_molecules: 0,
            };
            let mut reader = match needletail::parse_fastx_file(path) {
                Ok(r) => Some(r),
                Err(e) if e.kind == ParseErrorKind::EmptyFile => None,
                Err(e) => {
                    return Err(MageckError::Fastq(format!(
                        "Failed to open FASTQ file '{}': {e}",
                        path.display()
                    )))
                }
            };
            if let Some(reader) = reader.as_mut() {
                while let Some(record) = reader.next() {
                    let record = record.map_err(|e| {
                        MageckError::Fastq(format!(
                            "Error parsing FASTQ record in '{}': {e}",
                            path.display()
                        ))
                    })?;
                    increment(&mut summary.total_reads)?;
                    let error = |reason: String| {
                        MageckError::Fastq(format!(
                            "'{}' record {}: {reason}",
                            path.display(),
                            summary.total_reads
                        ))
                    };
                    let quality = record.qual().ok_or_else(|| error("Expected FASTQ quality values; FASTA is unsupported for UMI/unmapped workflows".into()))?;
                    let seq = record.seq();
                    let (result, offset) = if config.trim_5 == Trim5Mode::Auto {
                        let mut selected = (MatchResult::NoMatch, 0);
                        for offset in 0..=40 {
                            let Some(window) = extract_window(&seq, offset, config.sgrna_len)
                            else {
                                break;
                            };
                            if let result @ MatchResult::Exact(_) = index.lookup(window) {
                                selected = (result, offset);
                                break;
                            }
                        }
                        selected
                    } else if let Some(window) = extract_window(&seq, guide_start, config.sgrna_len)
                    {
                        (index.lookup(window), guide_start)
                    } else {
                        increment(&mut summary.short_reads)?;
                        if let Some(writer) = unmapped.as_mut() {
                            write_unmapped(*writer, record.id(), &seq, quality)?;
                        }
                        continue;
                    };
                    let mapped = match result {
                        MatchResult::Exact(i) => {
                            increment(&mut summary.exact_matches)?;
                            Some(i)
                        }
                        MatchResult::OneMismatch(i) if config.test_mismatch == 1 => {
                            increment(&mut summary.one_mismatch_matches)?;
                            Some(i)
                        }
                        MatchResult::Ambiguous => {
                            increment(&mut summary.ambiguous_reads)?;
                            None
                        }
                        _ => {
                            increment(&mut summary.unmapped_reads)?;
                            None
                        }
                    };
                    if let Some(i) = mapped {
                        increment(&mut summary.mapped_reads)?;
                        seen_guides[i] = true;
                        let is_new =
                            if let (Some(umi), Some(deduper)) = (&config.umi, deduper.as_mut()) {
                                let end = umi.start.checked_add(umi.length).ok_or_else(|| {
                                    error("UMI coordinate addition overflow".into())
                                })?;
                                let guide_end =
                                    offset.checked_add(config.sgrna_len).ok_or_else(|| {
                                        error("Guide coordinate addition overflow".into())
                                    })?;
                                if umi.start < guide_end && offset < end {
                                    return Err(error(
                                        "UMI window overlaps the selected guide window".into(),
                                    ));
                                }
                                let sequence = extract_window(&seq, umi.start, umi.length)
                                    .ok_or_else(|| {
                                        error("Mapped read is too short for the UMI window".into())
                                    })?;
                                let is_new = deduper
                                    .observe_capped(i, sequence, umi.max_keys)
                                    .map_err(|e| error(e.to_string()))?;
                                increment(&mut molecules.valid_mapped_reads)?;
                                if is_new {
                                    increment(&mut molecules.new_molecules)?;
                                } else {
                                    increment(&mut molecules.duplicate_molecules)?;
                                }
                                is_new
                            } else {
                                true
                            };
                        if is_new {
                            if counts[i] >= MAX_SAFE_INT_F64 {
                                return Err(error(
                                    "Molecule/read count exceeds exact f64 integer range".into(),
                                ));
                            }
                            counts[i] += 1.0;
                        }
                    } else if let Some(writer) = unmapped.as_mut() {
                        write_unmapped(*writer, record.id(), &seq, quality)?;
                    }
                }
            }
            summaries.push(PhysicalFastqSummary {
                path: path.clone(),
                label: group.label.clone(),
                summary,
                zerocounts: seen_guides.iter().filter(|&&seen| !seen).count(),
            });
            if config.umi.is_some() {
                umi_summaries.push(molecules);
            }
        }
        sample_counts.push(counts);
    }
    let counts = assemble_count_matrix(&sample_counts, index.len())
        .map_err(|e| MageckError::Fastq(e.to_string()))?;
    Ok((GroupedCountResult { counts, summaries }, umi_summaries))
}

fn write_unmapped(
    writer: &mut dyn std::io::Write,
    id: &[u8],
    seq: &[u8],
    quality: &[u8],
) -> Result<(), MageckError> {
    for bytes in [b"@".as_slice(), id, b"\n", seq, b"\n+\n", quality, b"\n"] {
        writer.write_all(bytes)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Write};

    use crate::count::matcher::{build_index, parse_library, Library};
    use tempfile::NamedTempFile;

    #[test]
    fn unmapped_writer_io_failure_propagates() {
        struct RejectWrite;
        impl std::io::Write for RejectWrite {
            fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
                Err(std::io::Error::other("injected write failure"))
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let error = write_unmapped(&mut RejectWrite, b"id", b"ACGT", b"IIII").unwrap_err();
        assert!(matches!(error, MageckError::Io(_)));
        assert!(error.to_string().contains("injected write failure"));
    }

    const STD_LIB_TSV: &str = "id\tsequence\tgene\n\
    g0\tACGTACGTACGTACGTACGT\tGENE_A\n\
    g1\tTGCATGCATGCATGCATGCA\tGENE_B\n\
    g2\tGGGGCCCCAAAATTTTGGGG\tGENE_C\n";

    fn make_library(tsv: &str) -> Library {
        parse_library(Cursor::new(tsv.as_bytes())).expect("parse library")
    }

    fn make_temp_fastq(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().expect("create temp fastq");
        file.write_all(content.as_bytes())
            .expect("write temp fastq");
        file.flush().expect("flush temp fastq");
        file
    }

    #[test]
    fn auto_trim_matches_variable_offsets_including_bounds_and_lowercase() {
        let index = build_index(&make_library(STD_LIB_TSV), false).unwrap();
        let mut data = String::new();
        for offset in [0, 5, 23, 40] {
            let seq = format!("{}acgtacgtacgtacgtacgt", "N".repeat(offset));
            data += &format!("@r\n{seq}\n+\n{}\n", "I".repeat(seq.len()));
        }
        let fastq = make_temp_fastq(&data);
        let (counts, summary) =
            process_one_fastq_with_trim(fastq.path(), &index, Trim5Mode::Auto, 20, false).unwrap();
        assert_eq!(counts, vec![4.0, 0.0, 0.0]);
        assert_eq!(summary.mapped_reads, 4);
    }

    #[test]
    fn auto_trim_uses_first_exact_match_and_counts_once() {
        let index = build_index(&make_library(STD_LIB_TSV), true).unwrap();
        let seq = b"ACGTACGTACGTACGTACGTTGCATGCATGCATGCATGCA";
        assert_eq!(match_auto_window(seq, &index, 20), MatchResult::Exact(0));
        let seq = b"ACGTACGTACGTACGTACGATGCATGCATGCATGCATGCA";
        assert_eq!(match_auto_window(seq, &index, 20), MatchResult::Exact(1));
        let fastq = make_temp_fastq("@r\nACGTACGTACGTACGTACGTTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
        let (counts, summary) =
            process_one_fastq_with_trim(fastq.path(), &index, Trim5Mode::Auto, 20, true).unwrap();
        assert_eq!(counts, vec![1.0, 0.0, 0.0]);
        assert_eq!(summary.mapped_reads, 1);
    }

    #[test]
    fn auto_trim_no_exact_match_is_unmapped_and_offset_41_is_excluded() {
        let index = build_index(&make_library(STD_LIB_TSV), true).unwrap();
        let seq = format!("{}ACGTACGTACGTACGTACGT", "N".repeat(41));
        let data = format!("@outside\n{seq}\n+\n{}\n@short\nACG\n+\nIII\n@mismatch\nACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIII\n", "I".repeat(seq.len()));
        let fastq = make_temp_fastq(&data);
        let (counts, summary) =
            process_one_fastq_with_trim(fastq.path(), &index, Trim5Mode::Auto, 20, true).unwrap();
        assert_eq!(counts, vec![0.0, 0.0, 0.0]);
        assert_eq!(summary.unmapped_reads, 3);
        assert_eq!(summary.short_reads, 0);
        assert_eq!(summary.one_mismatch_matches, 0);
    }

    #[test]
    fn fixed_trim_extracts_guide_and_none_preserves_zero_offset() {
        let index = build_index(&make_library(STD_LIB_TSV), false).unwrap();
        let fastq =
            make_temp_fastq("@r\nNNNNNACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIII\n");
        let fixed = process_one_fastq(fastq.path(), &index, 5, 20, false).unwrap();
        assert_eq!(fixed.0, vec![1.0, 0.0, 0.0]);
        assert_eq!(fixed.1.mapped_reads, 1);
        assert_eq!(
            process_one_fastq_with_trim(fastq.path(), &index, Trim5Mode::None, 20, false).unwrap(),
            process_one_fastq(fastq.path(), &index, 0, 20, false).unwrap()
        );
        assert_eq!(
            process_one_fastq(fastq.path(), &index, 0, 20, false)
                .unwrap()
                .1
                .mapped_reads,
            0
        );
    }

    #[test]
    fn test_extract_window_direct() {
        let seq = b"NNNNNACGTACGTACGTACGTACGT"; // 5 nt prefix + 20 nt guide = 25 nt

        // Offset 0, length 5
        assert_eq!(extract_window(seq, 0, 5), Some(&b"NNNNN"[..]));

        // Offset 5, length 20 (ends at exact end of slice)
        assert_eq!(
            extract_window(seq, 5, 20),
            Some(&b"ACGTACGTACGTACGTACGT"[..])
        );

        // Off-by-one: offset 5, length 21 (exceeds slice length 25)
        assert_eq!(extract_window(seq, 5, 21), None);

        // Offset beyond slice
        assert_eq!(extract_window(seq, 26, 1), None);

        // Integer overflow prevention
        assert_eq!(extract_window(seq, usize::MAX, 1), None);
        assert_eq!(extract_window(seq, 10, usize::MAX), None);

        // Empty sequence
        assert_eq!(extract_window(b"", 0, 1), None);
        assert_eq!(extract_window(b"", 0, 0), Some(&b""[..]));
    }
    #[test]
    fn test_process_fastq_stream_cursor() {
        let lib = make_library(STD_LIB_TSV);
        let index = build_index(&lib, false).expect("build index");

        let fq_data = b"@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
        let cursor = Cursor::new(fq_data);

        let (counts, summary) =
            process_fastq_stream(cursor, &index, 0, 20, false).expect("process fastq stream");

        assert_eq!(counts, vec![1.0, 0.0, 0.0]);
        assert_eq!(summary.exact_matches, 1);
    }
    #[test]
    fn test_assemble_count_matrix() {
        let s0 = vec![1.0, 2.0, 3.0];
        let s1 = vec![4.0, 5.0, 6.0];

        let matrix = assemble_count_matrix(&[s0, s1], 3).expect("assemble matrix");

        assert_eq!(matrix.n_samples(), 2);
        assert_eq!(matrix.n_guides(), 3);
        assert_eq!(matrix.sample_slice(0), &[1.0, 2.0, 3.0]);
        assert_eq!(matrix.sample_slice(1), &[4.0, 5.0, 6.0]);
    }

    #[test]
    fn test_summary_invariants_direct() {
        let lib_csv = "id,sequence,gene\n\
                       g1,AAAAAAAAAA,gene1\n\
                       g2,TAAAAAAAAA,gene2\n\
                       g3,CCCCCCCCCC,gene3\n";
        let lib = parse_library(Cursor::new(lib_csv)).unwrap();
        let index = build_index(&lib, true).unwrap();
        let fq = "@r1\nAAAAAAAAAA\n+\nIIIIIIIIII\n\
                  @r2\nCCCCCCCCTC\n+\nIIIIIIIIII\n\
                  @r3\nCAAAAAAAAA\n+\nIIIIIIIIII\n\
                  @r4\nGGGGGGGGGG\n+\nIIIIIIIIII\n\
                  @r5\nATGC\n+\nIIII\n";
        let tmp = make_temp_fastq(fq);
        let (counts, summary) = process_one_fastq(tmp.path(), &index, 0, 10, true).unwrap();

        assert_eq!(summary.total_reads, 5);
        assert_eq!(summary.exact_matches, 1);
        assert_eq!(summary.one_mismatch_matches, 1);
        assert_eq!(summary.ambiguous_reads, 1);
        assert_eq!(summary.unmapped_reads, 1);
        assert_eq!(summary.short_reads, 1);
        assert_eq!(
            summary.mapped_reads,
            summary.exact_matches + summary.one_mismatch_matches
        );
        assert_eq!(
            summary.total_reads,
            summary.mapped_reads
                + summary.ambiguous_reads
                + summary.unmapped_reads
                + summary.short_reads
        );

        assert_eq!(counts, vec![1.0, 0.0, 1.0]);
    }

    #[test]
    fn test_accumulate_counts_into_slice_success() {
        let mut target = vec![10.0, 20.0, 30.0];
        let file = vec![1.0, 2.0, 3.0];
        accumulate_counts_into_slice(&mut target, &file).expect("accumulate success");
        assert_eq!(target, vec![11.0, 22.0, 33.0]);
    }

    #[test]
    fn test_accumulate_counts_into_slice_dimension_mismatch() {
        let mut target = vec![10.0, 20.0];
        let file = vec![1.0, 2.0, 3.0];
        let err = accumulate_counts_into_slice(&mut target, &file).unwrap_err();
        match err {
            MageckError::Fastq(msg) => {
                assert!(msg.contains("Count length mismatch"));
            }
            other => panic!("expected Fastq error, got {other:?}"),
        }
    }

    #[test]
    fn test_accumulate_counts_into_slice_rejects_invalid_values() {
        // Negative count
        let mut target = vec![0.0];
        let file = vec![-1.0];
        assert!(matches!(
            accumulate_counts_into_slice(&mut target, &file),
            Err(MageckError::Fastq(_))
        ));

        // NaN
        let mut target = vec![0.0];
        let file = vec![f64::NAN];
        assert!(matches!(
            accumulate_counts_into_slice(&mut target, &file),
            Err(MageckError::Fastq(_))
        ));

        // Infinity
        let mut target = vec![0.0];
        let file = vec![f64::INFINITY];
        assert!(matches!(
            accumulate_counts_into_slice(&mut target, &file),
            Err(MageckError::Fastq(_))
        ));

        // Non-integral
        let mut target = vec![0.0];
        let file = vec![1.5];
        assert!(matches!(
            accumulate_counts_into_slice(&mut target, &file),
            Err(MageckError::Fastq(_))
        ));

        // Exceeding MAX_SAFE_INT_F64
        let mut target = vec![MAX_SAFE_INT_F64];
        let file = vec![1.0];
        assert!(matches!(
            accumulate_counts_into_slice(&mut target, &file),
            Err(MageckError::Fastq(_))
        ));

        // Exactly MAX_SAFE_INT_F64 accepted
        let mut target = vec![MAX_SAFE_INT_F64 - 1.0];
        let file = vec![1.0];
        assert!(accumulate_counts_into_slice(&mut target, &file).is_ok());
        assert_eq!(target, vec![MAX_SAFE_INT_F64]);

        // Above MAX_SAFE_INT_F64 rejected
        let mut target = vec![0.0];
        let file = vec![MAX_SAFE_INT_F64 + 2.0];
        assert!(matches!(
            accumulate_counts_into_slice(&mut target, &file),
            Err(MageckError::Fastq(_))
        ));
    }

    #[test]
    fn test_accumulate_counts_into_slice_atomicity() {
        // Prove that the target slice is unchanged if a later element is invalid
        let mut target = vec![10.0, 20.0, 30.0];
        let file = vec![1.0, f64::NAN, 3.0]; // First element valid, second invalid
        let err = accumulate_counts_into_slice(&mut target, &file).unwrap_err();
        assert!(matches!(err, MageckError::Fastq(_)));

        // Target should be completely unmodified, even the first element
        assert_eq!(target, vec![10.0, 20.0, 30.0]);
    }

    #[test]
    fn test_count_grouped_fastq_unit() {
        let lib = make_library(STD_LIB_TSV);
        let index = build_index(&lib, false).expect("build index");

        let fq1 = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
        let fq2 = "@r2\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n";
        let tmp1 = make_temp_fastq(fq1);
        let tmp2 = make_temp_fastq(fq2);

        let group = FastqSampleGroup {
            fastq_paths: vec![tmp1.path().to_path_buf(), tmp2.path().to_path_buf()],
            label: "group1".to_string(),
        };

        let result = count_grouped_fastq(&[group], &index, 0, 20, false).expect("count grouped");

        assert_eq!(result.counts.n_samples(), 1);
        assert_eq!(result.counts.n_guides(), 3);
        assert_eq!(result.counts.sample_slice(0), &[1.0, 1.0, 0.0]);
        assert_eq!(result.summaries.len(), 2);
        assert_eq!(result.summaries[0].label, "group1");
        assert_eq!(result.summaries[1].label, "group1");
        assert_eq!(result.summaries[0].path, tmp1.path());
        assert_eq!(result.summaries[1].path, tmp2.path());
        assert_eq!(result.summaries[0].zerocounts, 2);
        assert_eq!(result.summaries[1].zerocounts, 2);
    }
}
