//! sgRNA count-table parser, validated CountTable representation, and deterministic sample resolution.
//!
//! Implements Phase 4B deliverables:
//! - [`CountTable`]: Validated in-memory count table with metadata and resolved sample indices.
//! - [`parse_count_table`]: Streaming TSV count-table parser.
//! - [`parse_count_table_file`]: File-based count-table parser.
//! - [`resolve_sample_indices`]: Deterministic sample resolution adhering to spec hierarchy.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use crate::data::count_matrix::CountMatrix;
use crate::error::MageckError;

/// Parsed and validated sgRNA count table with metadata and resolved sample indices.
#[derive(Debug, Clone, PartialEq)]
pub struct CountTable {
    sgrna_ids: Vec<String>,
    genes: Vec<String>,
    sample_labels: Vec<String>,
    counts: CountMatrix,
    treatment_indices: Vec<usize>,
    control_indices: Vec<usize>,
}

/// Decomposed components of a [`CountTable`]:
/// `(sgrna_ids, genes, sample_labels, counts, treatment_indices, control_indices)`.
pub type CountTableParts = (
    Vec<String>,
    Vec<String>,
    Vec<String>,
    CountMatrix,
    Vec<usize>,
    Vec<usize>,
);

impl CountTable {
    /// Constructs a validated [`CountTable`] from individual validated components.
    pub fn from_parts(
        sgrna_ids: Vec<String>,
        genes: Vec<String>,
        sample_labels: Vec<String>,
        counts: CountMatrix,
        treatment_indices: Vec<usize>,
        control_indices: Vec<usize>,
    ) -> Result<Self, MageckError> {
        let n_guides = sgrna_ids.len();
        let n_samples = sample_labels.len();

        if n_guides == 0 {
            return Err(MageckError::Config(
                "Count table must contain at least one sgRNA record".to_string(),
            ));
        }
        if genes.len() != n_guides {
            return Err(MageckError::Config(format!(
                "Gene identifier count ({}) does not match sgRNA count ({n_guides})",
                genes.len()
            )));
        }
        if counts.n_guides() != n_guides {
            return Err(MageckError::Config(format!(
                "Count matrix guide count ({}) does not match sgRNA count ({n_guides})",
                counts.n_guides()
            )));
        }
        if n_samples == 0 {
            return Err(MageckError::Config(
                "Count table must contain at least one sample column".to_string(),
            ));
        }
        if counts.n_samples() != n_samples {
            return Err(MageckError::Config(format!(
                "Count matrix sample count ({}) does not match sample label count ({n_samples})",
                counts.n_samples()
            )));
        }

        // Validate sample labels
        let mut seen_labels = HashSet::with_capacity(n_samples);
        for (i, label) in sample_labels.iter().enumerate() {
            if label.is_empty() {
                return Err(MageckError::Config(format!(
                    "Sample label at index {i} must not be empty"
                )));
            }
            if !seen_labels.insert(label.clone()) {
                return Err(MageckError::Config(format!(
                    "Duplicate sample column label: '{label}'"
                )));
            }
        }

        // Validate sgRNAs and genes
        let mut seen_sgrnas = HashMap::with_capacity(n_guides);
        for (i, (sgrna, gene)) in sgrna_ids.iter().zip(genes.iter()).enumerate() {
            let row = i + 1;
            if sgrna.is_empty() {
                return Err(MageckError::Config(format!(
                    "Row {row}: sgRNA identifier is empty"
                )));
            }
            if gene.is_empty() {
                return Err(MageckError::Config(format!(
                    "Row {row}: gene identifier is empty for sgRNA '{sgrna}'"
                )));
            }
            if let Some(&first_row) = seen_sgrnas.get(sgrna) {
                return Err(MageckError::Config(format!(
                    "Row {row}: duplicate sgRNA identifier '{sgrna}' (first seen at row {first_row})"
                )));
            }
            seen_sgrnas.insert(sgrna.clone(), row);
        }

        // Validate matrix values: non-negative and finite
        for &val in counts.as_slice() {
            if !val.is_finite() || val < 0.0 {
                return Err(MageckError::Config(format!(
                    "Count matrix contains invalid count {val}: must be finite and non-negative"
                )));
            }
        }

        // Validate treatment indices
        if treatment_indices.is_empty() {
            return Err(MageckError::Config(
                "Treatment sample list must not be empty".to_string(),
            ));
        }
        let mut seen_treat = HashSet::with_capacity(treatment_indices.len());
        for &idx in &treatment_indices {
            if idx >= n_samples {
                return Err(MageckError::Config(format!(
                    "Treatment sample index {idx} out of range (0..{n_samples})"
                )));
            }
            if !seen_treat.insert(idx) {
                let label = match sample_labels.get(idx) {
                    Some(s) => s.as_str(),
                    None => "",
                };
                return Err(MageckError::Config(format!(
                    "Duplicate treatment sample selection: sample '{label}' (index {idx}) selected more than once"
                )));
            }
        }

        // Validate control indices
        if control_indices.is_empty() {
            return Err(MageckError::Config(
                "Control sample list must not be empty".to_string(),
            ));
        }
        let mut seen_ctrl = HashSet::with_capacity(control_indices.len());
        for &idx in &control_indices {
            if idx >= n_samples {
                return Err(MageckError::Config(format!(
                    "Control sample index {idx} out of range (0..{n_samples})"
                )));
            }
            if !seen_ctrl.insert(idx) {
                let label = match sample_labels.get(idx) {
                    Some(s) => s.as_str(),
                    None => "",
                };
                return Err(MageckError::Config(format!(
                    "Duplicate control sample selection: sample '{label}' (index {idx}) selected more than once"
                )));
            }
        }

        // Validate no overlap
        for &idx in &treatment_indices {
            if seen_ctrl.contains(&idx) {
                let label = match sample_labels.get(idx) {
                    Some(s) => s.as_str(),
                    None => "",
                };
                return Err(MageckError::Config(format!(
                    "Sample '{label}' (index {idx}) is selected in both treatment and control groups"
                )));
            }
        }

        Ok(Self {
            sgrna_ids,
            genes,
            sample_labels,
            counts,
            treatment_indices,
            control_indices,
        })
    }

    /// Constructs a validated [`CountTable`] by resolving treatment and control tokens against sample labels.
    pub fn new<S1: AsRef<str>, S2: AsRef<str>>(
        sgrna_ids: Vec<String>,
        genes: Vec<String>,
        sample_labels: Vec<String>,
        counts: CountMatrix,
        treatment_ids: &[S1],
        control_ids: &[S2],
    ) -> Result<Self, MageckError> {
        let (treatment_indices, control_indices) =
            resolve_sample_indices(&sample_labels, treatment_ids, control_ids)?;
        Self::from_parts(
            sgrna_ids,
            genes,
            sample_labels,
            counts,
            treatment_indices,
            control_indices,
        )
    }

    /// Returns the ordered sgRNA identifiers.
    #[inline]
    pub fn sgrna_ids(&self) -> &[String] {
        &self.sgrna_ids
    }

    /// Alias for [`sgrna_ids`](Self::sgrna_ids), returning the ordered sgRNA identifiers.
    #[inline]
    pub fn sgrnas(&self) -> &[String] {
        &self.sgrna_ids
    }

    /// Returns the ordered gene identifiers, aligned 1-to-1 with guides.
    #[inline]
    pub fn genes(&self) -> &[String] {
        &self.genes
    }

    /// Returns the ordered sample labels from the table header.
    #[inline]
    pub fn sample_labels(&self) -> &[String] {
        &self.sample_labels
    }

    /// Returns a reference to the underlying [`CountMatrix`].
    #[inline]
    pub fn counts(&self) -> &CountMatrix {
        &self.counts
    }

    /// Returns the resolved 0-based sample indices for treatment samples in user-specified order.
    #[inline]
    pub fn treatment_indices(&self) -> &[usize] {
        &self.treatment_indices
    }

    /// Returns the resolved 0-based sample indices for control samples in user-specified order.
    #[inline]
    pub fn control_indices(&self) -> &[usize] {
        &self.control_indices
    }

    /// Returns the number of guides (rows) in the table.
    #[inline]
    pub fn n_guides(&self) -> usize {
        self.sgrna_ids.len()
    }

    /// Returns the number of samples (columns) in the table.
    #[inline]
    pub fn n_samples(&self) -> usize {
        self.sample_labels.len()
    }

    /// Returns `true` if the table contains zero guides.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sgrna_ids.is_empty()
    }

    /// Returns the 0-based column index of a sample label, if found.
    pub fn sample_index(&self, label: &str) -> Option<usize> {
        self.sample_labels.iter().position(|l| l == label)
    }

    /// Consumes the table, returning its components:
    /// `(sgrna_ids, genes, sample_labels, counts, treatment_indices, control_indices)`.
    pub fn into_parts(self) -> CountTableParts {
        (
            self.sgrna_ids,
            self.genes,
            self.sample_labels,
            self.counts,
            self.treatment_indices,
            self.control_indices,
        )
    }

    /// Consumes the table, returning the underlying [`CountMatrix`].
    #[inline]
    pub fn into_counts(self) -> CountMatrix {
        self.counts
    }
}

/// Resolves a single sample selection token against sample labels according to the spec hierarchy.
fn resolve_single_token(
    token: &str,
    sample_labels: &[String],
    group_name: &str,
) -> Result<usize, MageckError> {
    if token.is_empty() {
        return Err(MageckError::Config(format!(
            "{group_name} sample selection token must not be empty"
        )));
    }

    // 1. Resolve each user token against an exact sample-column label first.
    // 2. Matching is case-sensitive.
    // 12. A numeric sample label such as "0" takes precedence over interpreting "0" as index zero.
    if let Some(pos) = sample_labels.iter().position(|l| l == token) {
        return Ok(pos);
    }

    // 14. Leading or trailing whitespace in a CLI selection token must not silently change which sample is selected.
    if token.starts_with(char::is_whitespace) || token.ends_with(char::is_whitespace) {
        return Err(MageckError::Config(format!(
            "{group_name} sample selection token '{token}' contains leading or trailing whitespace and does not match any sample label"
        )));
    }

    // 3. If and only if no exact sample-label match exists, parse the complete token as a base-10 usize.
    if token.chars().all(|c| c.is_ascii_digit()) {
        if let Ok(idx) = token.parse::<usize>() {
            // 4. Numeric indices are zero-based within sample columns only.
            // 5. Required metadata columns such as sgRNA and Gene are not included in index numbering.
            // 6. Reject an index outside the sample-column range.
            if idx < sample_labels.len() {
                return Ok(idx);
            } else {
                return Err(MageckError::Config(format!(
                    "{group_name} sample index {idx} is out of range (count table has {} sample columns, valid indices 0..{})",
                    sample_labels.len(),
                    sample_labels.len().saturating_sub(1)
                )));
            }
        }
    }

    // 7. Reject a token that is neither an exact label nor a valid in-range index.
    Err(MageckError::Config(format!(
        "{group_name} sample selection token '{token}' does not match any sample label or valid 0-based sample column index"
    )))
}

/// Resolves treatment and control selection tokens against sample column labels.
///
/// Implements the deterministic resolution hierarchy:
/// 1. Exact sample-column label match first (case-sensitive).
/// 2. Numeric sample label takes precedence over 0-based column index.
/// 3. If no label match, parses complete token as base-10 `usize` (0-based within sample columns).
/// 4. Rejects out-of-range indices and invalid tokens.
/// 5. Rejects duplicate treatment selections and duplicate control selections.
/// 6. Rejects any sample selected in both treatment and control groups.
/// 7. Preserves user-specified order in resolved index vectors.
pub fn resolve_sample_indices<S1: AsRef<str>, S2: AsRef<str>>(
    sample_labels: &[String],
    treatment_ids: &[S1],
    control_ids: &[S2],
) -> Result<(Vec<usize>, Vec<usize>), MageckError> {
    if treatment_ids.is_empty() {
        return Err(MageckError::Config(
            "Treatment sample list must not be empty".to_string(),
        ));
    }
    if control_ids.is_empty() {
        return Err(MageckError::Config(
            "Control sample list must not be empty".to_string(),
        ));
    }

    let mut treatment_indices = Vec::with_capacity(treatment_ids.len());
    let mut seen_treatment = HashSet::with_capacity(treatment_ids.len());

    for token in treatment_ids {
        let t_str = token.as_ref();
        let idx = resolve_single_token(t_str, sample_labels, "Treatment")?;
        if !seen_treatment.insert(idx) {
            let label = match sample_labels.get(idx) {
                Some(s) => s.as_str(),
                None => "",
            };
            return Err(MageckError::Config(format!(
                "Duplicate treatment sample selection: sample '{label}' (index {idx}) was selected more than once"
            )));
        }
        treatment_indices.push(idx);
    }

    let mut control_indices = Vec::with_capacity(control_ids.len());
    let mut seen_control = HashSet::with_capacity(control_ids.len());

    for token in control_ids {
        let c_str = token.as_ref();
        let idx = resolve_single_token(c_str, sample_labels, "Control")?;
        if !seen_control.insert(idx) {
            let label = match sample_labels.get(idx) {
                Some(s) => s.as_str(),
                None => "",
            };
            return Err(MageckError::Config(format!(
                "Duplicate control sample selection: sample '{label}' (index {idx}) was selected more than once"
            )));
        }
        control_indices.push(idx);
    }

    for &idx in &treatment_indices {
        if seen_control.contains(&idx) {
            let label = match sample_labels.get(idx) {
                Some(s) => s.as_str(),
                None => "",
            };
            return Err(MageckError::Config(format!(
                "Sample '{label}' (index {idx}) cannot be selected in both treatment and control groups"
            )));
        }
    }

    Ok((treatment_indices, control_indices))
}

/// Parses an sgRNA count table from a reader implementing [`BufRead`].
///
/// Enforces:
/// - Exact TSV format (tab-separated `\t`).
/// - Non-empty header with column 1 = `sgRNA` and column 2 = `Gene`.
/// - At least one sample column after `Gene`.
/// - Unique and non-empty sample column labels.
/// - Deterministic resolution of treatment and control selections.
/// - Exact field count per data row matching header.
/// - Non-empty sgRNA and Gene identifiers on every data row.
/// - Unique sgRNA identifiers.
/// - Finite, non-negative `f64` count values (rejects NaN, +inf, -inf, negative).
/// - Preserves guide row order and sample column order.
/// - Rejects empty data bodies (must have at least one data row).
/// - Constructs [`CountMatrix`] via its validated public API.
pub fn parse_count_table<R: BufRead, S1: AsRef<str>, S2: AsRef<str>>(
    reader: R,
    treatment_ids: &[S1],
    control_ids: &[S2],
) -> Result<CountTable, MageckError> {
    let parsed = crate::data::ParsedCountTable::from_reader(reader)?;
    let (treatment_indices, control_indices) =
        resolve_sample_indices(parsed.sample_labels(), treatment_ids, control_ids)?;
    let (sgrna_ids, genes, sample_labels, counts) = parsed.into_parts();
    Ok(CountTable {
        sgrna_ids,
        genes,
        sample_labels,
        counts,
        treatment_indices,
        control_indices,
    })
}

/// Parses an sgRNA count table from a file on disk.
///
/// Returns [`MageckError::Io`] if the file cannot be opened or read.
pub fn parse_count_table_file<P: AsRef<Path>, S1: AsRef<str>, S2: AsRef<str>>(
    path: P,
    treatment_ids: &[S1],
    control_ids: &[S2],
) -> Result<CountTable, MageckError> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref).map_err(|err| {
        MageckError::Io(io::Error::new(
            err.kind(),
            format!(
                "Failed to open count table file '{}': {err}",
                path_ref.display()
            ),
        ))
    })?;
    let reader = BufReader::new(file);
    parse_count_table(reader, treatment_ids, control_ids)
}
