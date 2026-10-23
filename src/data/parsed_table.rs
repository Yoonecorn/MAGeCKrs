//! Role-neutral parsed count-table representation and streaming TSV parser.
//!
//! Implements intrinsic count-table validation without assigning sample roles.

use std::collections::{HashMap, HashSet};
use std::io::BufRead;

use crate::data::count_matrix::CountMatrix;
use crate::error::MageckError;

/// Parsed sgRNA count table with physical source-column ordering preserved.
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedCountTable {
    sgrna_ids: Vec<String>,
    genes: Vec<String>,
    sample_labels: Vec<String>,
    counts: CountMatrix,
}

impl ParsedCountTable {
    /// Parses and intrinsically validates an sgRNA count table from a streaming reader.
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self, MageckError> {
        let mut lines = reader.lines();
        let mut line_num = 0usize;
        let header_line = match lines.next() {
            Some(res) => {
                line_num += 1;
                res.map_err(MageckError::Io)?
            }
            None => {
                return Err(MageckError::Config(
                    "Count table is empty: missing header row".to_string(),
                ));
            }
        };

        if header_line.is_empty() {
            return Err(MageckError::Config(
                "Count table header row is empty".to_string(),
            ));
        }

        let header_cols: Vec<&str> = header_line.split('\t').collect();
        if header_cols.len() < 2 {
            return Err(MageckError::Config(
                "Count table header requires at least 'sgRNA' and 'Gene' columns".to_string(),
            ));
        }

        let col0 = header_cols.first().copied().unwrap_or("");
        if col0 != "sgRNA" {
            return Err(MageckError::Config(format!(
                "Count table header column 1 must be 'sgRNA', got '{col0}'"
            )));
        }

        let col1 = header_cols.get(1).copied().unwrap_or("");
        if col1 != "Gene" {
            return Err(MageckError::Config(format!(
                "Count table header column 2 must be 'Gene', got '{col1}'"
            )));
        }

        if header_cols.len() < 3 {
            return Err(MageckError::Config(
                "Count table contains no sample columns (expected at least one sample column after 'sgRNA' and 'Gene')"
                    .to_string(),
            ));
        }

        let n_cols = header_cols.len();
        let sample_labels: Vec<String> = header_cols[2..].iter().map(|&s| s.to_string()).collect();
        let n_samples = sample_labels.len();
        let mut seen_samples = HashSet::with_capacity(n_samples);
        for label in &sample_labels {
            if label.is_empty() {
                return Err(MageckError::Config(
                    "Count table header contains an empty sample column label".to_string(),
                ));
            }
            if !seen_samples.insert(label.clone()) {
                return Err(MageckError::Config(format!(
                    "Duplicate sample column label in count table header: '{label}'"
                )));
            }
        }

        let mut sgrna_ids = Vec::new();
        let mut genes = Vec::new();
        let mut seen_sgrnas = HashMap::new();
        let mut counts_by_sample = vec![Vec::new(); n_samples];
        for line_res in lines {
            line_num += 1;
            let line = line_res.map_err(MageckError::Io)?;
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != n_cols {
                return Err(MageckError::Config(format!(
                    "Line {line_num}: field count mismatch (expected {n_cols} columns to match header, got {})",
                    fields.len()
                )));
            }

            let sgrna = fields.first().copied().unwrap_or("");
            if sgrna.is_empty() {
                return Err(MageckError::Config(format!(
                    "Line {line_num}: sgRNA identifier is empty"
                )));
            }
            if let Some(&first_line) = seen_sgrnas.get(sgrna) {
                return Err(MageckError::Config(format!(
                    "Line {line_num}: duplicate sgRNA identifier '{sgrna}' (first seen on line {first_line})"
                )));
            }
            seen_sgrnas.insert(sgrna.to_string(), line_num);

            let gene = fields.get(1).copied().unwrap_or("");
            if gene.is_empty() {
                return Err(MageckError::Config(format!(
                    "Line {line_num}: gene identifier is empty for sgRNA '{sgrna}'"
                )));
            }

            for (sample_index, sample_slot) in counts_by_sample.iter_mut().enumerate() {
                let field_val = fields.get(2 + sample_index).copied().unwrap_or("");
                let val = match field_val.parse::<f64>() {
                    Ok(value) => value,
                    Err(_) => {
                        return Err(MageckError::Config(format!(
                            "Line {line_num}: sample '{}' has non-numeric count '{field_val}'",
                            sample_labels[sample_index]
                        )));
                    }
                };
                if !val.is_finite() {
                    return Err(MageckError::Config(format!(
                        "Line {line_num}: sample '{}' has non-finite count '{field_val}'",
                        sample_labels[sample_index]
                    )));
                }
                if val < 0.0 {
                    return Err(MageckError::Config(format!(
                        "Line {line_num}: sample '{}' has negative count '{field_val}'",
                        sample_labels[sample_index]
                    )));
                }
                sample_slot.push(if val == 0.0 { 0.0 } else { val });
            }
            sgrna_ids.push(sgrna.to_string());
            genes.push(gene.to_string());
        }

        let n_guides = sgrna_ids.len();
        if n_guides == 0 {
            return Err(MageckError::Config(
                "Count table contains no data rows".to_string(),
            ));
        }
        let total = n_samples.checked_mul(n_guides).ok_or_else(|| {
            MageckError::Config(
                "Count table dimensions overflow usize (n_samples * n_guides)".to_string(),
            )
        })?;
        let mut flat_data = Vec::with_capacity(total);
        for sample_vec in counts_by_sample {
            flat_data.extend(sample_vec);
        }
        let counts = CountMatrix::from_vec(flat_data, n_samples, n_guides).map_err(|error| {
            MageckError::Config(format!("Failed to construct count matrix: {error}"))
        })?;

        Ok(Self {
            sgrna_ids,
            genes,
            sample_labels,
            counts,
        })
    }

    /// Validated neutral handoff from COUNT. No roles, transformations or normalization.
    /// Visibility is limited to orchestration inside this crate; the TSV parser is unchanged.
    pub(crate) fn from_parts(
        sgrna_ids: Vec<String>,
        genes: Vec<String>,
        sample_labels: Vec<String>,
        mut counts: CountMatrix,
    ) -> Result<Self, MageckError> {
        if sgrna_ids.is_empty()
            || sample_labels.is_empty()
            || genes.len() != sgrna_ids.len()
            || counts.n_guides() != sgrna_ids.len()
            || counts.n_samples() != sample_labels.len()
        {
            return Err(MageckError::Config(
                "Neutral count-table dimensions must agree and contain guides and samples".into(),
            ));
        }
        for (fields, label, unique) in [
            (&sgrna_ids, "sgRNA", true),
            (&genes, "Gene", false),
            (&sample_labels, "sample", true),
        ] {
            let mut seen = HashSet::new();
            for field in fields {
                if field.is_empty() || field.contains(['\t', '\n']) {
                    return Err(MageckError::Config(format!(
                        "Invalid neutral {label} identifier {field:?}"
                    )));
                }
                if unique && !seen.insert(field) {
                    return Err(MageckError::Config(format!(
                        "Duplicate neutral {label} identifier {field:?}"
                    )));
                }
            }
        }
        for value in counts.as_mut_slice() {
            if !value.is_finite() || *value < 0.0 {
                return Err(MageckError::Config(
                    "Neutral counts must be finite and non-negative".into(),
                ));
            }
            // The existing neutral TSV parser canonicalizes signed zero too.
            if *value == 0.0 {
                *value = 0.0;
            }
        }
        Ok(Self {
            sgrna_ids,
            genes,
            sample_labels,
            counts,
        })
    }

    /// Returns the ordered sgRNA identifiers.
    #[inline]
    pub fn sgrna_ids(&self) -> &[String] {
        &self.sgrna_ids
    }

    /// Returns the ordered gene identifiers aligned to sgRNAs.
    #[inline]
    pub fn genes(&self) -> &[String] {
        &self.genes
    }

    /// Returns source-order sample labels.
    #[inline]
    pub fn sample_labels(&self) -> &[String] {
        &self.sample_labels
    }

    /// Returns the sample-major count matrix.
    #[inline]
    pub fn counts(&self) -> &CountMatrix {
        &self.counts
    }

    /// Returns the number of guide rows.
    #[inline]
    pub fn n_guides(&self) -> usize {
        self.sgrna_ids.len()
    }

    /// Returns the number of sample columns.
    #[inline]
    pub fn n_samples(&self) -> usize {
        self.sample_labels.len()
    }

    /// Returns whether the table contains no guide rows.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sgrna_ids.is_empty()
    }

    /// Returns a source-order sample index by exact label.
    pub fn sample_index(&self, label: &str) -> Option<usize> {
        self.sample_labels.iter().position(|sample| sample == label)
    }

    /// Consumes the table without cloning its validated components.
    pub fn into_parts(self) -> (Vec<String>, Vec<String>, Vec<String>, CountMatrix) {
        (self.sgrna_ids, self.genes, self.sample_labels, self.counts)
    }
}

#[cfg(test)]
mod campaign_constructor_tests {
    use super::*;
    fn construct(
        ids: &[&str],
        genes: &[&str],
        labels: &[&str],
        values: Vec<f64>,
    ) -> Result<ParsedCountTable, MageckError> {
        let matrix = CountMatrix::from_vec(values, 2, 2).unwrap();
        ParsedCountTable::from_parts(
            ids.iter().map(|s| s.to_string()).collect(),
            genes.iter().map(|s| s.to_string()).collect(),
            labels.iter().map(|s| s.to_string()).collect(),
            matrix,
        )
    }
    #[test]
    fn preserves_order_raw_fractional_values_and_matches_tsv() {
        let table = construct(
            &["second", "first"],
            &["B", "A"],
            &["z", "a"],
            vec![4.5, -0.0, 2.0, 3.25],
        )
        .unwrap();
        let parsed = ParsedCountTable::from_reader(
            "sgRNA\tGene\tz\ta\nsecond\tB\t4.5\t2\nfirst\tA\t-0\t3.25\n".as_bytes(),
        )
        .unwrap();
        assert_eq!(table, parsed);
        assert_eq!(table.counts().as_slice()[1].to_bits(), 0.0f64.to_bits());
        assert_eq!(table.sgrna_ids(), &["second", "first"]);
        // Shared gene IDs and space-containing nominal identifiers remain valid.
        assert!(construct(&[" a ", "b"], &["A", "A"], &["z", "a"], vec![1.0; 4]).is_ok());
    }
    #[test]
    fn embedded_carriage_return_is_preserved_like_existing_tsv_parser() {
        let table = construct(&["a\rb", "c"], &["A", "B"], &["z", "a"], vec![1.0; 4]).unwrap();
        let parsed = ParsedCountTable::from_reader(
            "sgRNA\tGene\tz\ta\na\rb\tA\t1\t1\nc\tB\t1\t1\n".as_bytes(),
        )
        .unwrap();
        assert_eq!(table, parsed);
    }

    #[test]
    fn rejects_dimensions_empty_duplicate_or_unrepresentable_identifiers() {
        for (ids, genes, labels) in [
            (vec!["a"], vec!["A", "B"], vec!["z", "a"]),
            (vec!["a", "b"], vec!["A"], vec!["z", "a"]),
            (vec!["a", "b"], vec!["A", "B"], vec!["z"]),
            (vec!["a", "a"], vec!["A", "B"], vec!["z", "a"]),
            (vec!["a", "b"], vec!["A", "B"], vec!["z", "z"]),
            (vec!["", "b"], vec!["A", "B"], vec!["z", "a"]),
            (vec!["a", "b"], vec!["", "B"], vec!["z", "a"]),
            (vec!["a", "b"], vec!["A", "B"], vec!["", "a"]),
            (vec!["a\tb", "b"], vec!["A", "B"], vec!["z", "a"]),
        ] {
            assert!(construct(&ids, &genes, &labels, vec![1.0; 4]).is_err());
        }
        assert!(
            ParsedCountTable::from_parts(vec![], vec![], vec![], CountMatrix::new(0, 0)).is_err()
        );
    }
    #[test]
    fn rejects_negative_nonfinite_counts_without_integer_restriction() {
        for value in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(construct(
                &["a", "b"],
                &["A", "B"],
                &["z", "a"],
                vec![value, 1.0, 2.0, 3.0]
            )
            .is_err());
        }
    }
}
