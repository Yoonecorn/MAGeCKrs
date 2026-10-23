//! Validated private design matrix for the MLE ingest boundary.

use std::collections::HashSet;
use std::io::BufRead;

use crate::error::MageckError;

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct DesignMatrix {
    pub(crate) sample_labels: Vec<String>,
    pub(crate) condition_labels: Vec<String>,
    pub(crate) matrix: Vec<f64>,
}

impl DesignMatrix {
    pub(crate) fn from_reader<R: BufRead>(reader: R) -> Result<Self, MageckError> {
        let mut lines = reader.lines();
        let header = lines
            .next()
            .ok_or_else(|| {
                MageckError::DesignMatrix("Design matrix is empty: missing header row".into())
            })?
            .map_err(MageckError::Io)?;
        if header.is_empty() {
            return Err(MageckError::DesignMatrix(
                "Line 1: design matrix header is empty".into(),
            ));
        }
        let fields: Vec<&str> = header.split('\t').collect();
        if fields.len() < 3 {
            return Err(MageckError::DesignMatrix(format!(
                "Line 1: design matrix requires at least 3 TSV fields (sample, intercept, and a fitted condition), got {}",
                fields.len()
            )));
        }
        for (index, name) in ["sample", "intercept"].iter().enumerate() {
            if fields[index].is_empty() {
                return Err(MageckError::DesignMatrix(format!(
                    "Line 1, column {}: {name} header label is empty",
                    index + 1
                )));
            }
        }
        let condition_labels: Vec<String> = fields[2..]
            .iter()
            .map(|value| (*value).to_owned())
            .collect();
        let mut seen_conditions = HashSet::new();
        for (index, label) in condition_labels.iter().enumerate() {
            if label.is_empty() {
                return Err(MageckError::DesignMatrix(format!(
                    "Line 1, column {}: condition label is empty",
                    index + 3
                )));
            }
            if !seen_conditions.insert(label.clone()) {
                return Err(MageckError::DesignMatrix(format!(
                    "Line 1, column {}: duplicate condition label '{label}'",
                    index + 3
                )));
            }
        }
        let mut sample_labels = Vec::new();
        let mut matrix = Vec::new();
        let mut seen_samples = HashSet::new();
        for (line_offset, line) in lines.enumerate() {
            let line_number = line_offset + 2;
            let line = line.map_err(MageckError::Io)?;
            let row: Vec<&str> = line.split('\t').collect();
            if row.len() != fields.len() {
                return Err(MageckError::DesignMatrix(format!(
                    "Line {line_number}: field count mismatch (expected {} columns, got {})",
                    fields.len(),
                    row.len()
                )));
            }
            let sample = row[0];
            if sample.is_empty() {
                return Err(MageckError::DesignMatrix(format!(
                    "Line {line_number}, column 1: sample label is empty"
                )));
            }
            if !seen_samples.insert(sample.to_owned()) {
                return Err(MageckError::Config(format!(
                    "Line {line_number}, column 1: duplicate sample row label '{sample}' in design matrix"
                )));
            }
            // Complete each field before visiting the next: syntax, finiteness,
            // canonical zero, then its intercept or first-row baseline invariant.
            for (index, field) in row.iter().enumerate().skip(1) {
                let column = index + 1;
                let value = field.parse::<f64>().map_err(|_| {
                    MageckError::DesignMatrix(format!(
                        "Line {line_number}, column {column}: non-numeric coefficient '{field}'"
                    ))
                })?;
                if !value.is_finite() {
                    return Err(MageckError::DesignMatrix(format!(
                        "Line {line_number}, column {column}: non-finite coefficient '{field}'"
                    )));
                }
                let value = if value == 0.0 { 0.0 } else { value };
                if index == 1 && value != 1.0 {
                    return Err(MageckError::DesignMatrix(format!(
                        "Line {line_number}, column {column}: intercept must evaluate to 1.0, got {value}"
                    )));
                }
                if index > 1 && sample_labels.is_empty() && value != 0.0 {
                    return Err(MageckError::Config(format!(
                        "Line {line_number}, column {column}: row 0 must be a baseline sample with condition '{}' equal to 0.0, got {value}",
                        condition_labels[index - 2]
                    )));
                }
                matrix.push(value);
            }
            sample_labels.push(sample.to_owned());
        }
        if sample_labels.len() < 2 {
            return Err(MageckError::DesignMatrix(format!(
                "Design matrix requires at least 2 sample rows, got {}",
                sample_labels.len()
            )));
        }
        validate_full_column_rank(
            &matrix,
            sample_labels.len(),
            fields.len() - 1,
            &condition_labels,
        )?;
        Ok(Self {
            sample_labels,
            condition_labels,
            matrix,
        })
    }

    pub(crate) fn from_str(input: &str) -> Result<Self, MageckError> {
        Self::from_reader(input.as_bytes())
    }
    pub(crate) fn sample_labels(&self) -> &[String] {
        &self.sample_labels
    }
    pub(crate) fn condition_labels(&self) -> &[String] {
        &self.condition_labels
    }
    pub(crate) fn matrix(&self) -> &[f64] {
        &self.matrix
    }
    pub(crate) fn n_samples(&self) -> usize {
        self.sample_labels.len()
    }
    pub(crate) fn n_conditions(&self) -> usize {
        self.condition_labels.len()
    }
    pub(crate) fn total_columns(&self) -> usize {
        self.n_conditions() + 1
    }
    pub(crate) fn entry(&self, row: usize, column: usize) -> f64 {
        self.matrix[row * self.total_columns() + column]
    }
}

fn validate_full_column_rank(
    matrix: &[f64],
    n_rows: usize,
    n_columns: usize,
    condition_labels: &[String],
) -> Result<(), MageckError> {
    if n_columns > n_rows {
        return Err(MageckError::DesignMatrix(format!(
            "Design matrix is not full column rank: m={n_rows}, n={n_columns}"
        )));
    }

    let frobenius_norm = scaled_l2(matrix);
    let threshold = 64.0 * f64::EPSILON * n_rows.max(n_columns) as f64 * frobenius_norm.max(1.0);
    let mut q_columns: Vec<Vec<f64>> = Vec::with_capacity(n_columns);

    for column in 0..n_columns {
        let mut residual = Vec::with_capacity(n_rows);
        for row in 0..n_rows {
            residual.push(matrix[row * n_columns + column]);
        }

        for _pass in 0..2 {
            for q_column in &q_columns {
                let mut dot = 0.0;
                for row in 0..n_rows {
                    dot += q_column[row] * residual[row];
                }
                for row in 0..n_rows {
                    residual[row] -= dot * q_column[row];
                }
            }
        }

        let residual_norm = scaled_l2(&residual);
        if residual_norm <= threshold {
            let label = if column == 0 {
                "intercept"
            } else {
                &condition_labels[column - 1]
            };
            return Err(MageckError::DesignMatrix(format!(
                "Design matrix is not full column rank: m={n_rows}, n={n_columns}; \
                 dependent column index {column}, label '{label}'; residual norm \
                 {residual_norm:e}; threshold {threshold:e}"
            )));
        }

        let mut q_column = Vec::with_capacity(n_rows);
        for &value in &residual {
            q_column.push(value / residual_norm);
        }
        q_columns.push(q_column);
    }

    Ok(())
}

fn scaled_l2(values: &[f64]) -> f64 {
    let mut scale = 0.0;
    let mut ssq = 1.0;
    for &value in values {
        let absolute = value.abs();
        if absolute != 0.0 {
            if scale < absolute {
                let ratio = scale / absolute;
                ssq = 1.0 + ssq * ratio * ratio;
                scale = absolute;
            } else {
                let ratio = absolute / scale;
                ssq += ratio * ratio;
            }
        }
    }
    if scale == 0.0 {
        0.0
    } else {
        scale * ssq.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::DesignMatrix;
    use crate::error::MageckError;

    fn rank_error(input: &str) -> String {
        match DesignMatrix::from_str(input) {
            Err(MageckError::DesignMatrix(message)) => message,
            result => panic!("expected design-matrix rank error, got {result:?}"),
        }
    }

    #[test]
    fn fixture_04a_multiple_baseline_samples_valid() {
        let design =
            DesignMatrix::from_str("sample\tbaseline\tbeta\nDay0\t1\t0\nDay1\t1\t0\nTreat\t1\t1\n")
                .unwrap();
        assert_eq!(design.entry(1, 1), 0.0);
    }

    #[test]
    fn fixture_04b_row_zero_non_baseline_rejected() {
        assert!(
            DesignMatrix::from_str("sample\tbaseline\tbeta\nTreat\t1\t1\nDay0\t1\t0\n").is_err()
        );
    }

    #[test]
    fn fixture_04c_column_zero_non_one_rejected() {
        assert!(
            DesignMatrix::from_str("sample\tbaseline\tbeta\nDay0\t0\t0\nTreat\t1\t1\n").is_err()
        );
    }

    #[test]
    fn fixture_05b_duplicate_design_matrix_label_error() {
        assert!(
            DesignMatrix::from_str("sample\tbaseline\tbeta\nDay0\t1\t0\nDay0\t1\t1\n").is_err()
        );
    }

    #[test]
    fn rank_oracle_o1_valid_full_rank_design() {
        let design = DesignMatrix::from_str(
            "sample\tbaseline\ta\tb\nD0\t1\t0\t0\nS1\t1\t1\t0\nS2\t1\t0\t1\nS3\t1\t1\t1\n",
        )
        .unwrap();
        assert_eq!(design.total_columns(), 3);
    }

    #[test]
    fn rank_oracles_o2_and_o4_reject_first_dependent_column() {
        let duplicate =
            rank_error("sample\tbaseline\ta\tb\nD0\t1\t0\t0\nS1\t1\t1\t1\nS2\t1\t0\t0\n");
        assert!(duplicate.contains("not full column rank"));
        assert!(duplicate.contains("m=3, n=3"));
        assert!(duplicate.contains("index 2"));
        assert!(duplicate.contains("label 'b'"));
        assert!(duplicate.contains("residual norm"));
        assert!(duplicate.contains("threshold"));

        let linear_combination =
            rank_error("sample\tbaseline\ta\tb\nD0\t1\t0\t0\nS1\t1\t1\t2\nS2\t1\t2\t4\n");
        assert!(linear_combination.contains("index 2"));
        assert!(linear_combination.contains("label 'b'"));
    }

    #[test]
    fn rank_oracle_o3_rejects_zero_condition_column() {
        let message = rank_error("sample\tbaseline\ta\nD0\t1\t0\nS1\t1\t0\nS2\t1\t0\n");
        assert!(message.contains("index 1"));
        assert!(message.contains("label 'a'"));
    }

    #[test]
    fn rank_oracle_o5_rejects_wide_matrix_before_qr() {
        let message = rank_error("sample\tbaseline\ta\tb\nD0\t1\t0\t0\nS1\t1\t1\t0\n");
        assert!(message.contains("not full column rank"));
        assert!(message.contains("m=2, n=3"));
        assert!(!message.contains("residual norm"));
        assert!(!message.contains("dependent column"));
    }

    #[test]
    fn rank_oracles_o6_through_o10_accept_independent_and_scaled_conditions() {
        for input in [
            "sample\tbaseline\ta\nD0\t1\t0\nS1\t1\t1\n",
            "sample\tbaseline\ta\nD0\t1\t0\nS1\t1\t0\nS2\t1\t1e-13\n",
            "sample\tbaseline\ta\nD0\t1\t0\nS1\t1\t0\nS2\t1\t1e-10\n",
            "sample\tbaseline\ta\nD0\t1\t0\nS1\t1\t0\nS2\t1\t1e-7\n",
            "sample\tbaseline\ta\nD0\t1\t0\nS1\t1\t0\nS2\t1\t1e-4\n",
        ] {
            assert!(DesignMatrix::from_str(input).is_ok(), "{input}");
        }
    }

    #[test]
    fn rank_oracles_o11_and_o12_preserve_non_finite_precedence() {
        for value in ["NaN", "inf"] {
            match DesignMatrix::from_str(&format!(
                "sample\tbaseline\ta\nD0\t1\t0\nS1\t1\t{value}\n"
            )) {
                Err(MageckError::DesignMatrix(message)) => {
                    assert!(message.contains("non-finite"));
                    assert!(!message.contains("not full column rank"));
                }
                result => panic!("expected non-finite parser error, got {result:?}"),
            }
        }
    }

    #[test]
    fn rank_diagnostics_and_results_are_deterministic() {
        let valid = "sample\tbaseline\ta\tb\nD0\t1\t0\t0\nS1\t1\t1\t0\nS2\t1\t0\t1\nS3\t1\t1\t1\n";
        let duplicate = "sample\tbaseline\ta\tb\nD0\t1\t0\t0\nS1\t1\t1\t1\nS2\t1\t0\t0\n";
        let expected = rank_error(duplicate);
        for _ in 0..100 {
            assert!(DesignMatrix::from_str(valid).is_ok());
            assert_eq!(rank_error(duplicate), expected);
        }
    }
    // Exact assertions pin category, text, physical line and one-based TSV token
    // column. Rank assertions elsewhere deliberately pin components only.
    fn assert_exact(input: &str, config: bool, expected: &str) {
        let message = match DesignMatrix::from_str(input) {
            Err(MageckError::Config(message)) if config => message,
            Err(MageckError::DesignMatrix(message)) if !config => message,
            result => panic!("wrong error category for {input:?}: {result:?}"),
        };
        assert_eq!(message, expected, "input: {input:?}");
    }

    #[test]
    fn dm1b_exact_header_and_structural_diagnostics() {
        for (input, expected) in [
            ("", "Design matrix is empty: missing header row"),
            ("\n", "Line 1: design matrix header is empty"),
            ("sample\nDay0\tnot_a_num\n", "Line 1: design matrix requires at least 3 TSV fields (sample, intercept, and a fitted condition), got 1"),
            ("sample\tintercept\n", "Line 1: design matrix requires at least 3 TSV fields (sample, intercept, and a fitted condition), got 2"),
            ("\tintercept\ta\n", "Line 1, column 1: sample header label is empty"),
            ("sample\t\ta\n", "Line 1, column 2: intercept header label is empty"),
            ("sample\tintercept\t\n", "Line 1, column 3: condition label is empty"),
            ("sample\tintercept\ta\ta\n", "Line 1, column 4: duplicate condition label 'a'"),
            ("sample\tintercept\ta\nD0\t1\t0\nD0\tbad\n", "Line 3: field count mismatch (expected 3 columns, got 2)"),
            ("sample\tintercept\ta\n\tNaN\tbad\n", "Line 2, column 1: sample label is empty"),
            ("sample\tintercept\ta\n", "Design matrix requires at least 2 sample rows, got 0"),
            ("sample\tintercept\ta\nD0\t1\t0\n", "Design matrix requires at least 2 sample rows, got 1"),
        ] { assert_exact(input, false, expected); }
    }

    #[test]
    fn dm1b_exact_adversarial_field_precedence() {
        for token in ["bad", "NaN"] {
            assert_exact(
                &format!("sample\ti\ta\nD0\t1\t0\nD0\t1\t{token}\n"),
                true,
                "Line 3, column 1: duplicate sample row label 'D0' in design matrix",
            );
        }
        for (input, config, expected) in [
            ("sample\ti\ta\tb\nD0\t1\t0\t0\nS1\t1\tbad_float\tNaN\n", false, "Line 3, column 3: non-numeric coefficient 'bad_float'"),
            ("sample\ti\ta\nD0\t1\t0\nS1\tNaN\tbad_float\n", false, "Line 3, column 2: non-finite coefficient 'NaN'"),
            ("sample\ti\ta\nD0\t1\t0\nS1\t2\tNaN\n", false, "Line 3, column 2: intercept must evaluate to 1.0, got 2"),
            ("sample\ti\ta\nD0\t1\t0\nS1\t2\tbad_float\n", false, "Line 3, column 2: intercept must evaluate to 1.0, got 2"),
            ("sample\ti\ta\nD0\t2\t1\n", false, "Line 2, column 2: intercept must evaluate to 1.0, got 2"),
            ("sample\ti\ta\nD0\t1\t1\n", true, "Line 2, column 3: row 0 must be a baseline sample with condition 'a' equal to 0.0, got 1"),
            ("sample\ti\ta\nD0\t1\tbad_float\n", false, "Line 2, column 3: non-numeric coefficient 'bad_float'"),
            ("sample\ti\ta\tb\nD0\t1\t1\tbad_float\n", true, "Line 2, column 3: row 0 must be a baseline sample with condition 'a' equal to 0.0, got 1"),
            ("sample\ti\ta\tb\nD0\t1\t1\tNaN\n", true, "Line 2, column 3: row 0 must be a baseline sample with condition 'a' equal to 0.0, got 1"),
            ("sample\ti\ta\nD0\t 1\t0\n", false, "Line 2, column 2: non-numeric coefficient ' 1'"),
            ("sample\ti\ta\nD0\t1\t0\nS1\t1\t1e999\n", false, "Line 3, column 3: non-finite coefficient '1e999'"),
        ] { assert_exact(input, config, expected); }
        for token in ["inf", "-inf", "NaN"] {
            assert_exact(
                &format!("sample\ti\ta\nD0\t1\t0\nS1\t1\t{token}\n"),
                false,
                &format!("Line 3, column 3: non-finite coefficient '{token}'"),
            );
        }
    }

    #[test]
    fn dm1b_signed_zero_bits_and_literal_header_grammar() {
        for zero in ["-0", "-0.0", "-0e10", "-1e-999", "+0"] {
            let input = format!(" sample \ta\ta\t a \n D0 \t1e0\t{zero}\t{zero}\nS1\t1.0\t-1\t{zero}\nS2\t1\t{zero}\t1\n");
            let design = DesignMatrix::from_str(&input).unwrap();
            assert_eq!(design.sample_labels()[0], " D0 ");
            assert_eq!(design.condition_labels(), &["a", " a "]);
            for index in [1, 2, 5, 7] {
                assert_eq!(design.matrix()[index].to_bits(), 0.0_f64.to_bits());
            }
            assert_eq!(design.entry(1, 1).to_bits(), (-1.0_f64).to_bits());
            assert_exact(
                &format!("s\ti\ta\nD0\t{zero}\t0\n"),
                false,
                "Line 2, column 2: intercept must evaluate to 1.0, got 0",
            );
        }
    }

    #[test]
    fn dm1b_reader_failures_remain_io_and_obey_fail_fast() {
        use std::io::{self, BufReader, Cursor, Read};
        struct Failure;
        impl Read for Failure {
            fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
                Err(io::Error::other("injected reader failure"))
            }
        }
        for prefix in ["", "sample\ti\ta\n", "sample\ti\ta\nD0\t1\t0\n"] {
            let reader = BufReader::new(Cursor::new(prefix.as_bytes()).chain(Failure));
            match DesignMatrix::from_reader(reader) {
                Err(MageckError::Io(error)) => {
                    assert_eq!(error.to_string(), "injected reader failure")
                }
                result => panic!("expected Io before EOF/rank: {result:?}"),
            }
        }
        let reader = BufReader::new(Cursor::new(b"sample\ti\ta\nD0\t2\tbad\n").chain(Failure));
        match DesignMatrix::from_reader(reader) {
            Err(MageckError::DesignMatrix(message)) => assert_eq!(
                message,
                "Line 2, column 2: intercept must evaluate to 1.0, got 2"
            ),
            result => panic!("expected earlier intercept failure: {result:?}"),
        }
        assert!(matches!(
            DesignMatrix::from_reader(&b"sample\ti\ta\nD0\t1\t\xff\n"[..]),
            Err(MageckError::Io(_))
        ));
    }
}
