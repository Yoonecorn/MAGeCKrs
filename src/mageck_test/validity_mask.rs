//! Comparison-scoped guide validity derived from normalized counts.

#![allow(dead_code)] // Foundation is intentionally disconnected until the later RRA integration task.

use crate::cli::RemoveZeroMode;
use crate::data::count_matrix::CountMatrix;
use crate::error::MageckError;
use crate::mageck_test::CountTable;

/// Immutable, row-aligned validity state for one comparison.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ComparisonValidity {
    control_means: Vec<f64>,
    treatment_means: Vec<f64>,
    mask: Vec<bool>,
    n_valid: usize,
    mode: RemoveZeroMode,
    threshold: f64,
}

impl ComparisonValidity {
    pub(crate) fn control_means(&self) -> &[f64] {
        &self.control_means
    }

    pub(crate) fn treatment_means(&self) -> &[f64] {
        &self.treatment_means
    }

    pub(crate) fn mask(&self) -> &[bool] {
        &self.mask
    }

    pub(crate) fn n_valid(&self) -> usize {
        self.n_valid
    }

    pub(crate) fn mode(&self) -> RemoveZeroMode {
        self.mode
    }

    pub(crate) fn threshold(&self) -> f64 {
        self.threshold
    }
}

fn config_error(comparison: &str, invariant: &str) -> MageckError {
    MageckError::Config(format!(
        "category=configuration stage=validity-mask comparison={comparison:?} invariant={invariant}"
    ))
}

fn input_error(comparison: &str, invariant: &str) -> MageckError {
    MageckError::Config(format!(
        "category=input-data stage=validity-mask comparison={comparison:?} invariant={invariant}"
    ))
}

fn math_error(comparison: &str, category: &str, invariant: &str) -> MageckError {
    MageckError::Math(format!(
        "category={category} stage=validity-mask comparison={comparison:?} invariant={invariant}"
    ))
}

fn validate_indices(
    indices: &[usize],
    group: &str,
    n_samples: usize,
    comparison: &str,
) -> Result<(), MageckError> {
    if indices.is_empty() {
        return Err(config_error(
            comparison,
            &format!("{group} group is nonempty"),
        ));
    }
    for (position, &index) in indices.iter().enumerate() {
        if index >= n_samples {
            return Err(config_error(
                comparison,
                &format!("{group} sample index {index} is within 0..{n_samples}"),
            ));
        }
        if indices[..position].contains(&index) {
            return Err(config_error(
                comparison,
                &format!("{group} sample index {index} is unique"),
            ));
        }
    }
    Ok(())
}

fn selected_mean(
    normalized_counts: &CountMatrix,
    indices: &[usize],
    guide: usize,
    group: &str,
    comparison: &str,
) -> Result<f64, MageckError> {
    let mut sum = 0.0;
    for &sample in indices {
        let value = normalized_counts.get(sample, guide).ok_or_else(|| {
            input_error(
                comparison,
                &format!("{group} sample index {sample} and guide index {guide} are addressable"),
            )
        })?;
        if !value.is_finite() || value < 0.0 {
            return Err(input_error(
                comparison,
                &format!(
                    "{group} normalized value at sample index {sample}, guide index {guide} is finite and non-negative (got {value})"
                ),
            ));
        }
        sum += value;
        if !sum.is_finite() {
            return Err(math_error(
                comparison,
                "mathematical-computation",
                &format!("{group} mean accumulation for guide index {guide} is finite"),
            ));
        }
    }
    let mean = sum / indices.len() as f64;
    if !mean.is_finite() {
        return Err(math_error(
            comparison,
            "mathematical-computation",
            &format!("{group} mean for guide index {guide} is finite"),
        ));
    }
    Ok(mean)
}

/// Builds one immutable, physical-row-aligned comparison validity mask.
pub(crate) fn build_comparison_validity(
    normalized_counts: &CountMatrix,
    table: &CountTable,
    mode: RemoveZeroMode,
    threshold: f64,
    comparison: &str,
) -> Result<ComparisonValidity, MageckError> {
    if comparison.is_empty() {
        return Err(config_error(
            comparison,
            "comparison identifier is nonempty",
        ));
    }
    if !threshold.is_finite() || threshold < 0.0 {
        return Err(config_error(
            comparison,
            "threshold is finite and non-negative",
        ));
    }
    if normalized_counts.n_guides() != table.n_guides() {
        return Err(input_error(
            comparison,
            "normalized matrix guide dimension matches count table",
        ));
    }
    if normalized_counts.n_samples() != table.n_samples() {
        return Err(input_error(
            comparison,
            "normalized matrix sample dimension matches count table",
        ));
    }

    let control_indices = table.control_indices();
    let treatment_indices = table.treatment_indices();
    validate_indices(control_indices, "control", table.n_samples(), comparison)?;
    validate_indices(
        treatment_indices,
        "treatment",
        table.n_samples(),
        comparison,
    )?;
    if control_indices
        .iter()
        .any(|index| treatment_indices.contains(index))
    {
        return Err(config_error(
            comparison,
            "control and treatment sample selections do not overlap",
        ));
    }

    let mut control_means = Vec::with_capacity(table.n_guides());
    let mut treatment_means = Vec::with_capacity(table.n_guides());
    let mut mask = Vec::with_capacity(table.n_guides());
    for guide in 0..table.n_guides() {
        let control_mean = selected_mean(
            normalized_counts,
            control_indices,
            guide,
            "control",
            comparison,
        )?;
        let treatment_mean = selected_mean(
            normalized_counts,
            treatment_indices,
            guide,
            "treatment",
            comparison,
        )?;
        let control_pass = control_mean > threshold;
        let treatment_pass = treatment_mean > threshold;
        let valid = match mode {
            RemoveZeroMode::None => true,
            RemoveZeroMode::Control => control_pass,
            RemoveZeroMode::Treatment => treatment_pass,
            RemoveZeroMode::Any => control_pass && treatment_pass,
            RemoveZeroMode::Both => control_pass || treatment_pass,
        };
        control_means.push(control_mean);
        treatment_means.push(treatment_mean);
        mask.push(valid);
    }

    if control_means.len() != table.n_guides()
        || treatment_means.len() != table.n_guides()
        || mask.len() != table.n_guides()
    {
        return Err(math_error(
            comparison,
            "internal-invariant",
            "mean and mask vectors align with physical guide rows",
        ));
    }
    let n_valid = mask.iter().filter(|&&is_valid| is_valid).count();
    if n_valid == 0 {
        return Err(math_error(
            comparison,
            "analysis-state",
            &format!(
                "N_valid > 0 (mode={mode}, threshold={threshold}, finalized valid-guide universe)"
            ),
        ));
    }

    Ok(ComparisonValidity {
        control_means,
        treatment_means,
        mask,
        n_valid,
        mode,
        threshold,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table() -> CountTable {
        CountTable::from_parts(
            vec!["g1".into(), "g2".into(), "g3".into(), "g4".into()],
            vec!["A".into(), "B".into(), "C".into(), "D".into()],
            vec!["c1".into(), "c2".into(), "t1".into(), "t2".into()],
            CountMatrix::from_vec(vec![1.0; 16], 4, 4).unwrap(),
            vec![2, 3],
            vec![0, 1],
        )
        .unwrap()
    }

    fn normalized(values: Vec<f64>) -> CountMatrix {
        CountMatrix::from_vec(values, 4, 4).unwrap()
    }

    fn fixture() -> CountMatrix {
        // Control means [2, 0, 2, 0]; treatment means [2, 2, 0, 0].
        normalized(vec![
            2.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0,
        ])
    }

    #[test]
    fn none_keeps_every_finite_physical_row() {
        let result =
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::None, 0.0, "cmp")
                .unwrap();
        assert_eq!(result.mask(), [true, true, true, true]);
    }

    #[test]
    fn control_requires_strict_control_pass() {
        assert_eq!(
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Control, 0.0, "cmp")
                .unwrap()
                .mask(),
            [true, false, true, false]
        );
    }

    #[test]
    fn treatment_requires_strict_treatment_pass() {
        assert_eq!(
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Treatment, 0.0, "cmp")
                .unwrap()
                .mask(),
            [true, true, false, false]
        );
    }

    #[test]
    fn any_requires_both_groups() {
        assert_eq!(
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Any, 0.0, "cmp")
                .unwrap()
                .mask(),
            [true, false, false, false]
        );
    }

    #[test]
    fn both_requires_either_group() {
        assert_eq!(
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Both, 0.0, "cmp")
                .unwrap()
                .mask(),
            [true, true, true, false]
        );
    }

    #[test]
    fn strict_threshold_rejects_equality() {
        let counts = normalized(vec![
            0.5, 1.0, 1.5, 2.0, 0.5, 1.0, 1.5, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
        ]);
        assert_eq!(
            build_comparison_validity(&counts, &table(), RemoveZeroMode::Control, 1.0, "cmp")
                .unwrap()
                .mask(),
            [false, false, true, true]
        );
    }

    #[test]
    fn negative_zero_is_finite_and_evaluated_by_threshold() {
        let mut counts = fixture();
        *counts.get_mut(0, 0).unwrap() = -0.0;
        assert!(
            build_comparison_validity(&counts, &table(), RemoveZeroMode::None, 0.0, "cmp").is_ok()
        );
    }

    #[test]
    fn comparison_specific_group_selection_changes_means_and_mask() {
        let alternate = CountTable::from_parts(
            vec!["g1".into(), "g2".into(), "g3".into(), "g4".into()],
            vec!["A".into(), "B".into(), "C".into(), "D".into()],
            vec!["c1".into(), "c2".into(), "t1".into(), "t2".into()],
            CountMatrix::from_vec(vec![1.0; 16], 4, 4).unwrap(),
            vec![0, 1],
            vec![2, 3],
        )
        .unwrap();
        assert_ne!(
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Control, 0.0, "first")
                .unwrap()
                .mask(),
            build_comparison_validity(
                &fixture(),
                &alternate,
                RemoveZeroMode::Control,
                0.0,
                "second"
            )
            .unwrap()
            .mask()
        );
    }

    #[test]
    fn physical_row_count_and_order_are_preserved() {
        let result =
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Both, 0.0, "cmp")
                .unwrap();
        assert_eq!(result.mask().len(), table().n_guides());
        assert_eq!(result.control_means(), [2.0, 0.0, 2.0, 0.0]);
    }

    #[test]
    fn mean_and_mask_vectors_align_with_table_rows() {
        let result =
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Both, 0.0, "cmp")
                .unwrap();
        assert_eq!(result.control_means().len(), result.treatment_means().len());
        assert_eq!(result.mask().len(), result.control_means().len());
    }

    #[test]
    fn n_valid_equals_true_entry_count() {
        let result =
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Both, 0.0, "cmp")
                .unwrap();
        assert_eq!(
            result.n_valid(),
            result.mask().iter().filter(|&&value| value).count()
        );
    }

    #[test]
    fn repeated_construction_is_identical() {
        let first =
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Both, 0.0, "cmp")
                .unwrap();
        let second =
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Both, 0.0, "cmp")
                .unwrap();
        assert_eq!(first, second);
        assert_eq!(first.mode(), RemoveZeroMode::Both);
        assert_eq!(first.threshold(), 0.0);
    }

    #[test]
    fn empty_comparison_identifier_is_configuration_error() {
        assert!(matches!(
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::None, 0.0, ""),
            Err(MageckError::Config(_))
        ));
    }

    #[test]
    fn missing_control_group_is_configuration_error() {
        let error = CountTable::from_parts(
            vec!["g".into()],
            vec!["G".into()],
            vec!["c".into()],
            CountMatrix::from_vec(vec![1.0], 1, 1).unwrap(),
            vec![0],
            vec![],
        )
        .unwrap_err();
        assert!(matches!(error, MageckError::Config(_)));
    }

    #[test]
    fn missing_treatment_group_is_configuration_error() {
        let error = CountTable::from_parts(
            vec!["g".into()],
            vec!["G".into()],
            vec!["c".into()],
            CountMatrix::from_vec(vec![1.0], 1, 1).unwrap(),
            vec![],
            vec![0],
        )
        .unwrap_err();
        assert!(matches!(error, MageckError::Config(_)));
    }

    #[test]
    fn invalid_duplicate_overlapping_or_out_of_range_indices_are_configuration_errors() {
        for (treatment, control) in [
            (vec![0, 0], vec![1]),
            (vec![0], vec![0]),
            (vec![2], vec![0]),
        ] {
            let result = CountTable::from_parts(
                vec!["g".into()],
                vec!["G".into()],
                vec!["c".into(), "t".into()],
                CountMatrix::from_vec(vec![1.0, 1.0], 2, 1).unwrap(),
                treatment,
                control,
            );
            assert!(matches!(result, Err(MageckError::Config(_))));
        }
    }

    #[test]
    fn matrix_guide_dimension_mismatch_is_input_data_error() {
        assert!(matches!(
            build_comparison_validity(
                &CountMatrix::from_vec(vec![1.0; 12], 4, 3).unwrap(),
                &table(),
                RemoveZeroMode::None,
                0.0,
                "cmp"
            ),
            Err(MageckError::Config(_))
        ));
    }

    #[test]
    fn matrix_sample_dimension_mismatch_is_input_data_error() {
        assert!(matches!(
            build_comparison_validity(
                &CountMatrix::from_vec(vec![1.0; 12], 3, 4).unwrap(),
                &table(),
                RemoveZeroMode::None,
                0.0,
                "cmp"
            ),
            Err(MageckError::Config(_))
        ));
    }

    #[test]
    fn nonfinite_control_input_is_input_data_error() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut counts = fixture();
            *counts.get_mut(0, 0).unwrap() = value;
            assert!(matches!(
                build_comparison_validity(&counts, &table(), RemoveZeroMode::None, 0.0, "cmp"),
                Err(MageckError::Config(_))
            ));
        }
    }

    #[test]
    fn nonfinite_treatment_input_is_input_data_error() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut counts = fixture();
            *counts.get_mut(2, 0).unwrap() = value;
            assert!(matches!(
                build_comparison_validity(&counts, &table(), RemoveZeroMode::None, 0.0, "cmp"),
                Err(MageckError::Config(_))
            ));
        }
    }

    #[test]
    fn negative_normalized_input_is_input_data_error() {
        let mut counts = fixture();
        *counts.get_mut(0, 0).unwrap() = -1.0;
        assert!(matches!(
            build_comparison_validity(&counts, &table(), RemoveZeroMode::None, 0.0, "cmp"),
            Err(MageckError::Config(_))
        ));
    }

    #[test]
    fn finite_accumulation_overflow_is_mathematical_error() {
        let counts = normalized(vec![
            f64::MAX,
            1.0,
            1.0,
            1.0,
            f64::MAX,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
        ]);
        assert!(matches!(
            build_comparison_validity(&counts, &table(), RemoveZeroMode::None, 0.0, "cmp"),
            Err(MageckError::Math(_))
        ));
    }

    #[test]
    fn globally_empty_none_mode_is_unreachable_for_nonempty_valid_table() {
        assert!(
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::None, 99.0, "cmp")
                .is_ok()
        );
    }

    #[test]
    fn globally_empty_filtered_mask_is_analysis_state_math_error() {
        let error =
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Both, 2.0, "cmp")
                .unwrap_err();
        match error {
            MageckError::Math(message) => assert!(
                message.contains("category=analysis-state") && message.contains("N_valid > 0")
            ),
            other => panic!("expected Math error, got {other:?}"),
        }
    }

    #[test]
    fn error_context_names_category_stage_comparison_and_invariant() {
        let error = build_comparison_validity(
            &fixture(),
            &table(),
            RemoveZeroMode::Both,
            2.0,
            "comparison-a",
        )
        .unwrap_err();
        match error {
            MageckError::Math(message) => {
                assert!(message.contains("category=analysis-state"));
                assert!(message.contains("stage=validity-mask"));
                assert!(message.contains("comparison=\"comparison-a\""));
                assert!(message.contains("invariant="));
            }
            other => panic!("expected Math error, got {other:?}"),
        }
    }

    #[test]
    fn returned_type_contains_no_rank_score_alpha_or_random_state() {
        let result =
            build_comparison_validity(&fixture(), &table(), RemoveZeroMode::Both, 0.0, "cmp")
                .unwrap();
        assert!(!format!("{result:?}").contains("rank"));
        assert!(!format!("{result:?}").contains("score"));
        assert!(!format!("{result:?}").contains("alpha"));
    }
}
