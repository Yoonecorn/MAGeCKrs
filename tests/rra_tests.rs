use std::io::Cursor;

use mageckrs::error::MageckError;
use mageckrs::mageck_test::{
    compute_gene_lo_values, compute_gene_lo_values_from_count_table, parse_count_table,
    GeneLoValueResult,
};
use mageckrs::stats::DirectionalAlpha;

fn genes(names: &[&str]) -> Vec<String> {
    names.iter().map(|name| (*name).to_string()).collect()
}

fn alpha(negative: f64, positive: f64) -> DirectionalAlpha {
    DirectionalAlpha { negative, positive }
}

fn score(
    names: &[&str],
    negative: &[Option<f64>],
    positive: &[Option<f64>],
    thresholds: DirectionalAlpha,
    n_valid: usize,
) -> Vec<GeneLoValueResult> {
    compute_gene_lo_values(&genes(names), negative, positive, thresholds, n_valid).unwrap()
}

fn target_with_fillers(
    target: &str,
    negative: &[Option<f64>],
    positive: &[Option<f64>],
    thresholds: DirectionalAlpha,
    n_valid: usize,
) -> GeneLoValueResult {
    let mut all_genes = vec![target.to_string(); negative.len()];
    let mut all_negative = negative.to_vec();
    let mut all_positive = positive.to_vec();
    let present = negative.iter().filter(|rank| rank.is_some()).count();
    assert!(n_valid >= present);
    for index in present..n_valid {
        all_genes.push(format!("FILLER_{index}"));
        all_negative.push(Some(0.5));
        all_positive.push(Some(0.5));
    }
    compute_gene_lo_values(
        &all_genes,
        &all_negative,
        &all_positive,
        thresholds,
        n_valid,
    )
    .unwrap()
    .remove(0)
}

fn math_message(result: Result<Vec<GeneLoValueResult>, MageckError>) -> String {
    match result {
        Err(MageckError::Math(message)) => message,
        other => panic!("Expected MageckError::Math, got {other:?}"),
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "expected {expected:.15}, got {actual:.15}"
    );
}

// T01 / GS-9: scattered rows retain first-appearance order and physical counts.
#[test]
fn test_grouping_contiguous_and_scattered_genes() {
    let result = score(
        &["A", "B", "A", "C", "B"],
        &[Some(0.1), Some(0.2), Some(0.3), Some(0.4), Some(0.5)],
        &[Some(0.2), Some(0.3), Some(0.4), Some(0.5), Some(0.6)],
        alpha(0.25, 0.25),
        5,
    );
    assert_eq!(
        result
            .iter()
            .map(|r| r.gene_id.as_str())
            .collect::<Vec<_>>(),
        ["A", "B", "C"]
    );
    assert_eq!(
        result.iter().map(|r| r.num_sgrnas).collect::<Vec<_>>(),
        [2, 2, 1]
    );
    assert_eq!(
        result
            .iter()
            .map(|r| r.num_valid_sgrnas)
            .collect::<Vec<_>>(),
        [2, 2, 1]
    );
}

// T02 / GS-2: one valid guide is scored with the gene-local count.
#[test]
fn test_grouping_single_guide_genes() {
    let result = target_with_fillers("SINGLE", &[Some(0.05)], &[Some(0.5)], alpha(0.1, 0.1), 100);
    assert_eq!(result.num_sgrnas, 1);
    assert_eq!(result.num_valid_sgrnas, 1);
    near(result.neg_score, 0.05, 1e-14);
    near(result.pos_score, 0.5, 1e-12);
    assert_eq!((result.neg_goodsgrna, result.pos_goodsgrna), (1, 0));
}

// T03: all physical observations may target one gene.
#[test]
fn test_grouping_all_same_gene() {
    let result = score(
        &["GENE1"; 5],
        &[Some(0.1), Some(0.2), Some(0.3), Some(0.4), Some(0.5)],
        &[Some(0.5); 5],
        alpha(0.25, 0.25),
        5,
    );
    assert_eq!(result.len(), 1);
    assert_eq!((result[0].num_sgrnas, result[0].num_valid_sgrnas), (5, 5));
}

// T04 / GS-ERR-1: either directional slice may disagree with physical rows.
#[test]
fn test_identity_validation_mismatched_lengths() {
    let names = genes(&["A", "B", "C"]);
    for (negative, positive) in [
        (vec![Some(0.1); 2], vec![Some(0.2); 3]),
        (vec![Some(0.1); 3], vec![Some(0.2); 2]),
    ] {
        let message = math_message(compute_gene_lo_values(
            &names,
            &negative,
            &positive,
            alpha(0.2, 0.2),
            2,
        ));
        assert!(message.starts_with("Length mismatch: genes (3)"));
    }
}

// T05 / GS-ERR-5: aligned rank slices replace the old guide-index metadata.
#[test]
fn test_identity_validation_directional_mask_alignment() {
    let message = math_message(compute_gene_lo_values(
        &genes(&["A", "B"]),
        &[Some(0.1), None],
        &[Some(0.2), Some(0.3)],
        alpha(0.2, 0.2),
        1,
    ));
    assert_eq!(
        message,
        "Row 1: directional mask mismatch between negative and positive ranks"
    );
}

// T06 / GS-ERR-2.
#[test]
fn test_identity_validation_empty_inputs() {
    let message = math_message(compute_gene_lo_values(&[], &[], &[], alpha(0.2, 0.2), 0));
    assert_eq!(message, "Input guide slices must not be empty");
}

// T07 / GS-ERR-3: only a truly empty identifier fails.
#[test]
fn test_identity_validation_empty_gene_string() {
    let message = math_message(compute_gene_lo_values(
        &genes(&["A", ""]),
        &[Some(0.2), Some(0.3)],
        &[Some(0.2), Some(0.3)],
        alpha(0.2, 0.2),
        2,
    ));
    assert_eq!(message, "Row 1: gene identifier is empty");
}

// T08 / GS-ERR-6: no NaN or infinity may serve as a structural sentinel.
#[test]
fn test_rank_validation_nan_and_infinity() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for (negative, positive) in [(Some(bad), Some(0.5)), (Some(0.5), Some(bad))] {
            let message = math_message(compute_gene_lo_values(
                &genes(&["A"]),
                &[negative],
                &[positive],
                alpha(0.2, 0.2),
                1,
            ));
            assert!(message.contains("Row 0: non-finite rank fraction"));
        }
    }
}

// T09 / GS-ERR-7: rank zero and negative ranks are invalid in either direction.
#[test]
fn test_rank_validation_zero_and_negative() {
    for bad in [0.0, -0.01] {
        for (negative, positive) in [(Some(bad), Some(0.5)), (Some(0.5), Some(bad))] {
            let message = math_message(compute_gene_lo_values(
                &genes(&["A"]),
                &[negative],
                &[positive],
                alpha(0.2, 0.2),
                1,
            ));
            assert!(message.contains("rank fraction out of bounds"));
        }
    }
}

// T10: ranks above one fail in either direction.
#[test]
fn test_rank_validation_greater_than_one() {
    for (negative, positive) in [(Some(1.05), Some(0.5)), (Some(0.5), Some(1.05))] {
        let message = math_message(compute_gene_lo_values(
            &genes(&["A"]),
            &[negative],
            &[positive],
            alpha(0.2, 0.2),
            1,
        ));
        assert!(message.contains("rank fraction out of bounds"));
    }
}

// T11 / GS-ERR-9: directional alpha ranges and finiteness.
#[test]
fn test_alpha_validation_bounds() {
    for thresholds in [
        alpha(0.04, 0.2),
        alpha(0.51, 0.2),
        alpha(0.2, -0.1),
        alpha(0.2, 1.1),
        alpha(f64::NAN, 0.2),
        alpha(0.2, f64::INFINITY),
    ] {
        let message = math_message(compute_gene_lo_values(
            &genes(&["A"]),
            &[Some(0.5)],
            &[Some(0.5)],
            thresholds,
            1,
        ));
        assert!(message.contains("alpha"));
    }
}

// T12: changing positive ranks cannot change the negative observed score.
#[test]
fn test_directional_independence() {
    let before = target_with_fillers(
        "TARGET",
        &[Some(0.1), Some(0.2)],
        &[Some(0.3), Some(0.4)],
        alpha(0.15, 0.25),
        10,
    );
    let after = target_with_fillers(
        "TARGET",
        &[Some(0.1), Some(0.2)],
        &[Some(0.05), Some(0.08)],
        alpha(0.15, 0.25),
        10,
    );
    near(before.neg_score, 0.19, 1e-14);
    assert_eq!(before.neg_score, after.neg_score);
    assert_ne!(before.pos_score, after.pos_score);
}

// T13 / GS-8: scoring sorts each directional list.
#[test]
fn test_unsorted_ranks_produce_sorted_evaluation() {
    let result = target_with_fillers(
        "UNSORTED",
        &[Some(0.2), Some(0.1)],
        &[Some(0.4), Some(0.2)],
        alpha(0.15, 0.25),
        10,
    );
    near(result.neg_score, 0.19, 1e-14);
    near(result.pos_score, 0.36, 1e-14);
}

// T14: literal score distinguishes global n_valid from gene-local count.
#[test]
fn test_literal_hand_calculated_lo_value_fixture() {
    let result = target_with_fillers(
        "TARGET",
        &[Some(0.05), Some(0.1)],
        &[Some(0.5); 2],
        alpha(0.1, 0.2),
        100,
    );
    near(result.neg_score, 0.01, 1e-14);
    assert!((result.neg_score - 0.994079470779666).abs() > 0.9);
}

// T15: unrelated valid guides do not change the gene-local Beta shape.
#[test]
fn test_universe_size_does_not_change_gene_local_score() {
    let ten = target_with_fillers(
        "TARGET",
        &[Some(0.1), Some(0.2)],
        &[Some(0.5); 2],
        alpha(0.2, 0.2),
        10,
    );
    let two = score(
        &["TARGET", "TARGET"],
        &[Some(0.1), Some(0.2)],
        &[Some(0.5); 2],
        alpha(0.2, 0.2),
        2,
    );
    near(ten.neg_score, 0.04, 1e-14);
    near(two[0].neg_score, 0.04, 1e-14);
    assert_eq!(ten.neg_score, two[0].neg_score);
}

// T16: good-sgRNA comparisons are inclusive.
#[test]
fn test_alpha_boundary_exact_cutoff() {
    let result = score(
        &["TARGET"],
        &[Some(0.25)],
        &[Some(0.25)],
        alpha(0.25, 0.25),
        1,
    );
    assert_eq!((result[0].neg_goodsgrna, result[0].pos_goodsgrna), (1, 1));
}

// T17: directional alpha is bounded by 0.5; zero goods use the smallest rank.
#[test]
fn test_goodsgrna_zero_and_clamped_maximum() {
    let ranks = [Some(0.3), Some(0.4), Some(0.5), Some(0.6)];
    let low = score(&["GENE"; 4], &ranks, &ranks, alpha(0.25, 0.01), 4);
    let high = score(&["GENE"; 4], &ranks, &ranks, alpha(0.5, 0.5), 4);
    assert_eq!((low[0].neg_goodsgrna, low[0].pos_goodsgrna), (0, 0));
    near(low[0].neg_score, 0.7599, 1e-14);
    near(low[0].pos_score, 0.7599, 1e-14);
    assert_eq!((high[0].neg_goodsgrna, high[0].pos_goodsgrna), (3, 3));
}

// T18: keep the strongly depleted literal fixture.
#[test]
fn test_strongly_depleted_synthetic_gene() {
    let result = target_with_fillers(
        "DEPLETED",
        &[Some(0.0001), Some(0.0002), Some(0.0003), Some(0.0004)],
        &[Some(0.5); 4],
        alpha(0.05, 0.2),
        1000,
    );
    near(result.neg_score, 2.56e-14, 1e-25);
    assert!(result.neg_score < 0.001);
}

// T19: keep the strongly enriched literal fixture.
#[test]
fn test_strongly_enriched_synthetic_gene() {
    let result = target_with_fillers(
        "ENRICHED",
        &[Some(0.9996), Some(0.9997), Some(0.9998), Some(0.9999)],
        &[Some(0.0001), Some(0.0002), Some(0.0003), Some(0.0004)],
        alpha(0.05, 0.2),
        1000,
    );
    near(result.pos_score, 2.56e-14, 1e-25);
    assert!(result.neg_score > 0.999);
}

// T20 / GS-9: long interleaved sequence is deterministic.
#[test]
fn test_deterministic_ordering_preserves_input_sequence() {
    let names = [
        "GAPDH", "TP53", "MYC", "BRCA1", "TP53", "EGFR", "MYC", "PTEN", "GAPDH", "CDK4", "EGFR",
        "RB1", "ATM", "BRCA1", "PTEN", "MDM2", "CDK4", "RB1", "ATM", "MDM2",
    ];
    let result = score(
        &names,
        &[Some(0.5); 20],
        &[Some(0.5); 20],
        alpha(0.25, 0.25),
        20,
    );
    assert_eq!(
        result
            .iter()
            .map(|r| r.gene_id.as_str())
            .collect::<Vec<_>>(),
        ["GAPDH", "TP53", "MYC", "BRCA1", "EGFR", "PTEN", "CDK4", "RB1", "ATM", "MDM2"]
    );
}

// T21: CountTable wrapper delegates to the exact in-place API.
#[test]
fn test_from_count_table_convenience_wrapper() {
    let tsv =
        "sgRNA\tGene\tctrl\ttreat\ng1\tGENE_A\t10\t25\ng2\tGENE_B\t5\t50\ng3\tGENE_A\t100\t15\n";
    let table = parse_count_table(Cursor::new(tsv), &["treat"], &["ctrl"]).unwrap();
    let negative = [Some(0.1), None, Some(0.2)];
    let positive = [Some(0.9), None, Some(0.8)];
    let thresholds = alpha(0.2, 0.25);
    let direct =
        compute_gene_lo_values(table.genes(), &negative, &positive, thresholds, 2).unwrap();
    let wrapped =
        compute_gene_lo_values_from_count_table(&table, &negative, &positive, thresholds, 2)
            .unwrap();
    assert_eq!(direct, wrapped);
    assert_eq!(wrapped[1].num_valid_sgrnas, 0);
}

// T22: retained result type has exactly the seven accepted public fields.
#[test]
fn test_result_fields_preserve_baseline_metrics_and_add_valid_count() {
    let result = GeneLoValueResult {
        gene_id: "BRCA1".into(),
        num_sgrnas: 4,
        num_valid_sgrnas: 3,
        neg_score: 0.0012,
        pos_score: 0.9876,
        neg_goodsgrna: 3,
        pos_goodsgrna: 0,
    };
    assert_eq!((result.num_sgrnas, result.num_valid_sgrnas), (4, 3));
    near(result.neg_score, 0.0012, 1e-15);
    near(result.pos_score, 0.9876, 1e-15);
    assert_eq!((result.neg_goodsgrna, result.pos_goodsgrna), (3, 0));
}

// GS-1 and GS-4: physical and valid counts diverge only for masked rows.
#[test]
fn p4f5_gs1_gs4_standard_and_mixed() {
    let standard = target_with_fillers(
        "GENE_A",
        &[Some(0.1), Some(0.2)],
        &[Some(0.2), Some(0.4)],
        alpha(0.15, 0.25),
        10,
    );
    let mixed = target_with_fillers(
        "MIXED",
        &[Some(0.1), None, Some(0.2), None],
        &[Some(0.2), None, Some(0.4), None],
        alpha(0.15, 0.25),
        10,
    );
    for result in [&standard, &mixed] {
        near(result.neg_score, 0.19, 1e-14);
        near(result.pos_score, 0.36, 1e-14);
        assert_eq!((result.neg_goodsgrna, result.pos_goodsgrna), (1, 1));
    }
    assert_eq!((standard.num_sgrnas, standard.num_valid_sgrnas), (2, 2));
    assert_eq!((mixed.num_sgrnas, mixed.num_valid_sgrnas), (4, 2));
}

// GS-3: zero-valid genes remain visible with finite sentinels.
#[test]
fn p4f5_gs3_zero_valid_gene() {
    let result = target_with_fillers("ZERO_VALID", &[None; 3], &[None; 3], alpha(0.25, 0.25), 10);
    assert_eq!((result.num_sgrnas, result.num_valid_sgrnas), (3, 0));
    assert_eq!((result.neg_score, result.pos_score), (1.0, 1.0));
    assert_eq!((result.neg_goodsgrna, result.pos_goodsgrna), (0, 0));
}

// GS-5: one contributing negative guide; no positive goods uses the smallest rank.
#[test]
fn p4f5_gs5_corrected_directional_score() {
    let result = target_with_fillers(
        "ASYMM",
        &[Some(0.04), Some(0.08)],
        &[Some(0.1), Some(0.2)],
        alpha(0.05, 0.01),
        10,
    );
    near(result.neg_score, 0.0784, 1e-14);
    near(result.pos_score, 0.19, 1e-14);
    assert_eq!((result.neg_goodsgrna, result.pos_goodsgrna), (1, 0));
}

// GS-6: ties are separate observations in both scoring and goodsgrna.
#[test]
fn p4f5_gs6_tied_rank_preservation() {
    let result = target_with_fillers(
        "TIED",
        &[Some(0.15), Some(0.15), Some(0.3)],
        &[Some(0.2); 3],
        alpha(0.2, 0.2),
        10,
    );
    near(result.neg_score, 0.06075, 1e-14);
    near(result.pos_score, 0.008, 1e-14);
    assert_eq!((result.neg_goodsgrna, result.pos_goodsgrna), (2, 3));
}

// GS-7: one valid rank at 1.0 remains finite.
#[test]
fn p4f5_gs7_boundary_fraction() {
    let result = score(
        &["BOUNDARY"],
        &[Some(1.0)],
        &[Some(1.0)],
        alpha(0.5, 0.5),
        1,
    );
    assert_eq!((result[0].neg_score, result[0].pos_score), (1.0, 1.0));
    assert_eq!((result[0].neg_goodsgrna, result[0].pos_goodsgrna), (0, 0));
}

// GS-10: duplicate physical guides contribute two ranks, not one.
#[test]
fn p4f5_gs10_duplicate_physical_guides() {
    let result = score(
        &["DUP", "DUP"],
        &[Some(0.1); 2],
        &[Some(0.2); 2],
        alpha(0.1, 0.2),
        2,
    );
    assert_eq!((result[0].num_sgrnas, result[0].num_valid_sgrnas), (2, 2));
    assert_eq!((result[0].neg_goodsgrna, result[0].pos_goodsgrna), (2, 2));
    near(result[0].neg_score, 0.01, 1e-14);
}

#[test]
fn p4f5_opaque_identifiers_and_repeat_determinism() {
    let names = ["A,B", "A", "A;B", "A/B", "A,B", " A "];
    let negative = [Some(0.1), Some(0.2), None, Some(0.3), Some(0.4), Some(0.5)];
    let positive = [Some(0.4), Some(0.3), None, Some(0.2), Some(0.1), Some(0.5)];
    let first = score(&names, &negative, &positive, alpha(0.2, 0.2), 5);
    let second = score(&names, &negative, &positive, alpha(0.2, 0.2), 5);
    assert_eq!(first, second);
    assert_eq!(
        first.iter().map(|r| r.gene_id.as_str()).collect::<Vec<_>>(),
        ["A,B", "A", "A;B", "A/B", " A "]
    );
    assert_eq!(first[0].num_sgrnas, 2);
    assert_eq!(first[2].num_valid_sgrnas, 0);
}

#[test]
fn p4f5_gs_err4_err8_universe_validation() {
    let names = genes(&["A"]);
    let zero = math_message(compute_gene_lo_values(
        &names,
        &[Some(0.2)],
        &[Some(0.2)],
        alpha(0.2, 0.2),
        0,
    ));
    assert_eq!(zero, "Global valid guide universe N_valid cannot be 0");
    let mismatch = math_message(compute_gene_lo_values(
        &names,
        &[Some(0.2)],
        &[Some(0.2)],
        alpha(0.2, 0.2),
        2,
    ));
    assert_eq!(
        mismatch,
        "Observed valid count (1) does not match n_valid (2)"
    );
}

#[test]
fn p4f5_validation_precedence_is_stage_ordered() {
    let names = genes(&[""]);
    let dimension = math_message(compute_gene_lo_values(
        &names,
        &[],
        &[Some(f64::NAN)],
        alpha(f64::NAN, 0.0),
        0,
    ));
    assert!(dimension.starts_with("Length mismatch"));
    let identifier = math_message(compute_gene_lo_values(
        &names,
        &[Some(f64::NAN)],
        &[Some(f64::NAN)],
        alpha(f64::NAN, 0.0),
        0,
    ));
    assert_eq!(identifier, "Row 0: gene identifier is empty");
    let universe = math_message(compute_gene_lo_values(
        &genes(&["A"]),
        &[Some(f64::NAN)],
        &[None],
        alpha(f64::NAN, 0.0),
        0,
    ));
    assert_eq!(universe, "Global valid guide universe N_valid cannot be 0");
    let mask = math_message(compute_gene_lo_values(
        &genes(&["A"]),
        &[Some(f64::NAN)],
        &[None],
        alpha(f64::NAN, 0.0),
        1,
    ));
    assert!(mask.contains("directional mask mismatch"));
    let rank = math_message(compute_gene_lo_values(
        &genes(&["A"]),
        &[Some(f64::NAN)],
        &[Some(0.2)],
        alpha(f64::NAN, 0.0),
        1,
    ));
    assert!(rank.contains("non-finite rank fraction"));
    let count = math_message(compute_gene_lo_values(
        &genes(&["A"]),
        &[Some(0.2)],
        &[Some(0.2)],
        alpha(f64::NAN, 0.0),
        2,
    ));
    assert!(count.contains("Observed valid count"));
}

#[test]
fn p4f5_zero_physical_rows_cannot_name_a_gene() {
    let message = math_message(compute_gene_lo_values(&[], &[], &[], alpha(0.2, 0.2), 1));
    assert_eq!(message, "Input guide slices must not be empty");
}

#[test]
fn p4f5_smallest_positive_rank_is_preserved() {
    let result = score(
        &["A"],
        &[Some(f64::from_bits(1))],
        &[Some(0.5)],
        alpha(0.2, 0.2),
        1,
    );
    assert_eq!(result[0].neg_score, f64::from_bits(1));
    assert!(result[0].neg_score.is_finite());
}

// Exact zero returned by the shared Beta kernel survives both directions.
#[test]
fn test_zero_lo_value_boundary_is_preserved() {
    let ranks = vec![Some(0.01); 1000];
    let result = score(&["CONTROL"; 1000], &ranks, &ranks, alpha(0.25, 0.25), 1000);
    assert_eq!((result[0].neg_score, result[0].pos_score), (0.0, 0.0));
    assert_eq!(
        (result[0].neg_goodsgrna, result[0].pos_goodsgrna),
        (1000, 1000)
    );
    assert_eq!(result[0].num_valid_sgrnas, 1000);
}

#[test]
fn zero_goodsgrna_uses_smallest_valid_negative_rank() {
    let result = target_with_fillers(
        "NO_GOOD",
        &[Some(0.8), None, Some(0.4), Some(0.6)],
        &[Some(0.5), None, Some(0.5), Some(0.5)],
        alpha(0.25, 0.25),
        100,
    );
    assert_eq!((result.num_sgrnas, result.num_valid_sgrnas), (4, 3));
    assert_eq!(result.neg_goodsgrna, 0);
    near(result.neg_score, 0.784, 1e-14);
}
