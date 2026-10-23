//! P4F7 crate-private integration tests, following the P4F6 source-module pattern.

#[path = "../src/execution.rs"]
mod execution;

mod error {
    pub use mageckrs::error::*;
}
mod stats {
    pub use mageckrs::stats::*;
}
mod mageck_test {
    pub mod count_table {
        pub use mageckrs::mageck_test::count_table::*;
    }
}
#[allow(dead_code)]
#[path = "../src/prng.rs"]
mod prng;
#[allow(dead_code)]
#[path = "../src/mageck_test/rra.rs"]
mod rra;

use error::MageckError;
use rra::{
    apply_directional_bh, compute_empirical_significance, ControlResolutionDiagnostics,
    GeneEmpiricalResult, GeneLoValueResult, NullModelConfig,
};

fn gene(id: &str, physical: usize, valid: usize, neg: f64, pos: f64) -> GeneEmpiricalResult {
    GeneEmpiricalResult {
        gene_id: id.into(),
        num_sgrnas: physical,
        num_valid_sgrnas: valid,
        neg_score: if valid == 0 { 1.0 } else { 0.2 },
        pos_score: if valid == 0 { 1.0 } else { 0.8 },
        neg_pvalue: neg,
        pos_pvalue: pos,
        neg_goodsgrna: if valid == 0 { 0 } else { 1 },
        pos_goodsgrna: 0,
    }
}

fn approx(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
}

fn math_message(result: Result<Vec<rra::GeneFdrResult>, MageckError>) -> String {
    match result {
        Err(MageckError::Math(message)) => message,
        other => panic!("expected Math error, got {other:?}"),
    }
}

#[test]
fn p4f7_empty_and_all_zero_valid_families() {
    assert!(apply_directional_bh(&[]).unwrap().is_empty());
    let input = [
        gene("zero physical", 0, 0, 1.0, 1.0),
        gene("masked", 2, 0, 1.0, 1.0),
    ];
    let output = apply_directional_bh(&input).unwrap();
    assert_eq!(output.len(), 2);
    for (actual, source) in output.iter().zip(&input) {
        assert_eq!(&actual.empirical, source);
        assert_eq!((actual.neg_fdr, actual.pos_fdr), (1.0, 1.0));
    }
}

#[test]
fn p4f7_single_eligible_and_zero_valid_exclusion() {
    let input = [gene("zero", 4, 0, 1.0, 1.0), gene("only", 2, 1, 0.02, 0.7)];
    let output = apply_directional_bh(&input).unwrap();
    assert_eq!((output[0].neg_fdr, output[0].pos_fdr), (1.0, 1.0));
    assert_eq!((output[1].neg_fdr, output[1].pos_fdr), (0.02, 0.7));
    assert_eq!(output[1].empirical, input[1]);
}

#[test]
fn p4f7_directional_families_scatter_and_preserve_all_fields() {
    let input = [
        gene("B", 3, 2, 0.04, 0.01),
        gene("zero", 5, 0, 1.0, 1.0),
        gene("A", 2, 1, 0.01, 0.04),
        gene("C", 1, 1, 0.03, 0.02),
    ];
    let copy = input.clone();
    let output = apply_directional_bh(&input).unwrap();
    assert_eq!(input, copy);
    assert_eq!(
        output
            .iter()
            .map(|g| g.empirical.gene_id.as_str())
            .collect::<Vec<_>>(),
        ["B", "zero", "A", "C"]
    );
    for (actual, source) in output.iter().zip(&input) {
        assert_eq!(&actual.empirical, source);
    }
    for (actual, expected) in
        output
            .iter()
            .zip([(0.04, 0.03), (1.0, 1.0), (0.03, 0.04), (0.04, 0.03)])
    {
        approx(actual.neg_fdr, expected.0);
        approx(actual.pos_fdr, expected.1);
    }
    assert_eq!(apply_directional_bh(&input).unwrap(), output);
}

#[test]
fn p4f7_z_a_arithmetic_ties_duplicates_and_monotonicity() {
    // Z-A valid-input family: [0.01, 0.05, 0.1] -> [0.03, 0.075, 0.1].
    let z = [
        gene("z3", 1, 1, 0.1, 0.1),
        gene("z1", 1, 1, 0.01, 0.01),
        gene("z2", 1, 1, 0.05, 0.05),
    ];
    let out = apply_directional_bh(&z).unwrap();
    for (actual, expected) in out.iter().zip([0.1, 0.03, 0.075]) {
        approx(actual.neg_fdr, expected);
        approx(actual.pos_fdr, expected);
    }
    let tied = [
        gene("a", 1, 1, 0.01, 0.01),
        gene("b", 1, 1, 0.01, 0.01),
        gene("c", 1, 1, 0.05, 0.05),
    ];
    let out = apply_directional_bh(&tied).unwrap();
    for (actual, expected) in out.iter().zip([0.015, 0.015, 0.05]) {
        approx(actual.neg_fdr, expected);
    }
    assert!(out.windows(2).all(|w| w[0].neg_fdr <= w[1].neg_fdr));
}

#[test]
fn p4f7_exact_zero_and_all_one_boundaries() {
    let input = [gene("one", 1, 1, 1.0, 0.0), gene("zero", 1, 1, 0.0, 1.0)];
    let out = apply_directional_bh(&input).unwrap();
    assert_eq!((out[0].neg_fdr, out[0].pos_fdr), (1.0, 0.0));
    assert_eq!((out[1].neg_fdr, out[1].pos_fdr), (0.0, 1.0));
    let all_one = [gene("a", 1, 1, 1.0, 1.0), gene("b", 1, 1, 1.0, 1.0)];
    assert!(apply_directional_bh(&all_one)
        .unwrap()
        .iter()
        .all(|g| (g.neg_fdr, g.pos_fdr) == (1.0, 1.0)));
}

#[test]
fn p4f7_rejects_all_malformed_empirical_values_including_zero_valid() {
    for (direction, value) in [
        ("negative", f64::NAN),
        ("negative", f64::INFINITY),
        ("negative", f64::NEG_INFINITY),
        ("negative", -0.01),
        ("negative", 1.01),
        ("positive", f64::NAN),
        ("positive", f64::INFINITY),
        ("positive", f64::NEG_INFINITY),
        ("positive", -0.01),
        ("positive", 1.01),
    ] {
        let mut input = [
            gene("eligible", 1, 1, 0.1, 0.2),
            gene("masked", 2, 0, 1.0, 1.0),
        ];
        if direction == "negative" {
            input[1].neg_pvalue = value;
        } else {
            input[1].pos_pvalue = value;
        }
        let message = math_message(apply_directional_bh(&input));
        assert!(
            message.contains("empirical p-value validation")
                && message.contains("masked")
                && message.contains("source index 1")
                && message.contains(direction)
                && message.contains(&value.to_string()),
            "{message}"
        );
    }
}

#[test]
fn p4f7_validation_precedence_and_identity() {
    let mut input = [
        gene("a", 1, 1, f64::NAN, f64::NAN),
        gene("b", 1, 1, 0.2, 0.2),
    ];
    input[1].num_valid_sgrnas = 2;
    assert!(math_message(apply_directional_bh(&input)).contains("guide-count validation"));
    input[1].num_valid_sgrnas = 1;
    assert!(math_message(apply_directional_bh(&input)).contains("negative"));
    input[0].neg_pvalue = 0.1;
    assert!(math_message(apply_directional_bh(&input)).contains("positive"));
    input[0].gene_id.clear();
    assert!(math_message(apply_directional_bh(&input)).contains("gene identity validation"));
    input[0].gene_id = "b".into();
    assert!(math_message(apply_directional_bh(&input)).contains("duplicated"));
    input[0].gene_id = "a|opaque".into();
    input[0].pos_pvalue = 0.1;
    assert!(apply_directional_bh(&input).is_ok());
}

fn observed(id: &str, neg_score: f64, pos_score: f64) -> GeneLoValueResult {
    GeneLoValueResult {
        gene_id: id.into(),
        num_sgrnas: 1,
        num_valid_sgrnas: 1,
        neg_score,
        pos_score,
        neg_goodsgrna: 0,
        pos_goodsgrna: 0,
    }
}

#[test]
fn p4f7_ordinary_and_control_p4f6_results_share_the_same_api() {
    let observed = [observed("a", 0.2, 0.8), observed("b", 0.4, 0.6)];
    let ordinary = compute_empirical_significance(
        &observed,
        &NullModelConfig::Ordinary,
        20,
        7,
        2,
        test_alpha(),
    )
    .unwrap();
    let diagnostics = ControlResolutionDiagnostics {
        total_supplied: 2,
        uniquely_resolved: 2,
        unresolved: 0,
        duplicated: 0,
        validity_excluded: 0,
        final_valid: 2,
    };
    let control = NullModelConfig::Control {
        neg_control_percentiles: &[0.2, 0.6],
        pos_control_percentiles: &[0.3, 0.7],
        diagnostics: &diagnostics,
    };
    let controlled =
        compute_empirical_significance(&observed, &control, 20, 7, 2, test_alpha()).unwrap();
    for input in [&ordinary, &controlled] {
        let output = apply_directional_bh(input).unwrap();
        assert_eq!(output.len(), input.len());
        assert!(output.iter().zip(input).all(|(a, b)| &a.empirical == b));
        assert!(output
            .iter()
            .all(|a| (0.0..=1.0).contains(&a.neg_fdr) && (0.0..=1.0).contains(&a.pos_fdr)));
    }
}

fn test_alpha() -> stats::DirectionalAlpha {
    stats::DirectionalAlpha {
        negative: 0.25,
        positive: 0.25,
    }
}
