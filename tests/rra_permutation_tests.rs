//! Numerical P4F6 integration at the crate-private engine boundary.
//! Source modules are included so this test can exercise the accepted
//! crate-private API without making it a public production interface.

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
    compute_empirical_significance, ControlResolutionDiagnostics, GeneLoValueResult,
    NullModelConfig,
};

fn gene(id: &str, k: usize, neg_score: f64, pos_score: f64) -> GeneLoValueResult {
    GeneLoValueResult {
        gene_id: id.into(),
        num_sgrnas: k,
        num_valid_sgrnas: k,
        neg_score,
        pos_score,
        neg_goodsgrna: 0,
        pos_goodsgrna: 0,
    }
}

fn diagnostics() -> ControlResolutionDiagnostics {
    ControlResolutionDiagnostics {
        total_supplied: 2,
        uniquely_resolved: 2,
        unresolved: 0,
        duplicated: 0,
        validity_excluded: 0,
        final_valid: 2,
    }
}

#[test]
fn p4f6_integration_ordinary_seed_and_zero_valid_retention() {
    let mut zero = gene("masked", 0, 1.0, 1.0);
    zero.num_sgrnas = 3;
    let genes = [gene("A", 1, 0.2, 0.8), gene("B", 2, 0.3, 0.7), zero];
    let first = compute_empirical_significance(
        &genes,
        &NullModelConfig::Ordinary,
        200,
        0x5eed,
        4,
        test_alpha(),
    )
    .unwrap();
    let repeated = compute_empirical_significance(
        &genes,
        &NullModelConfig::Ordinary,
        200,
        0x5eed,
        4,
        test_alpha(),
    )
    .unwrap();
    assert_eq!(first, repeated);
    assert_eq!(
        first.iter().map(|g| g.gene_id.as_str()).collect::<Vec<_>>(),
        ["A", "B", "masked"]
    );
    assert_ne!(first[0].neg_pvalue, first[0].pos_pvalue);
    assert_eq!(first[2].num_valid_sgrnas, 0);
    assert_eq!((first[2].neg_pvalue, first[2].pos_pvalue), (1.0, 1.0));
}

#[test]
fn p4f6_integration_control_one_and_directional_pool_isolation() {
    let d = diagnostics();
    let config = NullModelConfig::Control {
        neg_control_percentiles: &[1.0],
        pos_control_percentiles: &[0.2, 0.6],
        diagnostics: &d,
    };
    let genes = [gene("A", 2, 0.5, 0.5)];
    let output = compute_empirical_significance(&genes, &config, 50, 7, 3, test_alpha()).unwrap();
    assert_eq!(output[0].neg_pvalue, 1.0 / 100.0);
    assert!(output[0].pos_pvalue > output[0].neg_pvalue);
    let empty = NullModelConfig::Control {
        neg_control_percentiles: &[],
        pos_control_percentiles: &[0.5],
        diagnostics: &d,
    };
    assert!(matches!(
        compute_empirical_significance(&genes, &empty, 50, 7, 3, test_alpha()),
        Err(MageckError::Config(_))
    ));
}

#[test]
fn p4f6_integration_empirical_epsilon_ties_and_duplicate_scores() {
    let d = diagnostics();
    let config = NullModelConfig::Control {
        neg_control_percentiles: &[0.5],
        pos_control_percentiles: &[0.5],
        diagnostics: &d,
    };
    let genes = [gene("equal", 1, 0.5, 0.5), gene("above", 1, 1.0, 1.0)];
    let output = compute_empirical_significance(&genes, &config, 5, 42, 1, test_alpha()).unwrap();
    assert_eq!(output[0].neg_pvalue, 0.5);
    assert_eq!(output[0].pos_pvalue, 0.5);
    assert_eq!(output[1].neg_pvalue, 0.9);
    assert_eq!(output[1].pos_pvalue, 0.9);
}

fn test_alpha() -> stats::DirectionalAlpha {
    stats::DirectionalAlpha {
        negative: 0.25,
        positive: 0.25,
    }
}
