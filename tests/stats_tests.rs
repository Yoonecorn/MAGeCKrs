use mageckrs::stats::*;
use statrs::distribution::{Beta, ContinuousCDF};
use std::f64;

#[test]
fn test_empirical_p_value_below_all() {
    assert_eq!(empirical_p_value(0.0, &[0.1, 0.2, 0.3], 3), Ok(1.0 / 6.0));
}

#[test]
fn test_empirical_p_value_above_all() {
    assert_eq!(empirical_p_value(1.0, &[0.1, 0.2, 0.3], 3), Ok(5.0 / 6.0));
}

#[test]
fn test_p4f5_e3_kernel_sorted_pool() {
    let pool = [0.1, 0.2, 0.3, 0.4, 0.5];
    assert_eq!(empirical_p_value(0.25, &pool, 5), Ok(2.0 / 5.0));
}

#[test]
fn test_p4f5_e4_exact_pool_member() {
    let pool = [0.1, 0.2, 0.3, 0.4, 0.5];
    assert_eq!(empirical_p_value(0.3, &pool, 5), Ok(1.0 / 2.0));
}

#[test]
fn test_p4f5_e5_midpoint_tie() {
    let pool = [5.0, 15.0];
    assert_eq!(empirical_p_value(10.0, &pool, 2), Ok(1.0 / 2.0));
}

#[test]
fn test_empirical_p_value_exact_nearest_distance_tie() {
    let pool = [0.0, 2.0, 4.0, 6.0];
    assert_eq!(empirical_p_value(3.0, &pool, pool.len()), Ok(0.5));
}

#[test]
fn test_empirical_p_value_lower_epsilon_near_boundary() {
    let pool = [0.0, 1.0, 2.0, 3.0];
    assert_eq!(empirical_p_value(1e-9, &pool, pool.len()), Ok(1.0 / 8.0));
}

#[test]
fn test_empirical_p_value_upper_epsilon_near_boundary() {
    let pool = [0.0, 1.0, 2.0, 3.0];
    assert_eq!(
        empirical_p_value(3.0 - 1e-9, &pool, pool.len()),
        Ok(7.0 / 8.0)
    );
}

#[test]
fn test_empirical_p_value_preserves_duplicate_observations() {
    let pool = [0.1, 0.2, 0.2, 0.2, 0.3];
    assert_eq!(empirical_p_value(0.2, &pool, pool.len()), Ok(0.5));
}

#[test]
fn test_empirical_p_value_singleton_pool() {
    for observed in [-100.0, 0.5, 100.0] {
        assert_eq!(empirical_p_value(observed, &[0.5], 1), Ok(0.5));
    }
}

#[test]
fn test_p4f5_e11_within_epsilon_of_pool_value() {
    let pool = [1.0, 2.0000000001, 3.0];
    assert_eq!(empirical_p_value(2.0, &pool, 3), Ok(1.0 / 2.0));
}

#[test]
fn test_p4f5_e12_rejects_legacy_empty_pool_pseudovalue() {
    assert_eq!(
        empirical_p_value(1.0, &[], 10),
        Err(EmpiricalPValueError::EmptyPool)
    );
}

#[test]
fn test_p4f5_e13_rejects_zero_p_before_search() {
    assert_eq!(
        empirical_p_value(1.0, &[1.0, 2.0], 0),
        Err(EmpiricalPValueError::ZeroPermutationCount)
    );
}

#[test]
fn test_p4f5_e14_rejects_pool_length_mismatch() {
    assert_eq!(
        empirical_p_value(2.0, &[1.0, 3.0], 5),
        Err(EmpiricalPValueError::PoolLengthMismatch {
            expected: 5,
            actual: 2,
        })
    );
}

#[test]
fn test_p4f5_validation_precedence() {
    assert_eq!(
        empirical_p_value(f64::NAN, &[], 0),
        Err(EmpiricalPValueError::NonFiniteObserved)
    );
    assert_eq!(
        empirical_p_value(1.0, &[], 0),
        Err(EmpiricalPValueError::ZeroPermutationCount)
    );
    assert_eq!(
        empirical_p_value(0.2, &[], 1),
        Err(EmpiricalPValueError::EmptyPool)
    );
    assert_eq!(
        empirical_p_value(0.2, &[0.1], 2),
        Err(EmpiricalPValueError::PoolLengthMismatch {
            expected: 2,
            actual: 1,
        })
    );
    assert_eq!(
        empirical_p_value(0.2, &[0.1, f64::NAN], 2),
        Err(EmpiricalPValueError::NonFinitePoolValue { index: 1 })
    );
    assert_eq!(
        empirical_p_value(0.2, &[0.2, 0.1], 2),
        Err(EmpiricalPValueError::PoolNotAscending { index: 1 })
    );
}

#[test]
fn test_empirical_p_value_successes_are_finite_and_strictly_in_range() {
    let pool = [0.1, 0.2, 0.3, 0.4];
    for observed in [-1.0, 0.15, 0.25, 0.35, 1.0] {
        let result = empirical_p_value(observed, &pool, pool.len()).unwrap();
        assert!(result.is_finite());
        assert!(0.0 < result && result < 1.0);
    }
}

#[test]
fn test_empirical_p_value_rejects_zero_permutation_count() {
    assert_eq!(
        empirical_p_value(0.2, &[0.1], 0),
        Err(EmpiricalPValueError::ZeroPermutationCount)
    );
}

#[test]
fn test_empirical_p_value_rejects_empty_pool() {
    assert_eq!(
        empirical_p_value(0.2, &[], 1),
        Err(EmpiricalPValueError::EmptyPool)
    );
}

#[test]
fn test_empirical_p_value_rejects_pool_length_mismatch() {
    assert_eq!(
        empirical_p_value(0.2, &[0.1, 0.2], 3),
        Err(EmpiricalPValueError::PoolLengthMismatch {
            expected: 3,
            actual: 2,
        })
    );
    assert_eq!(
        empirical_p_value(0.2, &[0.1, 0.2, 0.3], 2),
        Err(EmpiricalPValueError::PoolLengthMismatch {
            expected: 2,
            actual: 3,
        })
    );
}

#[test]
fn test_empirical_p_value_rejects_non_finite_observed_values() {
    for observed in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            empirical_p_value(observed, &[0.1], 1),
            Err(EmpiricalPValueError::NonFiniteObserved)
        );
    }
}

#[test]
fn test_empirical_p_value_rejects_first_non_finite_pool_value() {
    assert_eq!(
        empirical_p_value(0.2, &[0.1, f64::NAN, f64::INFINITY], 3),
        Err(EmpiricalPValueError::NonFinitePoolValue { index: 1 })
    );
    assert_eq!(
        empirical_p_value(0.2, &[f64::NEG_INFINITY], 1),
        Err(EmpiricalPValueError::NonFinitePoolValue { index: 0 })
    );
}

#[test]
fn test_empirical_p_value_rejects_first_descending_pair() {
    assert_eq!(
        empirical_p_value(0.2, &[0.1, 0.3, 0.2, 0.4], 4),
        Err(EmpiricalPValueError::PoolNotAscending { index: 2 })
    );
}

#[test]
fn test_beta_cdf_edges() {
    assert_eq!(beta_cdf(0.0, 1.0, 1.0), 0.0);
    assert_eq!(beta_cdf(1.0, 1.0, 1.0), 1.0);
    assert_eq!(beta_cdf(-0.1, 1.0, 1.0), 0.0);
    assert_eq!(beta_cdf(1.1, 1.0, 1.0), 1.0);
}

#[test]
fn test_beta_cdf_values() {
    // Beta(1, 1) is Uniform(0, 1) -> CDF(x) = x
    let cdf = beta_cdf(0.3, 1.0, 1.0);
    assert!((cdf - 0.3).abs() < 1e-14, "Expected 0.3, got {}", cdf);

    // Beta(2, 2) CDF at 0.5 is exactly 0.5 due to symmetry
    let cdf2 = beta_cdf(0.5, 2.0, 2.0);
    assert!((cdf2 - 0.5).abs() < 1e-14, "Expected 0.5, got {}", cdf2);

    // Asymmetric Beta(2, 5) CDF at 0.3 (direct statrs reference check)
    let b = Beta::new(2.0, 5.0).unwrap();
    let expected = b.cdf(0.3);
    let actual = beta_cdf(0.3, 2.0, 5.0);
    assert!((actual - expected).abs() < 1e-14);
}

#[test]
fn test_beta_cdf_invalid() {
    assert!(beta_cdf(0.5, -1.0, 1.0).is_nan());
    assert!(beta_cdf(0.5, 1.0, -1.0).is_nan());
    assert!(beta_cdf(0.5, 0.0, 1.0).is_nan());
    assert!(beta_cdf(f64::NAN, 1.0, 1.0).is_nan());
    assert!(beta_cdf(0.5, f64::NAN, 1.0).is_nan());
    assert!(beta_cdf(0.5, f64::INFINITY, 1.0).is_nan());
    assert!(beta_cdf(f64::INFINITY, 1.0, 1.0).is_nan());
    assert!(beta_cdf(f64::NEG_INFINITY, 1.0, 1.0).is_nan());
}

#[test]
fn test_lo_value_empty_single() {
    assert_eq!(lo_value(&[], 10), 1.0);
    // Single rank 0.5 with n_total = 1. j=1, Beta(1, 1). CDF(0.5) = 0.5
    let lv = lo_value(&[0.5], 1);
    assert!((lv - 0.5).abs() < 1e-14);
}

#[test]
fn test_lo_value_hand_calculated() {
    // Two contributing ranks: min(1 - 0.9^2, 0.2^2) = 0.04.
    let lv = lo_value(&[0.1, 0.2], 2);
    assert!((lv - 0.04).abs() < 1e-14);
}

#[test]
fn test_lo_value_n_includes_unsupplied_valid_gene_ranks() {
    let ranks = [0.1];
    assert!((lo_value(&ranks, ranks.len()) - 0.1).abs() < 1e-14);
    assert!((lo_value(&ranks, 10) - 0.6513215599).abs() < 1e-14);
}

#[test]
fn test_lo_value_rank_fraction_1() {
    let lv = lo_value(&[1.0], 1);
    assert_eq!(lv, 1.0);
}

#[test]
fn test_lo_value_borrowed_slice() {
    let ranks = [0.1, 0.5];
    let _ = lo_value(&ranks, ranks.len());
    assert_eq!(ranks, [0.1, 0.5]); // input unchanged
}

#[test]
fn test_bh_correction_empty_single() {
    assert_eq!(bh_correction(&[]), Vec::<f64>::new());
    assert_eq!(bh_correction(&[0.05]), vec![0.05]);
}

#[test]
fn test_bh_correction_already_sorted() {
    // p = [0.01, 0.05]
    // adj = [0.01 * 2 / 1 = 0.02, 0.05 * 2 / 2 = 0.05]
    let adj = bh_correction(&[0.01, 0.05]);
    assert!((adj[0] - 0.02).abs() < 1e-14);
    assert!((adj[1] - 0.05).abs() < 1e-14);
}

#[test]
fn test_bh_correction_reverse_sorted() {
    let adj = bh_correction(&[0.05, 0.01]);
    assert!((adj[0] - 0.05).abs() < 1e-14);
    assert!((adj[1] - 0.02).abs() < 1e-14);
}

#[test]
fn test_bh_correction_ties() {
    let adj = bh_correction(&[0.01, 0.01, 0.05]);
    // Sorted: 0.01, 0.01, 0.05. Ranks: 1, 2, 3
    // For p=0.05 (rank 3): adj = 0.05 * 3 / 3 = 0.05. min_adj = 0.05
    // For second 0.01 (rank 2): adj = 0.01 * 3 / 2 = 0.015. min_adj = 0.015
    // For first 0.01 (rank 1): adj = 0.01 * 3 / 1 = 0.03. But min_adj = 0.015!
    // So both 0.01 get 0.015.
    assert!((adj[0] - 0.015).abs() < 1e-14);
    assert!((adj[1] - 0.015).abs() < 1e-14);
    assert!((adj[2] - 0.05).abs() < 1e-14);
}

#[test]
fn test_bh_correction_monotonicity() {
    let pvals = [0.01, 0.02, 0.03, 0.04, 0.05];
    let adj = bh_correction(&pvals);
    for i in 0..adj.len() - 1 {
        assert!(adj[i] <= adj[i + 1]);
    }
}

#[test]
fn test_bh_correction_bounds() {
    let adj = bh_correction(&[0.0, 1.0, 0.9]);
    assert_eq!(adj[0], 0.0);
    assert_eq!(adj[1], 1.0);
    assert_eq!(adj[2], 1.0); // 0.9 * 3 / 2 = 1.35 -> capped to 1.0
}

#[test]
fn test_bh_correction_borrowed_input() {
    let p = [0.1, 0.2];
    let _ = bh_correction(&p);
    assert_eq!(p, [0.1, 0.2]);
}

#[test]
fn test_bh_correction_invalid_inputs() {
    let p = [0.01, f64::NAN, -0.1, 1.5, 0.05];
    let adj = bh_correction(&p);

    // n_valid = 2 (0.01, 0.05).
    // Sorted valid: 0.01 (rank 1), 0.05 (rank 2).
    // 0.05 -> 0.05 * 2 / 2 = 0.05. min_adj = 0.05.
    // 0.01 -> 0.01 * 2 / 1 = 0.02. min_adj = 0.02.
    assert!((adj[0] - 0.02).abs() < 1e-14);
    assert!(adj[1].is_nan());
    assert!(adj[2].is_nan());
    assert!(adj[3].is_nan());
    assert!((adj[4] - 0.05).abs() < 1e-14);
}

#[test]
fn test_assign_ranks_empty_single() {
    assert_eq!(assign_ranks(&[]), Vec::<f64>::new());
    assert_eq!(assign_ranks(&[0.5]), vec![1.0]);
}

#[test]
fn test_assign_ranks_distinct() {
    let r = assign_ranks(&[0.5, 0.1, 0.8]);
    // Sorted: 0.1, 0.5, 0.8
    // Ranks: 1, 2, 3 -> fractions: 1/3, 2/3, 1
    assert!((r[1] - 1.0 / 3.0).abs() < 1e-14);
    assert!((r[0] - 2.0 / 3.0).abs() < 1e-14);
    assert!((r[2] - 1.0).abs() < 1e-14);
}

#[test]
fn test_assign_ranks_ties() {
    let r = assign_ranks(&[0.5, 0.5, 0.1]);
    // Sorted pvals: 0.1, 0.5, 0.5
    // Ranks:
    // 0.1 -> 1/3
    // 0.5, 0.5 -> average rank = (2 + 3) / 2 = 2.5
    // fraction = 2.5 / 3 = 0.8333333333333334
    assert!((r[2] - 1.0 / 3.0).abs() < 1e-14);
    assert!((r[0] - 2.5 / 3.0).abs() < 1e-14);
    assert!((r[1] - 2.5 / 3.0).abs() < 1e-14);
}

#[test]
fn test_assign_ranks_bounds() {
    let r = assign_ranks(&[0.0, 1.0, 0.5]);
    assert!((r[0] - 1.0 / 3.0).abs() < 1e-14);
    assert!((r[2] - 2.0 / 3.0).abs() < 1e-14);
    assert_eq!(r[1], 1.0);
}

#[test]
fn test_assign_ranks_zero_based_failure() {
    // The requirement says: a test that would fail under zero-based rank fractions
    let r = assign_ranks(&[0.1]);
    assert_eq!(r[0], 1.0); // If 0-based, it would be 0.0 / 1 = 0.0
    assert!(r[0] > 0.0);
}

#[test]
fn test_assign_ranks_average_success() {
    // Spec explicitly requires mid-rank assignment in (0, 1]
    let r = assign_ranks(&[0.1, 0.1]);
    // With mid-rank, both receive (1 + 2) / 2 = 1.5 -> fraction 1.5 / 2 = 0.75
    assert!((r[0] - 0.75).abs() < 1e-14);
    assert!((r[1] - 0.75).abs() < 1e-14);
}

#[test]
fn test_assign_ranks_invalid_inputs() {
    let r = assign_ranks(&[0.1, f64::NAN, 0.5]);
    assert!((r[0] - 1.0 / 2.0).abs() < 1e-14);
    assert!(r[1].is_nan());
    assert!((r[2] - 2.0 / 2.0).abs() < 1e-14);
}

#[test]
fn test_permutation_seed_basic() {
    assert_eq!(permutation_seed(0, 0), 0);
    assert_eq!(permutation_seed(0, 1), 1);
    assert_eq!(permutation_seed(1, 0), 0x9e3779b97f4a7c15);
    let out1 = permutation_seed(42, 5);
    let out2 = permutation_seed(42, 5);
    assert_eq!(out1, out2);
}

#[test]
fn test_permutation_seed_wrapping() {
    // a test that would fail if checked or saturating arithmetic were used
    let base = u64::MAX;
    let seed = permutation_seed(base, u64::MAX);
    let expected = base.wrapping_mul(0x9e3779b97f4a7c15).wrapping_add(u64::MAX);
    assert_eq!(seed, expected);
}

#[test]
fn test_assign_ranks_infinities_and_zeros() {
    let pvals = [0.0, f64::INFINITY, -0.0, f64::NEG_INFINITY];
    let r = assign_ranks(&pvals);
    // valid_indexed ordered: NEG_INFINITY (idx 3), -0.0 (idx 2), 0.0 (idx 0), INFINITY (idx 1)
    // Ties:
    // idx 3: rank 1 -> fraction 1/4 = 0.25
    // idx 2 & idx 0 are tied (0.0 == -0.0). Base ranks 2 and 3. Average rank = 2.5 -> fraction 2.5/4 = 0.625
    // idx 1: rank 4 -> fraction 4/4 = 1.0
    assert!((r[3] - 0.25).abs() < 1e-14);
    assert!((r[2] - 0.625).abs() < 1e-14);
    assert!((r[0] - 0.625).abs() < 1e-14);
    assert!((r[1] - 1.0).abs() < 1e-14);
}

#[test]
fn test_assign_ranks_negative_zero_tie_explicit() {
    // Explicit test for [-0.0, 0.0, 1.0]
    // -0.0 and 0.0 occupy one-based positions 1 and 2 -> mid-rank 1.5 / 3 = 0.5
    // 1.0 occupies position 3 -> 3 / 3 = 1.0
    let r1 = assign_ranks(&[-0.0, 0.0, 1.0]);
    assert_eq!(r1.len(), 3);
    assert!((r1[0] - 0.5).abs() < 1e-14);
    assert!((r1[1] - 0.5).abs() < 1e-14);
    assert!((r1[2] - 1.0).abs() < 1e-14);

    // Permuted input order [0.0, -0.0, 1.0]
    let r2 = assign_ranks(&[0.0, -0.0, 1.0]);
    assert_eq!(r2.len(), 3);
    assert!((r2[0] - 0.5).abs() < 1e-14);
    assert!((r2[1] - 0.5).abs() < 1e-14);
    assert!((r2[2] - 1.0).abs() < 1e-14);
}

#[test]
fn test_assign_ranks_repeated_infinities() {
    // Repeated positive infinity
    let r_pos = assign_ranks(&[f64::INFINITY, f64::INFINITY]);
    assert_eq!(r_pos.len(), 2);
    assert!((r_pos[0] - 0.75).abs() < 1e-14);
    assert!((r_pos[1] - 0.75).abs() < 1e-14);

    // Repeated negative infinity
    let r_neg = assign_ranks(&[f64::NEG_INFINITY, f64::NEG_INFINITY]);
    assert_eq!(r_neg.len(), 2);
    assert!((r_neg[0] - 0.75).abs() < 1e-14);
    assert!((r_neg[1] - 0.75).abs() < 1e-14);
}

#[test]
fn test_lo_value_structural_inputs() {
    // 1. n_total = 0 with empty rank slice -> 1.0
    assert_eq!(lo_value(&[], 0), 1.0);

    // 2. n_total = 0 with nonempty rank slice -> NaN (ranks.len() > n_total)
    assert!(lo_value(&[0.5], 0).is_nan());

    // 3. ranks.len() > n_total -> NaN (e.g. ranks.len() = 3, n_total = 2)
    assert!(lo_value(&[0.1, 0.2, 0.3], 2).is_nan());

    // 4. ranks.len() = n_total -> valid finite f64
    let lv_eq = lo_value(&[0.1, 0.2], 2);
    assert!(!lv_eq.is_nan());
    // j=1: Beta(1, 2).cdf(0.1) = 1 - (1 - 0.1)^2 = 0.19
    // j=2: Beta(2, 1).cdf(0.2) = 0.2^2 = 0.04
    // min is 0.04
    assert!((lv_eq - 0.04).abs() < 1e-14);

    // 5. valid ranks mixed with NaN -> NaN
    assert!(lo_value(&[0.1, f64::NAN], 2).is_nan());

    // 6. valid ranks mixed with out-of-range ranks -> NaN
    assert!(lo_value(&[-0.1, 0.5], 2).is_nan());
    assert!(lo_value(&[0.5, 1.2], 2).is_nan());
}

#[test]
fn test_directional_alpha_positive_passing_all_but_no_valid_pass_is_clamped() {
    let result = calc_directional_alpha(&[0.5, 0.5], &[0.01, 0.5], &[false, true], 0.05).unwrap();
    assert_eq!(result.positive, 0.05);
}

#[test]
fn test_directional_alpha_positive_fallback_uses_passing_all() {
    let result = calc_directional_alpha(&[0.5, 0.5], &[0.1, 0.5], &[false, true], 0.05).unwrap();
    assert_eq!(result.positive, 0.01);
}

#[test]
fn test_directional_alpha_positive_valid_fraction() {
    let result = calc_directional_alpha(&[0.5, 0.5], &[0.01, 0.01], &[true, true], 0.05).unwrap();
    assert_eq!(result.positive, 0.5);
}

#[test]
fn test_clamp_alpha_below_lower_boundary() {
    assert_eq!(clamp_alpha(0.02), 0.05);
    let result = calc_directional_alpha(&[0.5], &[0.5], &[true], 0.05).unwrap();
    assert_eq!(result.negative, 0.05);
}

#[test]
fn test_clamp_alpha_exact_lower_boundary() {
    assert_eq!(clamp_alpha(0.05), 0.05);
}

#[test]
fn test_clamp_alpha_interior() {
    assert_eq!(clamp_alpha(0.25), 0.25);
}

#[test]
fn test_clamp_alpha_exact_upper_boundary() {
    assert_eq!(clamp_alpha(0.5), 0.5);
}

#[test]
fn test_clamp_alpha_above_upper_boundary() {
    assert_eq!(clamp_alpha(0.6), 0.5);
}

#[test]
fn test_directional_alpha_threshold_is_inclusive() {
    let result = calc_directional_alpha(&[0.05], &[0.5], &[true], 0.05).unwrap();
    assert_eq!(result.negative, 0.5);
}

#[test]
fn test_directional_alpha_positive_uses_directional_clamp() {
    let result = calc_directional_alpha(&[0.5; 3], &[0.01, 0.02, 0.5], &[true; 3], 0.05).unwrap();
    assert_eq!(result.positive, 0.5);
}

#[test]
fn test_directional_alpha_rejects_empty_valid_universe() {
    assert_eq!(
        calc_directional_alpha(&[], &[], &[], 0.05),
        Err(DirectionalAlphaError::EmptyValidUniverse)
    );
}

#[test]
fn test_directional_alpha_rejects_length_mismatch_with_actual_lengths() {
    assert_eq!(
        calc_directional_alpha(&[0.1], &[0.1, 0.2], &[true], 0.05),
        Err(DirectionalAlphaError::LengthMismatch {
            negative_len: 1,
            positive_len: 2,
            valid_mask_len: 1,
        })
    );
}

#[test]
fn test_directional_alpha_rejects_all_non_finite_thresholds() {
    for threshold in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            calc_directional_alpha(&[0.1], &[0.1], &[true], threshold),
            Err(DirectionalAlphaError::NonFiniteThreshold)
        );
    }
}

#[test]
fn test_directional_alpha_rejects_out_of_range_thresholds() {
    for threshold in [-0.1, 1.1] {
        assert_eq!(
            calc_directional_alpha(&[0.1], &[0.1], &[true], threshold),
            Err(DirectionalAlphaError::ThresholdOutOfRange)
        );
    }
}

#[test]
fn test_directional_alpha_rejects_first_invalid_negative_value_including_masked() {
    assert_eq!(
        calc_directional_alpha(&[0.1, f64::NAN, f64::INFINITY], &[0.1; 3], &[true; 3], 0.05),
        Err(DirectionalAlphaError::NonFiniteNegativePValue { index: 1 })
    );
    assert_eq!(
        calc_directional_alpha(&[-0.1, 1.1], &[0.1; 2], &[false, true], 0.05),
        Err(DirectionalAlphaError::NegativePValueOutOfRange { index: 0 })
    );
}

#[test]
fn test_directional_alpha_rejects_first_invalid_positive_value_including_masked() {
    assert_eq!(
        calc_directional_alpha(
            &[0.1; 3],
            &[0.1, f64::NEG_INFINITY, f64::NAN],
            &[true; 3],
            0.05
        ),
        Err(DirectionalAlphaError::NonFinitePositivePValue { index: 1 })
    );
    assert_eq!(
        calc_directional_alpha(&[0.1; 2], &[1.1, -0.1], &[false, true], 0.05),
        Err(DirectionalAlphaError::PositivePValueOutOfRange { index: 0 })
    );
}

#[test]
fn test_directional_alpha_validation_precedence() {
    assert_eq!(
        calc_directional_alpha(&[f64::NAN], &[], &[true], f64::NAN),
        Err(DirectionalAlphaError::LengthMismatch {
            negative_len: 1,
            positive_len: 0,
            valid_mask_len: 1,
        })
    );
    assert_eq!(
        calc_directional_alpha(&[f64::NAN], &[f64::NAN], &[true], f64::NAN),
        Err(DirectionalAlphaError::NonFiniteThreshold)
    );
    assert_eq!(
        calc_directional_alpha(&[f64::NAN], &[f64::NAN], &[true], 1.1),
        Err(DirectionalAlphaError::ThresholdOutOfRange)
    );
    assert_eq!(
        calc_directional_alpha(&[f64::NAN], &[f64::NAN], &[true], 0.05),
        Err(DirectionalAlphaError::NonFiniteNegativePValue { index: 0 })
    );
}
