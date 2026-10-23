use statrs::distribution::{Beta, ContinuousCDF};

use crate::error::MageckError;

#[allow(dead_code)] // CNV pipeline wiring is explicitly deferred beyond CNV2.
pub(crate) mod brent;
#[allow(dead_code)] // CNV pipeline wiring is deferred beyond CNV3.
pub(crate) mod cnv;
#[allow(unused_imports)]
pub(crate) use brent::{brent_minimize, OptimizationDiagnostics, OptimizationError, StopReason};

#[allow(dead_code)] // The defensive check is reached only through the deferred CNV consumer.
fn checked_type7_interpolation_result(result: f64) -> Result<f64, MageckError> {
    if !result.is_finite() {
        return Err(MageckError::Math(
            "Type 7 percentile interpolation produced a non-finite result".to_string(),
        ));
    }
    Ok(result)
}

/// Computes the Type 7 sample percentile without mutating the caller's slice.
///
/// Validation is deliberately ordered: quantile finiteness, quantile range,
/// input emptiness, then input finiteness in original index order.
#[allow(dead_code)] // CNV pipeline wiring is explicitly deferred to CNV2.
pub(crate) fn type7_percentile(values: &[f64], q: f64) -> Result<f64, MageckError> {
    if !q.is_finite() {
        return Err(MageckError::Math(format!(
            "Percentile quantile q must be finite, got {q}"
        )));
    }
    if !(0.0..=1.0).contains(&q) {
        return Err(MageckError::Math(format!(
            "Percentile quantile q must be in [0.0, 1.0], got {q}"
        )));
    }
    if values.is_empty() {
        return Err(MageckError::Math(
            "Percentile input slice cannot be empty".to_string(),
        ));
    }
    for (index, &value) in values.iter().enumerate() {
        if !value.is_finite() {
            return Err(MageckError::Math(format!(
                "Percentile input contains non-finite value at index {index}: {value}"
            )));
        }
    }

    let mut sorted = values.to_vec();
    sorted.sort_unstable_by(f64::total_cmp);

    let n = sorted.len();
    let h = (n - 1) as f64 * q;
    let i = h.floor() as usize;
    let gamma = h - i as f64;

    if i >= n - 1 {
        return Ok(sorted[n - 1]);
    }

    let a = sorted[i];
    let b = sorted[i + 1];
    if gamma == 0.0 {
        return Ok(a);
    }

    let result = if a < 0.0 && b > 0.0 {
        (1.0 - gamma) * a + gamma * b
    } else {
        a + gamma * (b - a)
    };
    checked_type7_interpolation_result(result)
}

/// Computes the fixed 99.9th percentile used for a future CNV upper bound.
#[allow(dead_code)] // CNV pipeline wiring is explicitly deferred to CNV2.
pub(crate) fn type7_percentile_999(values: &[f64]) -> Result<f64, MageckError> {
    type7_percentile(values, 0.999)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DirectionalAlpha {
    pub negative: f64,
    pub positive: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DirectionalAlphaError {
    LengthMismatch {
        negative_len: usize,
        positive_len: usize,
        valid_mask_len: usize,
    },
    NonFiniteThreshold,
    ThresholdOutOfRange,
    NonFiniteNegativePValue {
        index: usize,
    },
    NegativePValueOutOfRange {
        index: usize,
    },
    NonFinitePositivePValue {
        index: usize,
    },
    PositivePValueOutOfRange {
        index: usize,
    },
    EmptyValidUniverse,
    NonFiniteResult,
    ResultOutOfRange,
}

/// Clamps a finite directional alpha fraction to the accepted interval.
pub fn clamp_alpha(raw_alpha: f64) -> f64 {
    raw_alpha.clamp(0.05, 0.5)
}

/// Calculates independent guide-level alpha values for negative and positive selection.
pub fn calc_directional_alpha(
    negative_p: &[f64],
    positive_p: &[f64],
    valid_mask: &[bool],
    threshold: f64,
) -> Result<DirectionalAlpha, DirectionalAlphaError> {
    if negative_p.len() != positive_p.len() || negative_p.len() != valid_mask.len() {
        return Err(DirectionalAlphaError::LengthMismatch {
            negative_len: negative_p.len(),
            positive_len: positive_p.len(),
            valid_mask_len: valid_mask.len(),
        });
    }
    if !threshold.is_finite() {
        return Err(DirectionalAlphaError::NonFiniteThreshold);
    }
    if !(0.0..=1.0).contains(&threshold) {
        return Err(DirectionalAlphaError::ThresholdOutOfRange);
    }

    for (index, &value) in negative_p.iter().enumerate() {
        if !value.is_finite() {
            return Err(DirectionalAlphaError::NonFiniteNegativePValue { index });
        }
        if !(0.0..=1.0).contains(&value) {
            return Err(DirectionalAlphaError::NegativePValueOutOfRange { index });
        }
    }
    for (index, &value) in positive_p.iter().enumerate() {
        if !value.is_finite() {
            return Err(DirectionalAlphaError::NonFinitePositivePValue { index });
        }
        if !(0.0..=1.0).contains(&value) {
            return Err(DirectionalAlphaError::PositivePValueOutOfRange { index });
        }
    }

    let n_valid = valid_mask.iter().filter(|&&is_valid| is_valid).count();
    if n_valid == 0 {
        return Err(DirectionalAlphaError::EmptyValidUniverse);
    }

    let negative_passing_valid = negative_p
        .iter()
        .zip(valid_mask)
        .filter(|&(&p_value, &is_valid)| is_valid && p_value <= threshold)
        .count();
    let positive_passing_all = positive_p
        .iter()
        .filter(|&&p_value| p_value <= threshold)
        .count();
    let positive_passing_valid = positive_p
        .iter()
        .zip(valid_mask)
        .filter(|&(&p_value, &is_valid)| is_valid && p_value <= threshold)
        .count();

    let negative = clamp_alpha(negative_passing_valid as f64 / n_valid as f64);
    let positive = if positive_passing_all == 0 {
        0.01
    } else {
        clamp_alpha(positive_passing_valid as f64 / n_valid as f64)
    };

    if !negative.is_finite() || !positive.is_finite() {
        return Err(DirectionalAlphaError::NonFiniteResult);
    }
    if !(0.05..=0.5).contains(&negative) || !(0.01..=0.5).contains(&positive) {
        return Err(DirectionalAlphaError::ResultOutOfRange);
    }

    Ok(DirectionalAlpha { negative, positive })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmpiricalPValueError {
    NonFiniteObserved,
    ZeroPermutationCount,
    EmptyPool,
    PoolLengthMismatch { expected: usize, actual: usize },
    NonFinitePoolValue { index: usize },
    PoolNotAscending { index: usize },
    NonFiniteResult,
    ResultOutOfRange,
}

/// Computes the empirical p-value for an observed lo-value against a sorted null pool.
///
/// The two searches reproduce MAGeCK2's nearest-index binary search, including its
/// preference for the upper index when a value is equidistant from adjacent entries.
pub fn empirical_p_value(
    observed_lo: f64,
    null_pool: &[f64],
    permutation_count: usize,
) -> Result<f64, EmpiricalPValueError> {
    if !observed_lo.is_finite() {
        return Err(EmpiricalPValueError::NonFiniteObserved);
    }
    if permutation_count == 0 {
        return Err(EmpiricalPValueError::ZeroPermutationCount);
    }
    if null_pool.is_empty() {
        return Err(EmpiricalPValueError::EmptyPool);
    }
    if null_pool.len() != permutation_count {
        return Err(EmpiricalPValueError::PoolLengthMismatch {
            expected: permutation_count,
            actual: null_pool.len(),
        });
    }

    for (index, &value) in null_pool.iter().enumerate() {
        if !value.is_finite() {
            return Err(EmpiricalPValueError::NonFinitePoolValue { index });
        }
        if index > 0 && null_pool[index - 1] > value {
            return Err(EmpiricalPValueError::PoolNotAscending { index });
        }
    }

    fn nearest_index(value: f64, values: &[f64]) -> usize {
        let mut lo = 0;
        let mut hi = values.len() - 1;

        if value <= values[lo] {
            return lo;
        }
        if value >= values[hi] {
            return hi;
        }

        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if value >= values[mid] {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        if (values[hi] - value).abs() > (values[lo] - value).abs() {
            lo
        } else {
            hi
        }
    }

    let epsilon = 1e-9_f64;
    let left = nearest_index(observed_lo - epsilon, null_pool);
    let right = nearest_index(observed_lo + epsilon, null_pool);
    let result = (left as f64 + right as f64 + 1.0) / (2.0 * permutation_count as f64);

    if !result.is_finite() {
        return Err(EmpiricalPValueError::NonFiniteResult);
    }
    if !(0.0 < result && result < 1.0) {
        return Err(EmpiricalPValueError::ResultOutOfRange);
    }

    Ok(result)
}

/// Computes the exact API and behavior required by the specification for the Beta CDF.
///
/// Handles `x <= 0.0` (returns 0.0) and `x >= 1.0` (returns 1.0) mathematically.
/// Invalid shape parameters return `NaN` without panicking.
pub fn beta_cdf(x: f64, a: f64, b: f64) -> f64 {
    if x.is_nan()
        || !x.is_finite()
        || a.is_nan()
        || b.is_nan()
        || !a.is_finite()
        || !b.is_finite()
        || a <= 0.0
        || b <= 0.0
    {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }

    match Beta::new(a, b) {
        Ok(dist) => dist.cdf(x),
        Err(_) => f64::NAN,
    }
}

/// Computes the MAGeCK-compatible RRA lo-value.
///
/// The non-negotiable invariant is:
/// `lo_value = min_{j=1..k} Beta(j, n-j+1).CDF(r_j)`
/// where `n_total` is the number of sgRNAs targeting this specific gene
/// before alpha filtering, including valid guides whose ranks are not supplied,
/// NOT the global universe size.
///
/// `ranks` must already be sorted in ascending order. This function does not sort internally.
/// A filtered subset is accepted when `ranks.len() <= n_total`.
/// Empty ranks return 1.0.
///
/// Structurally invalid inputs return `f64::NAN`:
/// - Non-empty `ranks` with `ranks.len() > n_total` (including `n_total == 0`)
/// - Any rank fraction outside `[0.0, 1.0]` or equal to `NaN`
pub fn lo_value(ranks: &[f64], n_total: usize) -> f64 {
    if ranks.is_empty() {
        return 1.0;
    }

    if ranks.len() > n_total {
        return f64::NAN;
    }

    let n = n_total as f64;
    let mut min_val = 1.0;

    for (idx, &r_j) in ranks.iter().enumerate() {
        if !(0.0..=1.0).contains(&r_j) {
            return f64::NAN;
        }

        let j = (idx + 1) as f64;
        let a = j;
        let b = n - j + 1.0;

        let cdf = beta_cdf(r_j, a, b);
        if cdf.is_nan() {
            return f64::NAN;
        }
        if cdf < min_val {
            min_val = cdf;
        }
    }

    min_val
}

/// Applies Benjamini-Hochberg multiple-testing correction.
///
/// Ties are resolved deterministically by preserving original index order.
/// P-values are sorted ascending, keeping original index for ties.
/// The reverse cumulative minimum is applied in sorted order to ensure monotonicity.
pub fn bh_correction(pvalues: &[f64]) -> Vec<f64> {
    if pvalues.is_empty() {
        return vec![];
    }
    let n_total = pvalues.len();

    let mut valid_indexed: Vec<(usize, f64)> = pvalues
        .iter()
        .copied()
        .enumerate()
        .filter(|&(_, p)| (0.0..=1.0).contains(&p))
        .collect();

    // Sort by value ascending, then by original index ascending.
    valid_indexed.sort_unstable_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));

    let n_valid = valid_indexed.len();
    let mut adj = vec![f64::NAN; n_total];
    let mut min_adj = f64::INFINITY;

    // Traverse from largest p-value to smallest p-value (which is n_valid-1 down to 0)
    for i in (0..n_valid).rev() {
        let rank = (i + 1) as f64;
        let (orig_idx, pval) = valid_indexed[i];

        let mut cur_adj = pval * (n_valid as f64) / rank;
        if cur_adj > min_adj {
            cur_adj = min_adj;
        } else {
            min_adj = cur_adj;
        }

        if cur_adj > 1.0 {
            cur_adj = 1.0;
        }

        adj[orig_idx] = cur_adj;
    }

    adj
}

/// Assigns rank fractions in `(0, 1]` based on mid-rank.
///
/// Resolves ties using average mid-ranks. For identical values, all elements
/// receive the average of their would-be ordinal rank fractions.
///
/// Tie grouping uses standard IEEE 754 numeric equality (`==`), grouping `-0.0` and `+0.0`
/// as numeric ties into a single mid-rank group. `total_cmp` is used solely for deterministic
/// total sorting order before grouping.
///
/// Invalid inputs (`NaN`) are excluded from `n_valid` and return `NaN`.
/// Positive and negative infinity are treated as valid inputs (sorting last and first, respectively;
/// identical infinities form tie groups).
pub fn assign_ranks(pvalues: &[f64]) -> Vec<f64> {
    if pvalues.is_empty() {
        return vec![];
    }

    let n_total = pvalues.len();
    let mut valid_indexed: Vec<(usize, f64)> = pvalues
        .iter()
        .copied()
        .enumerate()
        .filter(|&(_, v)| !v.is_nan())
        .collect();

    // Sort by value ascending, then by original index ascending.
    valid_indexed.sort_unstable_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));

    let mut ranks = vec![f64::NAN; n_total];
    let n_valid = valid_indexed.len();

    if n_valid > 0 {
        let mut i = 0;
        while i < n_valid {
            let mut j = i;
            while j < n_valid && valid_indexed[j].1 == valid_indexed[i].1 {
                j += 1;
            }

            // Items from i to j-1 are tied.
            // Base ranks are (i+1) to j. Average rank is (i + 1 + j) / 2.0
            let mid_rank = (i + 1 + j) as f64 / 2.0;
            let rank_fraction = mid_rank / (n_valid as f64);

            for k in i..j {
                ranks[valid_indexed[k].0] = rank_fraction;
            }

            i = j;
        }
    }

    ranks
}

/// Assigns mid-rank fractions over only the structurally valid guides.
///
/// All supplied values, including masked values, must be finite. Results retain
/// source order, with `None` for every masked guide.
#[allow(dead_code)] // P4F4 pipeline wiring is deferred to a later task.
pub(crate) fn assign_masked_rank_fractions(
    pvalues: &[f64],
    valid_mask: &[bool],
) -> Result<Vec<Option<f64>>, crate::error::MageckError> {
    if pvalues.len() != valid_mask.len() {
        return Err(MageckError::Math(format!(
            "Masked rank fractions require equal lengths: {} p-values and {} mask entries",
            pvalues.len(),
            valid_mask.len()
        )));
    }

    for (index, &value) in pvalues.iter().enumerate() {
        if !value.is_finite() {
            return Err(MageckError::Math(format!(
                "Masked rank fractions require finite p-values; index {index} has {value}"
            )));
        }
    }

    let n_valid = valid_mask.iter().filter(|&&valid| valid).count();
    if n_valid == 0 {
        return Err(MageckError::Math(
            "Masked rank fractions require at least one valid guide".to_string(),
        ));
    }

    let mut indexed: Vec<(usize, f64)> = pvalues
        .iter()
        .copied()
        .enumerate()
        .filter(|(index, _)| valid_mask[*index])
        .collect();
    indexed.sort_unstable_by(|left, right| {
        left.1
            .total_cmp(&right.1)
            .then_with(|| left.0.cmp(&right.0))
    });

    let mut result = vec![None; pvalues.len()];
    let mut left = 0;
    while left < n_valid {
        let mut right = left + 1;
        while right < n_valid && indexed[right].1 == indexed[left].1 {
            right += 1;
        }

        // L = left, k = right - left. Convert before arithmetic to avoid
        // integer overflow in (2L + k + 1) / (2 * N_valid).
        let mid_rank = left as f64 + ((right - left) as f64 + 1.0) / 2.0;
        let fraction = mid_rank / n_valid as f64;
        if !fraction.is_finite() || !(0.0 < fraction && fraction <= 1.0) {
            return Err(MageckError::Math(format!(
                "Masked rank fraction is out of range for tie group at sorted positions {left}..{right}: {fraction}"
            )));
        }
        for &(source_index, _) in &indexed[left..right] {
            result[source_index] = Some(fraction);
        }
        left = right;
    }

    Ok(result)
}

/// Generates a deterministic permutation seed from a base and task ID.
pub fn permutation_seed(base: u64, task_id: u64) -> u64 {
    base.wrapping_mul(0x9e3779b97f4a7c15).wrapping_add(task_id)
}

#[cfg(test)]
mod tests {
    use super::{
        assign_masked_rank_fractions, checked_type7_interpolation_result, type7_percentile,
        type7_percentile_999,
    };
    use crate::error::MageckError;

    #[test]
    fn p4f15b_f1_filtered_subset_retains_full_gene_n() {
        let score = super::lo_value(&[0.01, 0.02], 3);
        assert!((score - 0.001184).abs() < 1e-14);
    }

    #[test]
    fn p4f15b_f2_full_vector_winner_is_second_rank() {
        let ranks = [0.2, 0.4, 0.8];
        let candidates: Vec<f64> = ranks
            .iter()
            .enumerate()
            .map(|(index, &rank)| super::beta_cdf(rank, (index + 1) as f64, (3 - index) as f64))
            .collect();
        assert!((candidates[0] - 0.488).abs() < 1e-14);
        assert!((candidates[1] - 0.352).abs() < 1e-14);
        assert!((candidates[2] - 0.512).abs() < 1e-14);
        assert_eq!(
            candidates
                .iter()
                .enumerate()
                .min_by(|a, b| a.1.total_cmp(b.1))
                .unwrap()
                .0,
            1
        );
        assert!((super::lo_value(&ranks, 3) - 0.352).abs() < 1e-14);
    }

    #[test]
    fn p4f15b_f3_first_rank_uses_full_gene_n() {
        assert!((super::lo_value(&[0.8], 2) - 0.96).abs() < 1e-14);
    }

    #[test]
    fn p4f15b_subset_contract_preserves_invalid_and_empty_boundaries() {
        assert_eq!(super::lo_value(&[], 0), 1.0);
        assert_eq!(super::lo_value(&[], 3), 1.0);
        assert!(super::lo_value(&[0.1], 0).is_nan());
        assert!(super::lo_value(&[0.1, 0.2], 1).is_nan());
        for rank in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.1, 1.1] {
            assert!(super::lo_value(&[rank], 3).is_nan());
        }
        assert_eq!(super::lo_value(&[0.0], 3), 0.0);
        assert_eq!(super::lo_value(&[1.0], 3), 1.0);
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1e-12,
            "expected {expected}, got {actual}"
        );
    }

    fn math_message(result: Result<f64, MageckError>) -> String {
        match result {
            Err(MageckError::Math(message)) => message,
            Err(other) => panic!("expected MageckError::Math, got {other:?}"),
            Ok(value) => panic!("expected error, got {value}"),
        }
    }

    fn masked_math_message(pvalues: &[f64], mask: &[bool]) -> String {
        match assign_masked_rank_fractions(pvalues, mask) {
            Err(MageckError::Math(message)) => message,
            Err(other) => panic!("expected MageckError::Math, got {other:?}"),
            Ok(value) => panic!("expected error, got {value:?}"),
        }
    }

    #[test]
    fn masked_ranks_reject_both_length_mismatches_before_value_validation() {
        for mask in [&[true][..], &[true, false, true][..]] {
            let message = masked_math_message(&[f64::NAN, 0.5], mask);
            assert!(message.contains("equal lengths"), "{message}");
        }
    }

    #[test]
    fn masked_ranks_reject_all_nonfinite_values_before_exclusion() {
        for (values, mask, index) in [
            (vec![f64::NAN, 1.0], vec![true, true], 0),
            (vec![1.0, f64::INFINITY], vec![true, true], 1),
            (vec![f64::NEG_INFINITY, 1.0], vec![true, true], 0),
            (vec![1.0, f64::NAN], vec![true, false], 1),
            (vec![f64::INFINITY, 1.0], vec![false, true], 0),
            (vec![1.0, f64::NEG_INFINITY], vec![true, false], 1),
        ] {
            let message = masked_math_message(&values, &mask);
            assert!(message.contains(&format!("index {index}")), "{message}");
        }
        let message = masked_math_message(&[1.0, f64::NAN, f64::INFINITY], &[true, false, true]);
        assert!(message.contains("index 1"), "{message}");
        let message = masked_math_message(&[f64::NAN], &[false]);
        assert!(message.contains("index 0"), "{message}");
    }

    #[test]
    fn masked_ranks_reject_empty_valid_universe() {
        for (values, mask) in [(vec![0.1, 0.2], vec![false, false]), (vec![], vec![])] {
            let message = masked_math_message(&values, &mask);
            assert!(message.contains("at least one valid guide"), "{message}");
        }
    }

    #[test]
    fn masked_ranks_single_valid_and_ascending_unique() {
        assert_eq!(
            assign_masked_rank_fractions(&[9.0, 0.3, -4.0], &[false, true, false]).unwrap(),
            vec![None, Some(1.0), None]
        );
        assert_eq!(
            assign_masked_rank_fractions(&[0.1, 0.2, 0.3, 0.4], &[true; 4]).unwrap(),
            vec![Some(0.25), Some(0.5), Some(0.75), Some(1.0)]
        );
    }

    #[test]
    fn masked_ranks_restore_source_order_and_ignore_masked_values() {
        let values = [0.3, -100.0, 0.1, 100.0, 0.2];
        let mask = [true, false, true, false, true];
        let expected = [Some(1.0), None, Some(1.0 / 3.0), None, Some(2.0 / 3.0)];
        for _ in 0..3 {
            assert_eq!(
                assign_masked_rank_fractions(&values, &mask).unwrap(),
                expected
            );
        }
    }

    #[test]
    fn masked_ranks_assign_exact_mid_ranks_to_all_tie_members() {
        for (values, expected) in [
            (vec![1.0, 2.0, 2.0, 4.0], vec![0.25, 0.625, 0.625, 1.0]),
            (vec![1.0, 2.0, 2.0, 2.0, 4.0], vec![0.2, 0.6, 0.6, 0.6, 1.0]),
            (vec![1.0, 1.0, 3.0], vec![0.5, 0.5, 1.0]),
            (vec![1.0, 3.0, 3.0], vec![1.0 / 3.0, 5.0 / 6.0, 5.0 / 6.0]),
            (vec![1.0, 1.0, 1.0], vec![2.0 / 3.0; 3]),
        ] {
            let mask = vec![true; values.len()];
            let actual = assign_masked_rank_fractions(&values, &mask).unwrap();
            assert_eq!(actual, expected.into_iter().map(Some).collect::<Vec<_>>());
            assert_eq!(actual.len(), values.len());
        }
    }

    #[test]
    fn masked_ranks_preserve_duplicates_and_finite_boundaries() {
        assert_eq!(
            assign_masked_rank_fractions(&[f64::MAX, -0.0, 0.0, -f64::MAX, 0.0], &[true; 5])
                .unwrap(),
            vec![Some(1.0), Some(0.6), Some(0.6), Some(0.2), Some(0.6)]
        );
        assert_eq!(
            assign_masked_rank_fractions(&[0.0, f64::from_bits(1)], &[true; 2]).unwrap(),
            vec![Some(0.5), Some(1.0)]
        );
    }

    #[test]
    fn masked_ranks_output_structure_and_domain_follow_mask() {
        let values = [7.0, 2.0, 2.0, -1.0, 9.0];
        let mask = [false, true, true, true, false];
        let result = assign_masked_rank_fractions(&values, &mask).unwrap();
        assert_eq!(result.len(), values.len());
        for (rank, &valid) in result.iter().zip(&mask) {
            assert_eq!(rank.is_some(), valid);
            if let Some(fraction) = rank {
                assert!(fraction.is_finite() && 0.0 < *fraction && *fraction <= 1.0);
            }
        }
        assert_eq!(
            result,
            vec![
                None,
                Some(5.0 / 6.0),
                Some(5.0 / 6.0),
                Some(1.0 / 3.0),
                None
            ]
        );
    }

    #[test]
    fn type7_f1_singleton() {
        assert_eq!(type7_percentile(&[42.0], 0.999).unwrap(), 42.0);
    }

    #[test]
    fn type7_f2_two_elements() {
        assert_close(type7_percentile(&[1.0, 10.0], 0.999).unwrap(), 9.991);
    }

    #[test]
    fn type7_f3_unsorted() {
        assert_close(type7_percentile(&[10.0, 1.0], 0.999).unwrap(), 9.991);
    }

    #[test]
    fn type7_f4_ten_elements() {
        let values: Vec<f64> = (1..=10).map(f64::from).collect();
        assert_close(type7_percentile(&values, 0.999).unwrap(), 9.991);
    }

    #[test]
    fn type7_f5_duplicates() {
        assert_eq!(type7_percentile(&[5.0, 5.0], 0.999).unwrap(), 5.0);
    }

    #[test]
    fn type7_f6_negative() {
        assert_close(type7_percentile(&[-10.0, -1.0], 0.999).unwrap(), -1.009);
    }

    #[test]
    fn type7_f7_exact_bound() {
        let values: Vec<f64> = (0..=1000).map(f64::from).collect();
        assert_eq!(type7_percentile(&values, 0.999).unwrap(), 999.0);
    }

    #[test]
    fn type7_f8_empty() {
        assert_eq!(
            math_message(type7_percentile(&[], 0.999)),
            "Percentile input slice cannot be empty"
        );
    }

    #[test]
    fn type7_f9_nan_index() {
        assert_eq!(
            math_message(type7_percentile(&[1.0, f64::NAN, 3.0], 0.999)),
            "Percentile input contains non-finite value at index 1: NaN"
        );
    }

    #[test]
    fn type7_f10_positive_infinity_index() {
        assert_eq!(
            math_message(type7_percentile(&[f64::INFINITY, 2.0], 0.999)),
            "Percentile input contains non-finite value at index 0: inf"
        );
    }

    #[test]
    fn type7_f11_negative_infinity_index() {
        assert_eq!(
            math_message(type7_percentile(&[1.0, 2.0, f64::NEG_INFINITY], 0.999)),
            "Percentile input contains non-finite value at index 2: -inf"
        );
    }

    #[test]
    fn type7_f12_first_nonfinite_precedence() {
        assert_eq!(
            math_message(type7_percentile(
                &[1.0, f64::NAN, 2.0, f64::INFINITY],
                0.999
            )),
            "Percentile input contains non-finite value at index 1: NaN"
        );
    }

    #[test]
    fn type7_f13_input_is_immutable_on_error() {
        let values = [3.0, f64::NAN, 1.0];
        let before: Vec<u64> = values.iter().map(|value| value.to_bits()).collect();
        let _ = type7_percentile(&values, 0.999);
        assert_eq!(
            values
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            before
        );
    }

    #[test]
    fn type7_f14_degenerate_bound() {
        assert_eq!(type7_percentile(&[1.0, 1.0, 1.0], 0.999).unwrap(), 1.0);
    }

    #[test]
    fn type7_f15_odd_length() {
        assert_close(type7_percentile(&[1.0, 2.0, 3.0], 0.999).unwrap(), 2.998);
    }

    #[test]
    fn type7_f16_minimum_boundary() {
        assert_eq!(type7_percentile(&[1.0, 10.0], 0.0).unwrap(), 1.0);
    }

    #[test]
    fn type7_f17_maximum_boundary() {
        assert_eq!(type7_percentile(&[1.0, 10.0], 1.0).unwrap(), 10.0);
    }

    #[test]
    fn type7_f18_median() {
        assert_eq!(type7_percentile(&[1.0, 10.0], 0.5).unwrap(), 5.5);
    }

    #[test]
    fn type7_f19_nan_quantile() {
        assert_eq!(
            math_message(type7_percentile(&[1.0, 10.0], f64::NAN)),
            "Percentile quantile q must be finite, got NaN"
        );
    }

    #[test]
    fn type7_f20_positive_infinite_quantile() {
        assert_eq!(
            math_message(type7_percentile(&[1.0, 10.0], f64::INFINITY)),
            "Percentile quantile q must be finite, got inf"
        );
    }

    #[test]
    fn type7_f21_negative_quantile() {
        assert_eq!(
            math_message(type7_percentile(&[1.0, 10.0], -0.001)),
            "Percentile quantile q must be in [0.0, 1.0], got -0.001"
        );
    }

    #[test]
    fn type7_f22_greater_than_one_quantile() {
        assert_eq!(
            math_message(type7_percentile(&[1.0, 10.0], 1.001)),
            "Percentile quantile q must be in [0.0, 1.0], got 1.001"
        );
    }

    #[test]
    fn type7_f23_negative_zero_and_positive_zero() {
        let result = type7_percentile(&[-0.0, 0.0], 0.999).unwrap();
        assert_eq!(result.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn type7_f24_input_is_immutable_on_success() {
        let values: [f64; 2] = [10.0, 1.0];
        let before: Vec<u64> = values.iter().map(|value| value.to_bits()).collect();
        let _ = type7_percentile(&values, 0.999).unwrap();
        assert_eq!(
            values
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            before
        );
    }

    #[test]
    fn type7_integer_gamma_preserves_negative_zero_bits() {
        assert_eq!(
            type7_percentile(&[-0.0, 1.0], 0.0).unwrap().to_bits(),
            (-0.0_f64).to_bits()
        );
    }

    #[test]
    fn type7_zero_and_positive_endpoint_use_difference_interpolation() {
        assert_eq!(type7_percentile(&[0.0, 10.0], 0.5).unwrap(), 5.0);
    }

    #[test]
    fn type7_negative_and_zero_endpoint_use_difference_interpolation() {
        assert_eq!(type7_percentile(&[-10.0, 0.0], 0.5).unwrap(), -5.0);
    }

    #[test]
    fn type7_negative_to_positive_uses_weighted_interpolation() {
        assert_eq!(type7_percentile(&[-f64::MAX, f64::MAX], 0.5).unwrap(), 0.0);
    }

    #[test]
    fn type7_extreme_finite_same_sign_endpoints_remain_finite() {
        let result = type7_percentile(&[f64::MAX / 2.0, f64::MAX], 0.5).unwrap();
        assert!(result.is_finite());
        assert!(result > f64::MAX / 2.0 && result < f64::MAX);
    }

    #[test]
    fn type7_nonfinite_interpolation_error_is_exact() {
        assert_eq!(
            math_message(checked_type7_interpolation_result(f64::INFINITY)),
            "Type 7 percentile interpolation produced a non-finite result"
        );
    }

    #[test]
    fn type7_validation_precedence_is_exact() {
        assert_eq!(
            math_message(type7_percentile(&[], f64::NAN)),
            "Percentile quantile q must be finite, got NaN"
        );
        assert_eq!(
            math_message(type7_percentile(&[], -1.0)),
            "Percentile quantile q must be in [0.0, 1.0], got -1"
        );
        assert_eq!(
            math_message(type7_percentile(&[f64::NAN], 0.5)),
            "Percentile input contains non-finite value at index 0: NaN"
        );
    }

    #[test]
    fn type7_fixed_wrapper_delegates_to_literal_999_quantile() {
        let values = [1.0, 10.0];
        assert_eq!(
            type7_percentile_999(&values).unwrap().to_bits(),
            type7_percentile(&values, 0.999).unwrap().to_bits()
        );
    }
}
