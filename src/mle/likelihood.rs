#![allow(dead_code, clippy::excessive_precision)]

use crate::error::MageckError;
use std::f64::consts::PI;

const LANCZOS_P: [f64; 9] = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278223245,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
];

/// The smallest-dispersion regime where the bounded Stirling evaluator avoids
/// cancellation between the two large finite-alpha log-gamma terms.
const SMALL_DISPERSION_THRESHOLD: f64 = 1e-4;
const LOG_SQRT_TWO_PI: f64 = 0.918_938_533_204_672_7;
const MAX_EXACT_LOG_FACTORIAL: f64 = 32.0;

// B_(2n) / (2n * (2n - 1)), n = 1..=10.
const STIRLING_BERNOULLI: [f64; 10] = [
    1.0 / 12.0,
    -1.0 / 360.0,
    1.0 / 1260.0,
    -1.0 / 1680.0,
    1.0 / 1188.0,
    -691.0 / 360_360.0,
    1.0 / 156.0,
    -3617.0 / 122_400.0,
    43_867.0 / 244_188.0,
    -174_611.0 / 125_400.0,
];

fn lgamma(mut z: f64) -> f64 {
    if z < 0.5 {
        return PI.ln() - (PI * z).sin().ln() - lgamma(1.0 - z);
    }
    z -= 1.0;
    let mut x = LANCZOS_P[0];
    for (i, coefficient) in LANCZOS_P.iter().enumerate().skip(1) {
        x += coefficient / (z + i as f64);
    }
    let t = z + 7.5;
    (z + 0.5) * t.ln() - t + (2.0 * PI).sqrt().ln() + x.ln()
}

fn check_finite(val: f64, name: &str) -> Result<(), MageckError> {
    if !val.is_finite() {
        return Err(MageckError::Math(format!("{} must be finite", name)));
    }
    Ok(())
}

fn checked_log1p(arg: f64, name: &str) -> Result<f64, MageckError> {
    if !arg.is_finite() || arg <= -1.0 {
        return Err(MageckError::Math(format!(
            "{} logarithm argument must be finite and greater than -1",
            name
        )));
    }
    let value = arg.ln_1p();
    check_finite(value, name)?;
    Ok(value)
}

fn checked_ln(arg: f64, name: &str) -> Result<f64, MageckError> {
    if !arg.is_finite() || arg <= 0.0 {
        return Err(MageckError::Math(format!(
            "{} logarithm argument must be finite and strictly greater than zero",
            name
        )));
    }
    let value = arg.ln();
    check_finite(value, name)?;
    Ok(value)
}

fn stirling_bernoulli_correction(x: f64) -> Result<f64, MageckError> {
    if !x.is_finite() || x <= 0.0 {
        return Err(MageckError::Math(
            "Stirling correction argument must be finite and strictly greater than zero"
                .to_string(),
        ));
    }

    let inverse_square = 1.0 / (x * x);
    check_finite(inverse_square, "Stirling inverse square")?;
    let mut inverse_odd = 1.0 / x;
    check_finite(inverse_odd, "Stirling inverse")?;
    let mut sum = 0.0;
    for coefficient in STIRLING_BERNOULLI {
        let term = coefficient * inverse_odd;
        check_finite(term, "Stirling correction term")?;
        sum += term;
        check_finite(sum, "Stirling correction")?;
        inverse_odd *= inverse_square;
        check_finite(inverse_odd, "Stirling inverse power")?;
    }
    Ok(sum)
}

fn bounded_log_factorial(k: f64) -> Result<f64, MageckError> {
    if k <= MAX_EXACT_LOG_FACTORIAL {
        let mut sum = 0.0;
        for value in 2..=k as u32 {
            sum += checked_ln(f64::from(value), "exact factorial term")?;
            check_finite(sum, "log factorial")?;
        }
        return Ok(sum);
    }

    let correction = stirling_bernoulli_correction(k)?;
    let value = (k + 0.5) * checked_ln(k, "factorial argument")? - k + LOG_SQRT_TWO_PI + correction;
    check_finite(value, "log factorial")?;
    Ok(value)
}

fn small_dispersion_logpmf(k: f64, mu: f64, alpha: f64) -> Result<f64, MageckError> {
    let r = 1.0 / alpha;
    let q = alpha * k;
    let t = alpha * mu;
    check_finite(r, "inverse dispersion")?;
    check_finite(q, "scaled count")?;
    check_finite(t, "scaled mean")?;

    if k == 0.0 {
        let value = -r * checked_log1p(t, "zero-count scaled mean")?;
        check_finite(value, "likelihood term")?;
        return Ok(value);
    }

    let correction_r = stirling_bernoulli_correction(r)?;
    let correction_rk = stirling_bernoulli_correction(r + k)?;
    let log_factorial = bounded_log_factorial(k)?;
    let stirling_factorial =
        (k + 0.5) * checked_ln(k, "factorial Stirling argument")? - k + LOG_SQRT_TWO_PI;
    check_finite(stirling_factorial, "factorial Stirling base")?;
    let r_k_correction = log_factorial - stirling_factorial;
    check_finite(r_k_correction, "factorial correction")?;

    let centered = (mu - k).abs() <= 0.1 * (mu + k);
    let value = if centered {
        let centered_gamma = (q - t) / (1.0 + t);
        let centered_mean = (mu - k) / k;
        check_finite(centered_gamma, "centered gamma logarithm argument")?;
        check_finite(centered_mean, "centered mean logarithm argument")?;
        let first = r * (1.0 + q) * checked_log1p(centered_gamma, "centered gamma")?;
        let second = k * checked_log1p(centered_mean, "centered mean")?;
        let normalizer = checked_log1p(q, "centered scaled count")?
            + checked_ln(2.0 * PI * k, "centered normalizer argument")?;
        check_finite(first, "centered first term")?;
        check_finite(second, "centered second term")?;
        check_finite(normalizer, "centered normalizer")?;
        first + second - 0.5 * normalizer + correction_rk - correction_r - r_k_correction
    } else {
        let constant = (r + k - 0.5) * checked_log1p(q, "noncentered scaled count")? - k
            + correction_rk
            - correction_r;
        let value =
            constant - log_factorial - (r + k) * checked_log1p(t, "noncentered scaled mean")?
                + k * checked_ln(mu, "noncentered mean")?;
        check_finite(constant, "noncentered constant")?;
        value
    };
    check_finite(value, "likelihood term")?;
    Ok(value)
}

fn logpmf(k: f64, mu: f64, alpha: f64) -> Result<f64, MageckError> {
    if alpha <= SMALL_DISPERSION_THRESHOLD {
        return small_dispersion_logpmf(k, mu, alpha);
    }

    let amu = alpha * mu;
    let r = 1.0 / alpha;
    check_finite(amu, "scaled mean")?;
    check_finite(r, "inverse dispersion")?;
    let ln_p = -checked_log1p(amu, "scaled mean")?;
    let ln_1mp = amu.ln() - checked_log1p(amu, "scaled mean")?;
    check_finite(ln_1mp, "success log probability")?;

    let gamma_kr = lgamma(k + r);
    let gamma_r = lgamma(r);
    let gamma_k = lgamma(k + 1.0);
    check_finite(gamma_kr, "log gamma count plus inverse dispersion")?;
    check_finite(gamma_r, "log gamma inverse dispersion")?;
    check_finite(gamma_k, "log gamma count")?;
    let res = gamma_kr - gamma_r - gamma_k + r * ln_p + k * ln_1mp;
    check_finite(res, "likelihood term")?;
    Ok(res)
}

pub(crate) fn nb_log_likelihood_terms(
    counts: &[f64],
    means: &[f64],
    dispersion: &[f64],
) -> Result<Vec<f64>, MageckError> {
    if counts.is_empty() {
        return Err(MageckError::Math("counts must not be empty".to_string()));
    }
    if counts.len() != means.len() {
        return Err(MageckError::Math(
            "counts and means must have equal length".to_string(),
        ));
    }
    if dispersion.len() != 1 && dispersion.len() != counts.len() {
        return Err(MageckError::Math(
            "dispersion length must be 1 or equal to counts length".to_string(),
        ));
    }

    let mut res = Vec::with_capacity(counts.len());
    for i in 0..counts.len() {
        let k = counts[i];
        let mu = means[i];
        let alpha = if dispersion.len() == 1 {
            dispersion[0]
        } else {
            dispersion[i]
        };

        if k < 0.0 || k.fract() != 0.0 || !k.is_finite() {
            return Err(MageckError::Math(
                "count must be a finite non-negative integer".to_string(),
            ));
        }
        if mu <= 0.0 || !mu.is_finite() {
            return Err(MageckError::Math(
                "mean must be finite and strictly greater than zero".to_string(),
            ));
        }
        if alpha <= 0.0 || !alpha.is_finite() {
            return Err(MageckError::Math(
                "dispersion must be finite and strictly greater than zero".to_string(),
            ));
        }

        res.push(logpmf(k, mu, alpha)?);
    }
    Ok(res)
}

pub(crate) fn nb_log_likelihood_sum(
    counts: &[f64],
    means: &[f64],
    dispersion: &[f64],
) -> Result<f64, MageckError> {
    let terms = nb_log_likelihood_terms(counts, means, dispersion)?;
    let mut sum = 0.0;
    for t in terms {
        sum += t;
    }
    check_finite(sum, "likelihood sum")?;
    Ok(sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        let abs_err = (actual - expected).abs();
        let tol = 2e-12 + 2e-12 * expected.abs();
        assert!(
            abs_err <= tol,
            "actual {} != expected {}, err {}",
            actual,
            expected,
            abs_err
        );
    }

    fn assert_slice_close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_close(*a, *e);
        }
    }

    fn assert_close_large_count(actual: f64, expected: f64) {
        const ABS_LOG_TOL_LARGE_COUNT: f64 = 1e-9;
        assert!(
            (actual - expected).abs() <= ABS_LOG_TOL_LARGE_COUNT,
            "actual {} != expected {}, err {}",
            actual,
            expected,
            (actual - expected).abs()
        );
    }

    fn assert_close_near_poisson(actual: f64, expected: f64) {
        const ABS_LOG_TOL_NEAR_POISSON: f64 = 2e-10;
        const RELATIVE_PROBABILITY_TOL_NEAR_POISSON: f64 = 5e-11;
        let log_error = (actual - expected).abs();
        assert!(
            log_error <= ABS_LOG_TOL_NEAR_POISSON,
            "actual {} != expected {}, err {}",
            actual,
            expected,
            log_error
        );
        if expected >= -700.0 {
            assert!(
                (actual - expected).exp_m1().abs() <= RELATIVE_PROBABILITY_TOL_NEAR_POISSON,
                "relative probability error too large for actual {} and expected {}",
                actual,
                expected
            );
        }
    }

    #[test]
    fn test_one_observation() {
        let counts = [3.0];
        let means = [2.5];
        let alpha = [0.2];
        let terms = nb_log_likelihood_terms(&counts, &means, &alpha).unwrap();
        let sum = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        assert_slice_close(&terms, &[-1.7678143450557364]);
        assert_close(sum, -1.7678143450557364);
    }

    #[test]
    fn test_multiple_observations() {
        let counts = [0.0, 2.0, 5.0, 10.0];
        let means = [0.5, 1.75, 4.25, 9.5];
        let alpha = [0.1];
        let terms = nb_log_likelihood_terms(&counts, &means, &alpha).unwrap();
        let sum = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        assert_slice_close(
            &terms,
            &[
                -0.4879016416943201,
                -1.4138231960382419,
                -1.9890057962226386,
                -2.4358762613461078,
            ],
        );
        assert_close(sum, -6.3266068953013086);
    }

    #[test]
    fn test_zero_counts() {
        let counts = [0.0];
        let means = [2.0];
        let alpha = [0.5];
        let terms = nb_log_likelihood_terms(&counts, &means, &alpha).unwrap();
        let sum = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        assert_slice_close(&terms, &[-1.3862943611198906]);
        assert_close(sum, -1.3862943611198906);
    }

    #[test]
    fn test_mixed_zero_and_positive() {
        let counts = [0.0, 4.0, 0.0, 12.0];
        let means = [0.25, 3.5, 8.0, 10.0];
        let alpha = [0.25];
        let terms = nb_log_likelihood_terms(&counts, &means, &alpha).unwrap();
        let sum = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        assert_slice_close(
            &terms,
            &[
                -0.24249848726573942,
                -2.0076467843876693,
                -4.3944491546724391,
                -2.9284212944850792,
            ],
        );
        assert_close(sum, -9.5730157208109272);
    }

    #[test]
    fn test_small_positive_means() {
        let counts = [0.0, 1.0, 2.0];
        let means = [1e-12, 1e-6, 0.01];
        let alpha = [0.2];
        let terms = nb_log_likelihood_terms(&counts, &means, &alpha).unwrap();
        let sum = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        let expected = [
            -0.00000000000099999999999990000000000001333333333333,
            -13.81551175796413507755523499927379874972877182540833,
            -9.73515201438086533825507590066789981380216506354374,
        ];
        assert_slice_close(&terms, &expected);
        assert_close(sum, expected.iter().sum());
    }

    #[test]
    fn test_large_counts() {
        let counts = [1000.0, 10000.0, 100000.0];
        let means = [975.0, 10200.0, 99500.0];
        let alpha = [0.05];
        let terms = nb_log_likelihood_terms(&counts, &means, &alpha).unwrap();
        let sum = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        let expected = [
            -6.34923242490806939821753773394546502027396715754014,
            -8.64046611503661115714836819324808207278032750180608,
            -10.93851579685615489337226573851550070880612149153796,
        ];
        for (actual, expected) in terms.iter().zip(expected) {
            assert_close_large_count(*actual, expected);
        }
        assert_close_large_count(sum, expected.iter().sum());
    }

    #[test]
    fn test_near_zero_dispersion() {
        let counts = [0.0, 3.0, 10.0];
        let means = [0.75, 3.2, 9.5];
        let alpha = [1e-12];
        let terms = nb_log_likelihood_terms(&counts, &means, &alpha).unwrap();
        let sum = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        let expected = [
            -0.74999999999971875000000014062499999992089843750005,
            -1.50230703981249241160796796322781304461405980932488,
            -2.09149458701543878930773326662096749358116780718275,
        ];
        for (actual, expected) in terms.iter().zip(expected) {
            assert_close_near_poisson(*actual, expected);
        }
        assert_close_near_poisson(sum, expected.iter().sum());
    }

    #[test]
    fn test_moderate_dispersion() {
        let counts = [1.0, 5.0, 20.0];
        let means = [1.2, 4.75, 18.0];
        let alpha = [0.5];
        let terms = nb_log_likelihood_terms(&counts, &means, &alpha).unwrap();
        let sum = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        assert_slice_close(
            &terms,
            &[
                -1.2276893309432522,
                -2.3980206136103739,
                -3.6678580614211875,
            ],
        );
        assert_close(sum, -7.293568005974814);
    }

    #[test]
    fn test_large_dispersion() {
        let counts = [0.0, 2.0, 30.0];
        let means = [1.0, 2.5, 25.0];
        let alpha = [100.0];
        let terms = nb_log_likelihood_terms(&counts, &means, &alpha).unwrap();
        let sum = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        assert_slice_close(
            &terms,
            &[
                -0.046151205168412597,
                -5.3516056076252614,
                -8.0570732437643144,
            ],
        );
        assert_close(sum, -13.454830056557988);
    }

    #[test]
    fn test_per_observation_dispersion() {
        let counts = [1.0, 2.0, 3.0];
        let means = [1.5, 2.5, 3.5];
        let alpha = [0.1, 0.2, 0.3];
        let terms = nb_log_likelihood_terms(&counts, &means, &alpha).unwrap();
        let sum = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        assert_slice_close(
            &terms,
            &[
                -1.1319162580185815,
                -1.5164999167748312,
                -1.8474213599807303,
            ],
        );
        assert_close(sum, -4.4958375347741431);
    }

    #[test]
    fn test_deterministic_ordered_sum_and_reproducibility() {
        let counts = [1000.0, 10000.0, 100000.0];
        let means = [975.0, 10200.0, 99500.0];
        let alpha = [0.05];
        let sum1 = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        let sum2 = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        assert_eq!(sum1, sum2); // Exact reproducibility
        assert_close_large_count(sum1, -25.928214336800835448738631665);
    }

    #[test]
    fn test_small_dispersion_routing_and_boundaries_are_finite_deterministic() {
        let below = f64::from_bits(SMALL_DISPERSION_THRESHOLD.to_bits() - 1);
        let at = SMALL_DISPERSION_THRESHOLD;
        let above = f64::from_bits(SMALL_DISPERSION_THRESHOLD.to_bits() + 1);
        for alpha in [below, at, above] {
            let first = logpmf(33.0, 33.0, alpha).unwrap();
            let second = logpmf(33.0, 33.0, alpha).unwrap();
            assert!(first.is_finite());
            assert_eq!(first, second);
        }

        for k in [0.0, 32.0, 33.0] {
            assert!(logpmf(k, k.max(1.0), at).unwrap().is_finite());
        }
    }

    #[test]
    fn test_centered_routing_boundary_and_tail_regimes_are_finite() {
        let k: f64 = 100.0;
        let centered_boundary = 11.0 * k / 9.0;
        let below = f64::from_bits(centered_boundary.to_bits() - 1);
        let above = f64::from_bits(centered_boundary.to_bits() + 1);
        for mu in [below, centered_boundary, above, k, 1e-6, 1e5] {
            let value = small_dispersion_logpmf(k, mu, 1e-8).unwrap();
            assert!(value.is_finite());
        }
    }

    #[test]
    fn test_small_dispersion_terms_preserve_ordered_sum() {
        let counts = [0.0, 32.0, 33.0, 100_000.0];
        let means = [0.75, 32.0, 36.0, 10.0];
        let alpha = [1e-13];
        let terms = nb_log_likelihood_terms(&counts, &means, &alpha).unwrap();
        let sum = nb_log_likelihood_sum(&counts, &means, &alpha).unwrap();
        let mut expected = 0.0;
        for term in terms {
            expected += term;
        }
        assert_eq!(sum, expected);
    }

    #[test]
    fn test_errors() {
        // Empty input
        assert!(nb_log_likelihood_sum(&[], &[], &[0.1]).is_err());

        // Mismatched lengths
        assert!(nb_log_likelihood_sum(&[1.0], &[1.0, 2.0], &[0.1]).is_err());

        // Invalid dispersion length
        assert!(nb_log_likelihood_sum(&[1.0, 2.0], &[1.0, 2.0], &[0.1, 0.2, 0.3]).is_err());

        // Fractional count
        assert!(nb_log_likelihood_sum(&[1.5], &[1.0], &[0.1]).is_err());

        // Negative count
        assert!(nb_log_likelihood_sum(&[-1.0], &[1.0], &[0.1]).is_err());

        // Non-finite count
        assert!(nb_log_likelihood_sum(&[f64::NAN], &[1.0], &[0.1]).is_err());
        assert!(nb_log_likelihood_sum(&[f64::INFINITY], &[1.0], &[0.1]).is_err());

        // Zero mean
        assert!(nb_log_likelihood_sum(&[1.0], &[0.0], &[0.1]).is_err());

        // Negative mean
        assert!(nb_log_likelihood_sum(&[1.0], &[-1.0], &[0.1]).is_err());

        // Non-finite mean
        assert!(nb_log_likelihood_sum(&[1.0], &[f64::NAN], &[0.1]).is_err());

        // Zero dispersion
        assert!(nb_log_likelihood_sum(&[1.0], &[1.0], &[0.0]).is_err());

        // Negative dispersion
        assert!(nb_log_likelihood_sum(&[1.0], &[1.0], &[-0.1]).is_err());

        // Non-finite dispersion
        assert!(nb_log_likelihood_sum(&[1.0], &[1.0], &[f64::INFINITY]).is_err());
    }
}
