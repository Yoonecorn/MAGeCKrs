//! Deterministic bounded scalar minimization using Brent's method.

/// Machine epsilon square root used by the SciPy-parity tolerance rule.
pub(crate) const SQRT_EPS: f64 = 1.4901161193847656e-8;
/// Golden-section contraction factor `(3 - sqrt(5)) / 2`.
pub(crate) const GOLDEN_MEAN: f64 = 0.38196601125010515;
/// Fixed absolute parameter tolerance for the production optimizer.
pub(crate) const DEFAULT_XATOL: f64 = 1e-5;
/// Fixed objective-evaluation ceiling for the production optimizer.
pub(crate) const DEFAULT_MAXFUN: usize = 500;

/// Stopping reason recorded upon successful optimization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum StopReason {
    IntervalToleranceMet,
    DegenerateIntervalBypass,
}

/// Structured diagnostics returned by the bounded scalar optimizer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct OptimizationDiagnostics {
    pub minimizer: f64,
    pub objective_value: f64,
    pub iterations: usize,
    pub func_evals: usize,
    pub final_bracket: (f64, f64),
    pub converged: bool,
    pub stop_reason: StopReason,
}

/// Errors produced at the crate-private optimizer boundary.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum OptimizationError {
    NonFiniteBounds(f64),
    InvalidBounds {
        lower: f64,
        upper: f64,
    },
    NonFiniteObjective {
        x: f64,
        fx: f64,
    },
    MaxEvaluationsExceeded {
        best_x: f64,
        best_fx: f64,
        evaluations: usize,
        final_bracket: (f64, f64),
    },
}

/// Minimizes `func` over the finite closed interval `[lower, upper]`.
pub(crate) fn brent_minimize<F>(
    func: F,
    lower: f64,
    upper: f64,
) -> Result<OptimizationDiagnostics, OptimizationError>
where
    F: FnMut(f64) -> f64,
{
    brent_minimize_kernel(func, lower, upper, DEFAULT_XATOL, DEFAULT_MAXFUN)
}

#[cfg(test)]
pub(crate) fn brent_minimize_with_budget<F>(
    func: F,
    lower: f64,
    upper: f64,
    xatol: f64,
    maxfun: usize,
) -> Result<OptimizationDiagnostics, OptimizationError>
where
    F: FnMut(f64) -> f64,
{
    assert!(xatol.is_finite(), "xatol must be finite, got {xatol}");
    assert!(
        xatol > 0.0,
        "xatol must be strictly positive (> 0.0), got {xatol}"
    );
    assert!(
        maxfun > 0,
        "maxfun must be strictly positive (> 0), got {maxfun}"
    );
    brent_minimize_kernel(func, lower, upper, xatol, maxfun)
}

fn brent_minimize_kernel<F>(
    mut func: F,
    lower: f64,
    upper: f64,
    xatol: f64,
    maxfun: usize,
) -> Result<OptimizationDiagnostics, OptimizationError>
where
    F: FnMut(f64) -> f64,
{
    if !lower.is_finite() {
        return Err(OptimizationError::NonFiniteBounds(lower));
    }
    if !upper.is_finite() {
        return Err(OptimizationError::NonFiniteBounds(upper));
    }
    if lower > upper {
        return Err(OptimizationError::InvalidBounds { lower, upper });
    }
    if lower == upper {
        let fx = func(lower);
        if !fx.is_finite() {
            return Err(OptimizationError::NonFiniteObjective { x: lower, fx });
        }
        return Ok(OptimizationDiagnostics {
            minimizer: lower,
            objective_value: fx,
            iterations: 0,
            func_evals: 1,
            final_bracket: (lower, upper),
            converged: true,
            stop_reason: StopReason::DegenerateIntervalBypass,
        });
    }

    let mut a = lower;
    let mut b = upper;
    let width = b - a;
    let xf_primary = a + GOLDEN_MEAN * width;
    let mut xf = if width.is_finite() && xf_primary.is_finite() {
        xf_primary
    } else {
        (1.0 - GOLDEN_MEAN) * a + GOLDEN_MEAN * b
    };
    if !xf.is_finite() {
        xf = 0.5 * a + 0.5 * b;
    }
    xf = xf.clamp(a, b);
    debug_assert!(xf.is_finite() && xf >= a && xf <= b);

    let mut fx = func(xf);
    let mut num = 1;
    if !fx.is_finite() {
        return Err(OptimizationError::NonFiniteObjective { x: xf, fx });
    }

    let mut fulc = xf;
    let mut nfc = xf;
    let mut ffulc = fx;
    let mut fnfc = fx;
    let mut rat: f64 = 0.0;
    let mut e: f64 = 0.0;
    let mut iterations = 0;

    loop {
        let xm_primary = 0.5 * (a + b);
        let xm = if xm_primary.is_finite() {
            xm_primary
        } else {
            0.5 * a + 0.5 * b
        }
        .clamp(a, b);
        let tol_1 = SQRT_EPS * xf.abs() + xatol / 3.0;
        let tol_2 = 2.0 * tol_1;
        let half_width_primary = 0.5 * (b - a);
        let half_width = if half_width_primary.is_finite() {
            half_width_primary
        } else {
            0.5 * b - 0.5 * a
        };

        if (xf - xm).abs() <= tol_2 - half_width {
            return Ok(OptimizationDiagnostics {
                minimizer: xf,
                objective_value: fx,
                iterations,
                func_evals: num,
                final_bracket: (a, b),
                converged: true,
                stop_reason: StopReason::IntervalToleranceMet,
            });
        }
        if num >= maxfun {
            return Err(OptimizationError::MaxEvaluationsExceeded {
                best_x: xf,
                best_fx: fx,
                evaluations: num,
                final_bracket: (a, b),
            });
        }

        let mut golden = true;
        if e.abs() > tol_1 {
            let mut r = (xf - nfc) * (fx - ffulc);
            let mut q = (xf - fulc) * (fx - fnfc);
            let mut p = (xf - fulc) * q - (xf - nfc) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }
            r = e;
            e = rat;
            let candidate = xf + p / q;
            if p.abs() < (0.5 * q * r).abs()
                && p > q * (a - xf)
                && p < q * (b - xf)
                && candidate.is_finite()
            {
                rat = p / q;
                if (candidate - a) < tol_2 || (b - candidate) < tol_2 {
                    let si = if xm >= xf { 1.0 } else { -1.0 };
                    rat = tol_1 * si;
                }
                golden = false;
            }
        }

        if golden {
            let (bound, e_primary) = if xf >= xm { (a, a - xf) } else { (b, b - xf) };
            e = e_primary;
            let rat_primary = GOLDEN_MEAN * e_primary;
            rat = if e_primary.is_finite() && rat_primary.is_finite() {
                rat_primary
            } else {
                GOLDEN_MEAN * bound - GOLDEN_MEAN * xf
            };
        }

        let si = if rat >= 0.0 { 1.0 } else { -1.0 };
        let step = si * rat.abs().max(tol_1);
        let mut u = xf + step;
        if !u.is_finite() {
            u = xm;
        }
        let u = u.clamp(a, b);
        debug_assert!(u.is_finite() && u >= a && u <= b);
        let fu = func(u);
        num += 1;
        iterations += 1;
        if !fu.is_finite() {
            return Err(OptimizationError::NonFiniteObjective { x: u, fx: fu });
        }

        if fu <= fx {
            if u >= xf {
                a = xf;
            } else {
                b = xf;
            }
            fulc = nfc;
            ffulc = fnfc;
            nfc = xf;
            fnfc = fx;
            xf = u;
            fx = fu;
        } else {
            if u < xf {
                a = u;
            } else {
                b = u;
            }
            if fu <= fnfc || nfc == xf {
                fulc = nfc;
                ffulc = fnfc;
                nfc = u;
                fnfc = fu;
            } else if fu <= ffulc || fulc == xf || fulc == nfc {
                fulc = u;
                ffulc = fu;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::{Cell, RefCell};

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual} != {expected}"
        );
    }

    #[test]
    fn o1_lower_boundary_minimum() {
        let result = brent_minimize(|x| (x - 0.5).powi(2), 1.0, 5.0).unwrap();
        assert_close(result.minimizer, 1.000_005_628_7, 1e-6);
        assert_eq!(result.func_evals, 28);
        assert_eq!(result.stop_reason, StopReason::IntervalToleranceMet);
    }

    #[test]
    fn o2_upper_boundary_minimum() {
        let result = brent_minimize(|x| (x - 6.0).powi(2), 1.0, 5.0).unwrap();
        assert_close(result.minimizer, 4.999_994_371_3, 1e-6);
        assert_eq!(result.func_evals, 28);
    }

    #[test]
    fn o3_interior_quadratic_minimum() {
        let result = brent_minimize(|x| (x - 2.5).powi(2), 1.0, 5.0).unwrap();
        assert_close(result.minimizer, 2.5, 1e-7);
        assert_eq!(result.func_evals, 6);
    }

    #[test]
    fn o4_flat_objective_favors_later_ties() {
        let result = brent_minimize(|_| 10.0, 1.0, 5.0).unwrap();
        assert_close(result.minimizer, 4.999_994_371_3, 1e-6);
        assert_eq!(result.func_evals, 28);
    }

    #[test]
    fn o5_concave_objective_reaches_an_endpoint() {
        let result = brent_minimize(|x| -(x - 3.0).powi(2), 1.0, 5.0).unwrap();
        assert!((result.minimizer - 1.0).abs() < 1e-5 || (result.minimizer - 5.0).abs() < 1e-5);
        assert_eq!(result.func_evals, 28);
    }

    #[test]
    fn o6_rejects_non_finite_initial_objective() {
        let error = brent_minimize(
            |x| if x > 2.2 { f64::NAN } else { (x - 1.5).powi(2) },
            1.0,
            5.0,
        )
        .unwrap_err();
        match error {
            OptimizationError::NonFiniteObjective { x, fx } => {
                assert_close(x, 2.527_864_045, 1e-9);
                assert!(fx.is_nan());
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn o7_bounds_errors_precede_objective_calls() {
        let calls = Cell::new(0);
        assert_eq!(
            brent_minimize(
                |_| {
                    calls.set(calls.get() + 1);
                    0.0
                },
                5.0,
                1.0
            ),
            Err(OptimizationError::InvalidBounds {
                lower: 5.0,
                upper: 1.0
            })
        );
        assert!(
            matches!(brent_minimize(|_| 0.0, 1.0, f64::NAN), Err(OptimizationError::NonFiniteBounds(x)) if x.is_nan())
        );
        assert_eq!(
            brent_minimize(|_| 0.0, f64::INFINITY, 5.0),
            Err(OptimizationError::NonFiniteBounds(f64::INFINITY))
        );
        assert_eq!(calls.get(), 0);
    }

    #[test]
    fn o8_budget_exhaustion_reports_exact_count() {
        let error =
            brent_minimize_with_budget(|x| (x - 2.5).powi(2), 1.0, 5.0, 1e-5, 3).unwrap_err();
        match error {
            OptimizationError::MaxEvaluationsExceeded {
                best_x,
                evaluations,
                ..
            } => {
                assert_close(best_x, 2.527_864_045, 1e-9);
                assert_eq!(evaluations, 3);
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn o9_tighter_tolerance_has_monotonic_evaluation_counts() {
        let counts: Vec<_> = [1e-2, 1e-4, 1e-6, 1e-8]
            .into_iter()
            .map(|xatol| {
                brent_minimize_with_budget(
                    |x| (x - std::f64::consts::E).powi(4) + (x - std::f64::consts::E).powi(2),
                    1.0,
                    5.0,
                    xatol,
                    500,
                )
                .unwrap()
                .func_evals
            })
            .collect();
        assert_eq!(counts, [7, 9, 10, 10]);
    }

    #[test]
    fn o10_is_bitwise_deterministic() {
        let expected = brent_minimize(|x| (x - 2.5).powi(2), 1.0, 5.0).unwrap();
        assert_eq!(expected.func_evals, 6);
        for _ in 0..100 {
            let result = brent_minimize(|x| (x - 2.5).powi(2), 1.0, 5.0).unwrap();
            assert_eq!(result.minimizer.to_bits(), expected.minimizer.to_bits());
            assert_eq!(result.func_evals, 6);
        }
    }

    #[test]
    fn d1_equal_bounds_bypass_with_one_evaluation() {
        let calls = Cell::new(0);
        let result = brent_minimize(
            |x| {
                calls.set(calls.get() + 1);
                x * x
            },
            1.0,
            1.0,
        )
        .unwrap();
        assert_eq!(calls.get(), 1);
        assert_eq!(result.minimizer, 1.0);
        assert_eq!(result.func_evals, 1);
        assert_eq!(result.iterations, 0);
        assert_eq!(result.stop_reason, StopReason::DegenerateIntervalBypass);
    }

    #[test]
    fn d2_equal_bounds_preserves_non_finite_objective_error() {
        let calls = Cell::new(0);
        let error = brent_minimize(
            |_| {
                calls.set(calls.get() + 1);
                f64::NAN
            },
            1.0,
            1.0,
        )
        .unwrap_err();
        assert_eq!(calls.get(), 1);
        assert!(
            matches!(error, OptimizationError::NonFiniteObjective { x: 1.0, fx } if fx.is_nan())
        );
    }

    fn extreme_trace(lower: f64, upper: f64) -> (OptimizationDiagnostics, Vec<u64>) {
        let calls = Cell::new(0);
        let points = RefCell::new(Vec::new());
        let result = brent_minimize(
            |x| {
                calls.set(calls.get() + 1);
                assert!(x.is_finite() && x >= lower && x <= upper);
                points.borrow_mut().push(x.to_bits());
                1.0
            },
            lower,
            upper,
        )
        .unwrap();
        assert!(calls.get() > 0);
        assert!(result.minimizer.is_finite());
        assert!(result.final_bracket.0 >= lower && result.final_bracket.1 <= upper);
        (result, points.into_inner())
    }

    fn assert_extreme_candidates(lower: f64, upper: f64) {
        let (first, first_trace) = extreme_trace(lower, upper);
        let (second, second_trace) = extreme_trace(lower, upper);
        assert_eq!(first.minimizer.to_bits(), second.minimizer.to_bits());
        assert_eq!(first_trace, second_trace);
    }

    #[test]
    fn e1_extreme_positive_bounds_remain_finite() {
        assert_extreme_candidates(1.0e308, 1.7e308);
    }

    #[test]
    fn e2_extreme_negative_bounds_remain_finite() {
        assert_extreme_candidates(-1.7e308, -1.0e308);
    }

    #[test]
    fn e3_extreme_opposite_sign_bounds_use_safe_arithmetic() {
        assert_extreme_candidates(-1.0e308, 1.0e308);
    }

    #[test]
    fn e4_near_maximum_candidates_are_clamped_and_finite() {
        assert_extreme_candidates(f64::MAX * 0.9, f64::MAX);
    }

    #[test]
    #[should_panic(expected = "xatol must be finite")]
    fn p1_nan_xatol_panics() {
        let _ = brent_minimize_with_budget(|_| 0.0, 1.0, 2.0, f64::NAN, 1);
    }

    #[test]
    #[should_panic(expected = "xatol must be finite")]
    fn p1_infinite_xatol_panics() {
        let _ = brent_minimize_with_budget(|_| 0.0, 1.0, 2.0, f64::INFINITY, 1);
    }

    #[test]
    #[should_panic(expected = "xatol must be finite")]
    fn p1_negative_infinite_xatol_panics() {
        let _ = brent_minimize_with_budget(|_| 0.0, 1.0, 2.0, f64::NEG_INFINITY, 1);
    }

    #[test]
    #[should_panic(expected = "xatol must be strictly positive")]
    fn p2_zero_xatol_panics() {
        let _ = brent_minimize_with_budget(|_| 0.0, 1.0, 2.0, 0.0, 1);
    }

    #[test]
    #[should_panic(expected = "xatol must be strictly positive")]
    fn p3_negative_xatol_panics() {
        let _ = brent_minimize_with_budget(|_| 0.0, 1.0, 2.0, -1e-5, 1);
    }

    #[test]
    #[should_panic(expected = "maxfun must be strictly positive")]
    fn p4_zero_maxfun_panics() {
        let _ = brent_minimize_with_budget(|_| 0.0, 1.0, 2.0, 1e-5, 0);
    }
}
