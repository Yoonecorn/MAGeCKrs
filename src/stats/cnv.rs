//! Deterministic in-memory copy-number model and score correction.

use std::fmt;

use crate::error::MageckError;
use crate::stats::brent::{brent_minimize, OptimizationDiagnostics, OptimizationError};
use crate::stats::type7_percentile_999;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CnvModel {
    pub(crate) breakpoint: f64,
    pub(crate) slope: f64,
    pub(crate) intercept: f64,
    pub(crate) sse: f64,
    pub(crate) mu_ref: f64,
    pub(crate) upper_bound: f64,
    pub(crate) diagnostics: OptimizationDiagnostics,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum CnvError {
    LengthMismatch { scores_len: usize, cn_len: usize },
    EmptyPopulation,
    NonFiniteScore { index: usize, value: f64 },
    NonFiniteCopyNumber { index: usize, value: f64 },
    NegativeCopyNumber { index: usize, value: f64 },
    InsufficientPopulation { count: usize },
    RankDeficientPopulation { variance: f64 },
    InsufficientLowSideObservations { bp: f64, count: usize },
    RankDeficientLowSide { bp: f64, variance: f64 },
    InvalidBounds { lower: f64, upper: f64 },
    OptimizationFailure(OptimizationError),
    NonFiniteLog2CopyNumber { index: usize, value: f64 },
    Log2TransformationOverflow { index: usize, value: f64 },
    PercentileCalculationFailure(String),
    NonFiniteCalculation(String),
}

impl fmt::Display for CnvError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LengthMismatch { scores_len, cn_len } => write!(
                f,
                "CNV scores length {scores_len} differs from copy-number length {cn_len}"
            ),
            Self::EmptyPopulation => write!(f, "CNV population is empty"),
            Self::NonFiniteScore { index, value } => {
                write!(f, "Non-finite CNV score at index {index}: {value}")
            }
            Self::NonFiniteCopyNumber { index, value } => {
                write!(f, "Non-finite copy number at index {index}: {value}")
            }
            Self::NegativeCopyNumber { index, value } => {
                write!(f, "Negative copy number at index {index}: {value}")
            }
            Self::InsufficientPopulation { count } => write!(
                f,
                "CNV population requires at least two observations, got {count}"
            ),
            Self::RankDeficientPopulation { variance } => {
                write!(f, "Rank-deficient CNV population (Sxx={variance})")
            }
            Self::InsufficientLowSideObservations { bp, count } => write!(
                f,
                "CNV breakpoint {bp} has only {count} low-side observations"
            ),
            Self::RankDeficientLowSide { bp, variance } => write!(
                f,
                "Rank-deficient CNV low side at breakpoint {bp} (Sxx={variance})"
            ),
            Self::InvalidBounds { lower, upper } => {
                write!(f, "Invalid CNV breakpoint bounds [{lower}, {upper}]")
            }
            Self::OptimizationFailure(error) => write!(f, "CNV optimization failed: {error:?}"),
            Self::NonFiniteLog2CopyNumber { index, value } => {
                write!(f, "Non-finite log2 copy number at index {index}: {value}")
            }
            Self::Log2TransformationOverflow { index, value } => write!(
                f,
                "Log2 copy-number transformation overflow at index {index}: {value}"
            ),
            Self::PercentileCalculationFailure(message) => {
                write!(f, "CNV percentile calculation failed: {message}")
            }
            Self::NonFiniteCalculation(message) => write!(f, "CNV calculation failed: {message}"),
        }
    }
}

impl std::error::Error for CnvError {}

pub(crate) fn transform_log2_cn(log2_cn: &[f64]) -> Result<Vec<f64>, CnvError> {
    let mut result = Vec::with_capacity(log2_cn.len());
    for (index, &value) in log2_cn.iter().enumerate() {
        if !value.is_finite() {
            return Err(CnvError::NonFiniteLog2CopyNumber { index, value });
        }
        let cn = value.exp2();
        if !cn.is_finite() {
            return Err(CnvError::Log2TransformationOverflow { index, value });
        }
        result.push(positive_zero(cn));
    }
    Ok(result)
}

fn positive_zero(value: f64) -> f64 {
    if value == 0.0 {
        0.0
    } else {
        value
    }
}

fn finite(value: f64, message: &'static str) -> Result<f64, CnvError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CnvError::NonFiniteCalculation(message.to_string()))
    }
}

// Two passes over the selected observations, always in original source order.
fn low_side_ols(scores: &[f64], cn: &[f64], bp: f64) -> Result<(f64, f64), CnvError> {
    let mut count = 0;
    let (mut sum_x, mut sum_y) = (0.0, 0.0);
    for (&y, &x) in scores.iter().zip(cn) {
        if x <= bp {
            count += 1;
            sum_x += x;
            sum_y += y;
        }
    }
    if count < 2 {
        return Err(CnvError::InsufficientLowSideObservations { bp, count });
    }
    let mean_x = finite(sum_x / count as f64, "Non-finite low-side mean")?;
    let mean_y = finite(sum_y / count as f64, "Non-finite low-side mean")?;
    let (mut s_xx, mut s_xy) = (0.0, 0.0);
    for (&y, &x) in scores.iter().zip(cn) {
        if x <= bp {
            s_xx += (x - mean_x) * (x - mean_x);
            s_xy += (x - mean_x) * (y - mean_y);
        }
    }
    if s_xx == 0.0 {
        return Err(CnvError::RankDeficientLowSide { bp, variance: 0.0 });
    }
    finite(s_xx, "Non-finite low-side variance")?;
    finite(s_xy, "Non-finite low-side covariance")?;
    let slope = s_xy / s_xx;
    let intercept = mean_y - slope * mean_x;
    if !slope.is_finite() || !intercept.is_finite() {
        return Err(CnvError::NonFiniteCalculation(
            "Non-finite OLS coefficients at lower bound bp=1.0".to_string(),
        ));
    }
    Ok((positive_zero(slope), positive_zero(intercept)))
}

fn objective(scores: &[f64], cn: &[f64], bp: f64) -> Result<(f64, f64, f64), CnvError> {
    let (slope, intercept) = low_side_ols(scores, cn, bp)?;
    let mut sse = 0.0;
    for (&y, &x) in scores.iter().zip(cn) {
        let prediction = if x <= bp {
            intercept + slope * x
        } else {
            intercept + slope * bp
        };
        let diff = y - prediction;
        sse += diff * diff;
    }
    let sse = finite(sse, "Non-finite SSE encountered")?;
    Ok((positive_zero(sse), slope, intercept))
}

pub(crate) fn fit_cnv_model(scores: &[f64], cn_values: &[f64]) -> Result<CnvModel, CnvError> {
    if scores.len() != cn_values.len() {
        return Err(CnvError::LengthMismatch {
            scores_len: scores.len(),
            cn_len: cn_values.len(),
        });
    }
    if scores.is_empty() {
        return Err(CnvError::EmptyPopulation);
    }
    for (index, &value) in scores.iter().enumerate() {
        if !value.is_finite() {
            return Err(CnvError::NonFiniteScore { index, value });
        }
    }
    for (index, &value) in cn_values.iter().enumerate() {
        if !value.is_finite() {
            return Err(CnvError::NonFiniteCopyNumber { index, value });
        }
    }
    for (index, &value) in cn_values.iter().enumerate() {
        if value < 0.0 {
            return Err(CnvError::NegativeCopyNumber { index, value });
        }
    }
    if scores.len() < 2 {
        return Err(CnvError::InsufficientPopulation {
            count: scores.len(),
        });
    }

    let mean_global = finite(
        cn_values.iter().sum::<f64>() / cn_values.len() as f64,
        "Non-finite population mean",
    )?;
    let s_xx_global = cn_values
        .iter()
        .map(|&x| (x - mean_global) * (x - mean_global))
        .sum::<f64>();
    if s_xx_global == 0.0 {
        return Err(CnvError::RankDeficientPopulation { variance: 0.0 });
    }
    finite(s_xx_global, "Non-finite population variance")?;
    low_side_ols(scores, cn_values, 1.0)?;

    let upper_bound = type7_percentile_999(cn_values).map_err(|error| match error {
        MageckError::Math(message) => CnvError::PercentileCalculationFailure(message),
        other => CnvError::PercentileCalculationFailure(other.to_string()),
    })?;
    if upper_bound < 1.0 {
        return Err(CnvError::InvalidBounds {
            lower: 1.0,
            upper: upper_bound,
        });
    }

    let mut calculation_error = None;
    let mut evaluate = |bp| match objective(scores, cn_values, bp) {
        Ok((sse, _, _)) => sse,
        Err(error) => {
            calculation_error = Some(error);
            f64::NAN
        }
    };
    let optimization = brent_minimize(&mut evaluate, 1.0, upper_bound);
    if let Some(error) = calculation_error {
        return Err(error);
    }
    let diagnostics = optimization.map_err(CnvError::OptimizationFailure)?;
    let breakpoint = diagnostics.minimizer;
    let (sse, slope, intercept) = objective(scores, cn_values, breakpoint)?;
    let mu_ref = finite(intercept + slope * 1.0, "Non-finite diploid reference")?;
    Ok(CnvModel {
        breakpoint,
        slope,
        intercept,
        sse,
        mu_ref: positive_zero(mu_ref),
        upper_bound,
        diagnostics,
    })
}

impl CnvModel {
    pub(crate) fn predict(&self, cn: f64) -> Result<f64, CnvError> {
        if !cn.is_finite() {
            return Err(CnvError::NonFiniteCopyNumber {
                index: 0,
                value: cn,
            });
        }
        if cn < 0.0 {
            return Err(CnvError::NegativeCopyNumber {
                index: 0,
                value: cn,
            });
        }
        let x = if cn < self.breakpoint {
            cn
        } else {
            self.breakpoint
        };
        finite(self.intercept + self.slope * x, "Non-finite CNV prediction")
    }

    pub(crate) fn adjustment(&self, cn: f64) -> Result<f64, CnvError> {
        let result = finite(self.mu_ref - self.predict(cn)?, "Non-finite CNV adjustment")?;
        Ok(positive_zero(result))
    }

    pub(crate) fn correct(&self, score: f64, cn: f64) -> Result<f64, CnvError> {
        if !score.is_finite() {
            return Err(CnvError::NonFiniteScore {
                index: 0,
                value: score,
            });
        }
        finite(score + self.adjustment(cn)?, "Non-finite corrected score")
    }

    pub(crate) fn correct_slice(
        &self,
        scores: &[f64],
        cn_values: &[f64],
    ) -> Result<Vec<f64>, CnvError> {
        if scores.len() != cn_values.len() {
            return Err(CnvError::LengthMismatch {
                scores_len: scores.len(),
                cn_len: cn_values.len(),
            });
        }
        scores
            .iter()
            .zip(cn_values)
            .enumerate()
            .map(|(index, (&score, &cn))| {
                self.correct(score, cn).map_err(|error| match error {
                    CnvError::NonFiniteScore { value, .. } => {
                        CnvError::NonFiniteScore { index, value }
                    }
                    CnvError::NonFiniteCopyNumber { value, .. } => {
                        CnvError::NonFiniteCopyNumber { index, value }
                    }
                    CnvError::NegativeCopyNumber { value, .. } => {
                        CnvError::NegativeCopyNumber { index, value }
                    }
                    other => other,
                })
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::brent::{StopReason, DEFAULT_XATOL};

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual} != {expected}"
        );
    }

    fn fixture() -> ([f64; 10], [f64; 10]) {
        let cn = transform_log2_cn(&[-1.0, 0.0, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.2]).unwrap();
        let scores = [-1.5, -0.8, -0.3, 0.2, 0.5, 0.9, 1.4, 1.8, 2.0, 2.0];
        (cn.try_into().unwrap(), scores)
    }

    #[test]
    fn validation_order_and_domains() {
        assert!(matches!(
            fit_cnv_model(&[f64::NAN], &[]),
            Err(CnvError::LengthMismatch {
                scores_len: 1,
                cn_len: 0
            })
        ));
        assert_eq!(fit_cnv_model(&[], &[]), Err(CnvError::EmptyPopulation));
        match fit_cnv_model(&[1.0, f64::NAN], &[f64::INFINITY, -1.0]) {
            Err(CnvError::NonFiniteScore { index: 1, value }) => assert!(value.is_nan()),
            other => panic!("{other:?}"),
        }
        match fit_cnv_model(&[1.0, 2.0], &[f64::NAN, -1.0]) {
            Err(CnvError::NonFiniteCopyNumber { index: 0, value }) => assert!(value.is_nan()),
            other => panic!("{other:?}"),
        }
        assert!(
            matches!(fit_cnv_model(&[1.0, 2.0], &[1.0, f64::INFINITY]), Err(CnvError::NonFiniteCopyNumber { index: 1, value }) if value == f64::INFINITY)
        );
        assert!(matches!(
            fit_cnv_model(&[1.0, 2.0], &[1.0, -0.5]),
            Err(CnvError::NegativeCopyNumber {
                index: 1,
                value: -0.5
            })
        ));
        assert_eq!(
            fit_cnv_model(&[1.0], &[1.0]),
            Err(CnvError::InsufficientPopulation { count: 1 })
        );
        assert!(fit_cnv_model(&[-2.5, -0.5, 1.2], &[0.5, 1.0, 2.0]).is_ok());
    }

    #[test]
    fn rank_and_lower_bound_preflight() {
        assert_eq!(
            fit_cnv_model(&[1.0, 2.0, 3.0], &[2.0; 3]),
            Err(CnvError::RankDeficientPopulation { variance: 0.0 })
        );
        assert_eq!(
            fit_cnv_model(&[1.0, 2.0, 3.0], &[0.5, 1.5, 2.0]),
            Err(CnvError::InsufficientLowSideObservations { bp: 1.0, count: 1 })
        );
        assert_eq!(
            fit_cnv_model(&[1.0, 2.0, 3.0], &[1.0, 1.0, 2.0]),
            Err(CnvError::RankDeficientLowSide {
                bp: 1.0,
                variance: 0.0
            })
        );
        if let Err(CnvError::RankDeficientLowSide { variance, .. }) =
            fit_cnv_model(&[1.0, 2.0, 3.0], &[1.0, 1.0, 2.0])
        {
            assert_eq!(variance.to_bits(), 0.0_f64.to_bits());
        }
    }

    #[test]
    fn log2_transform_boundaries_and_precedence() {
        assert_eq!(transform_log2_cn(&[]), Ok(vec![]));
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            match transform_log2_cn(&[1.0, value, 1024.0]) {
                Err(CnvError::NonFiniteLog2CopyNumber {
                    index: 1,
                    value: found,
                }) if value.is_nan() => assert!(found.is_nan()),
                Err(CnvError::NonFiniteLog2CopyNumber {
                    index: 1,
                    value: found,
                }) => assert_eq!(found, value),
                other => panic!("{other:?}"),
            }
        }
        assert_eq!(
            transform_log2_cn(&[1024.0]),
            Err(CnvError::Log2TransformationOverflow {
                index: 0,
                value: 1024.0
            })
        );
        assert!(transform_log2_cn(&[1023.0]).unwrap()[0].is_finite());
        assert_eq!(transform_log2_cn(&[0.0, -0.0]).unwrap(), vec![1.0, 1.0]);
        assert_eq!(
            transform_log2_cn(&[-1100.0]).unwrap()[0].to_bits(),
            0.0_f64.to_bits()
        );
    }

    #[test]
    fn interval_and_boundary_cases() {
        assert!(
            matches!(fit_cnv_model(&[1.0; 4], &[0.2, 0.4, 0.6, 0.8]), Err(CnvError::InvalidBounds { lower: 1.0, upper }) if (upper - 0.7994).abs() < 1e-12)
        );
        let model = fit_cnv_model(&[-1.0, 0.0, 0.5], &[0.5, 1.0, 1.0]).unwrap();
        assert_eq!(model.breakpoint, 1.0);
        assert_eq!(model.diagnostics.iterations, 0);
        assert_eq!(model.diagnostics.func_evals, 1);
        assert_eq!(model.diagnostics.final_bracket, (1.0, 1.0));
        assert_eq!(
            model.diagnostics.stop_reason,
            StopReason::DegenerateIntervalBypass
        );

        let linear = fit_cnv_model(
            &[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            &[0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        )
        .unwrap();
        close(linear.upper_bound, 2.9975, 1e-6);
        close(linear.breakpoint, 2.9975, 3.0 * DEFAULT_XATOL);
        close(linear.slope, 2.0, 1e-6);
        close(linear.intercept, -1.0, 1e-6);
        close(linear.sse, 2.5e-5, 1e-6);
        assert!(linear.diagnostics.final_bracket.0 <= linear.breakpoint);
        assert!(linear.diagnostics.final_bracket.1 >= linear.breakpoint);
    }

    #[test]
    fn known_piecewise_and_flat_scores() {
        let model = fit_cnv_model(
            &[0.25, 1.0, 1.75, 2.5, 2.8, 2.8, 2.8, 2.8],
            &[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
        )
        .unwrap();
        close(model.breakpoint, 2.2, 1e-5);
        close(model.slope, 1.5, 1e-5);
        close(model.intercept, -0.5, 1e-5);
        assert!(model.sse <= 1e-10);
        close(model.adjustment(2.5).unwrap(), -1.8, 1e-5);
        close(model.correct(2.8, 2.5).unwrap(), 1.0, 1e-5);
        assert_eq!(model.adjustment(1.0).unwrap().to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            model.predict(model.breakpoint + 0.5).unwrap().to_bits(),
            model.predict(model.breakpoint + 1.0).unwrap().to_bits()
        );

        let flat = fit_cnv_model(&[5.0; 5], &[0.5, 1.0, 1.5, 2.0, 2.5]).unwrap();
        close(flat.breakpoint, 2.498, 3.0 * DEFAULT_XATOL);
        assert_eq!(flat.slope.to_bits(), 0.0_f64.to_bits());
        assert_eq!(flat.sse.to_bits(), 0.0_f64.to_bits());
        assert_eq!(flat.correct(5.0, 2.0).unwrap(), 5.0);
    }

    #[test]
    fn repaired_direct_ten_point_oracle_and_determinism() {
        let (cn, scores) = fixture();
        assert_eq!(cn.iter().filter(|&&x| x <= 1.0).count(), 2);
        close(cn[9], 4.59479341998814, 1e-12);
        let model = fit_cnv_model(&scores, &cn).unwrap();
        close(model.upper_bound, 4.5894402792, 1e-6);
        close(model.breakpoint, 3.1637400632, 1e-5);
        close(model.slope, 1.2670796159, 1e-5);
        close(model.intercept, -2.0753772108, 1e-5);
        close(model.sse, 0.0527995599, 1e-5);
        close(model.mu_ref, -0.8082975949, 1e-5);
        close(model.diagnostics.minimizer, 3.1637400632, 1e-5);
        close(model.diagnostics.objective_value, 0.0527995599, 1e-5);
        assert_eq!(model.diagnostics.iterations, 7);
        assert_eq!(model.diagnostics.func_evals, 8);
        assert!(model.diagnostics.converged);
        assert_eq!(
            model.diagnostics.stop_reason,
            StopReason::IntervalToleranceMet
        );
        close(model.diagnostics.final_bracket.0, 3.1637366827, 1e-5);
        close(model.diagnostics.final_bracket.1, 3.1637434437, 1e-5);
        let expected = [
            -0.8664601920,
            -0.8,
            -0.8248415615,
            -0.7390341308,
            -0.7670796159,
            -0.7439049250,
            -0.9167627390,
            -0.9416309282,
            -0.7416309282,
            -0.7416309282,
        ];
        let corrected = model.correct_slice(&scores, &cn).unwrap();
        for (actual, expected) in corrected.iter().zip(expected) {
            close(*actual, expected, 1e-5);
        }
        for _ in 0..100 {
            let again = fit_cnv_model(&scores, &cn).unwrap();
            assert_eq!(again.upper_bound.to_bits(), model.upper_bound.to_bits());
            assert_eq!(again.breakpoint.to_bits(), model.breakpoint.to_bits());
            assert_eq!(again.slope.to_bits(), model.slope.to_bits());
            assert_eq!(again.intercept.to_bits(), model.intercept.to_bits());
            assert_eq!(again.sse.to_bits(), model.sse.to_bits());
            assert_eq!(again.mu_ref.to_bits(), model.mu_ref.to_bits());
            assert_eq!(
                again.diagnostics.minimizer.to_bits(),
                model.diagnostics.minimizer.to_bits()
            );
            assert_eq!(
                again.diagnostics.objective_value.to_bits(),
                model.diagnostics.objective_value.to_bits()
            );
            assert_eq!(
                again.diagnostics.final_bracket.0.to_bits(),
                model.diagnostics.final_bracket.0.to_bits()
            );
            assert_eq!(
                again.diagnostics.final_bracket.1.to_bits(),
                model.diagnostics.final_bracket.1.to_bits()
            );
            assert_eq!(again.diagnostics.iterations, model.diagnostics.iterations);
            assert_eq!(again.diagnostics.func_evals, model.diagnostics.func_evals);
        }
    }

    #[test]
    fn correction_validation_and_top_level_mapping() {
        let (cn, scores) = fixture();
        let model = fit_cnv_model(&scores, &cn).unwrap();
        assert!(matches!(
            model.correct_slice(&[1.0], &[]),
            Err(CnvError::LengthMismatch { .. })
        ));
        assert!(matches!(
            model.correct(f64::INFINITY, 1.0),
            Err(CnvError::NonFiniteScore { .. })
        ));
        assert!(matches!(
            model.adjustment(-1.0),
            Err(CnvError::NegativeCopyNumber { .. })
        ));
        match model.predict(f64::NAN) {
            Err(CnvError::NonFiniteCopyNumber { index: 0, value }) => assert!(value.is_nan()),
            other => panic!("{other:?}"),
        }
        let error: MageckError = CnvError::EmptyPopulation.into();
        assert!(
            matches!(error, MageckError::Math(ref message) if message == "CNV population is empty")
        );
    }

    #[test]
    fn correction_keeps_input_order_and_reports_original_index() {
        let (cn, scores) = fixture();
        let model = fit_cnv_model(&scores, &cn).unwrap();
        let reordered_cn = [cn[9], cn[1], cn[0], cn[7]];
        let reordered_scores = [scores[9], scores[1], scores[0], scores[7]];
        let corrected = model
            .correct_slice(&reordered_scores, &reordered_cn)
            .unwrap();
        for index in 0..corrected.len() {
            assert_eq!(
                corrected[index].to_bits(),
                model
                    .correct(reordered_scores[index], reordered_cn[index])
                    .unwrap()
                    .to_bits()
            );
        }
        match model.correct_slice(&[1.0, f64::NAN], &[1.0, 1.0]) {
            Err(CnvError::NonFiniteScore { index: 1, value }) => assert!(value.is_nan()),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn finite_inputs_with_nonfinite_arithmetic_are_rejected() {
        assert!(matches!(
            fit_cnv_model(&[1.0, 2.0, 3.0], &[0.0, 1.0, f64::MAX]),
            Err(CnvError::NonFiniteCalculation(_))
        ));
    }
}
