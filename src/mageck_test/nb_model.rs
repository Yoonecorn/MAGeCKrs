use crate::data::count_matrix::CountMatrix;
use crate::error::MageckError;

pub const POISSON_DISPERSION: f64 = 1.0e-4;
pub const POISSON_SIZE: f64 = 10_000.0;

/// Optional regression coefficients from a log-log OLS or WLS fit.
#[derive(Debug, Clone, PartialEq)]
pub struct RegressionCoefficients {
    pub intercept: f64,
    pub slope: f64,
}

/// A deterministic Negative Binomial mean-variance model.
#[derive(Debug, Clone, PartialEq)]
pub struct NbModel {
    mean: Vec<f64>,
    variance: Vec<f64>,
    raw_dispersion: Vec<Option<f64>>,
    fitted_dispersion: Vec<f64>,
    final_dispersion: Vec<f64>,
    size: Vec<f64>,
    eligible_guides_count: usize,
    regression_coefficients: Option<RegressionCoefficients>,
    median_dispersion: f64,
}

impl NbModel {
    /// Fits dispersion from the selected control and treatment columns while
    /// retaining the single control column as the comparison mean.
    ///
    /// Pooled ddof=1 moments determine the WLS trend, weighted by pooled mean.
    /// The trend is evaluated at the control mean with 100% shrinkage to fitted
    /// dispersion. Multi-control `fit_meanvar` retains its conservative maximum.
    pub fn fit_pooled_dispersion(
        counts: &CountMatrix,
        control_index: usize,
        pooled_indices: &[usize],
    ) -> Result<Self, MageckError> {
        if !pooled_indices.contains(&control_index) {
            return Err(MageckError::Math(
                "Pooled dispersion samples must include the control".into(),
            ));
        }
        let mut model = fit_meanvar(counts, pooled_indices)?;
        // Reuse the validated pooled moments, eligibility and median fallback;
        // leave the multi-control OLS implementation untouched.
        model.regression_coefficients = None;
        if model.eligible_guides_count >= 3 {
            let mut sum_w = 0.0;
            let mut sum_wx = 0.0;
            let mut sum_wy = 0.0;
            let mut sum_wxx = 0.0;
            for (g, raw) in model.raw_dispersion.iter().enumerate() {
                if let Some(raw) = raw {
                    let w = model.mean[g];
                    let x = w.ln();
                    let y = raw.ln();
                    sum_w += w;
                    sum_wx += w * x;
                    sum_wy += w * y;
                    sum_wxx += w * x * x;
                }
            }
            let xbar = sum_wx / sum_w;
            let ybar = sum_wy / sum_w;
            // Center before accumulating covariance: the uncentered normal
            // equations cancel large terms when eligible means are close.
            let mut denominator = 0.0;
            let mut numerator = 0.0;
            for (g, raw) in model.raw_dispersion.iter().enumerate() {
                if let Some(raw) = raw {
                    let w = model.mean[g];
                    let dx = w.ln() - xbar;
                    denominator += w * dx * dx;
                    numerator += w * dx * (raw.ln() - ybar);
                }
            }
            // The original identifiability guard divided by sum_w. Equal
            // pooled means still use the accepted median fallback.
            let scale = sum_wxx.abs().max((sum_wx * xbar).abs());
            if sum_w.is_finite()
                && sum_w > 0.0
                && denominator.is_finite()
                && denominator > f64::EPSILON * scale
            {
                let slope = numerator / denominator;
                let intercept = ybar - slope * xbar;
                if slope.is_finite() && intercept.is_finite() {
                    model.regression_coefficients =
                        Some(RegressionCoefficients { intercept, slope });
                }
            }
        }
        // fit_meanvar validates the control before this access. Comparison means
        // and p-value inputs stay control-only, including zero-count sentinels.
        model.mean = counts.sample_slice(control_index).to_vec();
        for (g, &control_mean) in model.mean.iter().enumerate() {
            let mut fitted = if control_mean == 0.0 {
                POISSON_DISPERSION
            } else {
                model.median_dispersion
            };
            if control_mean > 0.0 {
                if let Some(coeff) = &model.regression_coefficients {
                    let pred = (coeff.intercept + coeff.slope * control_mean.ln()).exp();
                    if pred.is_finite() && pred > 0.0 {
                        fitted = pred;
                    }
                }
            }
            let size = 1.0 / fitted;
            if !size.is_finite() || size <= 0.0 {
                return Err(MageckError::Math(
                    "Invalid size parameter calculated".into(),
                ));
            }
            model.fitted_dispersion[g] = fitted;
            model.final_dispersion[g] = fitted;
            model.size[g] = size;
        }
        Ok(model)
    }

    /// Ordered per-guide arithmetic means across selected control samples.
    #[inline]
    pub fn mean(&self) -> &[f64] {
        &self.mean
    }

    /// Ordered per-guide unbiased variances (denominator n - 1).
    #[inline]
    pub fn variance(&self) -> &[f64] {
        &self.variance
    }

    /// Ordered per-guide raw dispersions, None if ineligible.
    #[inline]
    pub fn raw_dispersion(&self) -> &[Option<f64>] {
        &self.raw_dispersion
    }

    /// Ordered per-guide smoothed/fallback dispersions.
    #[inline]
    pub fn fitted_dispersion(&self) -> &[f64] {
        &self.fitted_dispersion
    }

    /// Ordered final dispersions: conservative for multi-control, fitted-only for pooled.
    #[inline]
    pub fn final_dispersion(&self) -> &[f64] {
        &self.final_dispersion
    }

    /// Ordered per-guide Negative Binomial size parameters (1 / final_dispersion).
    #[inline]
    pub fn size(&self) -> &[f64] {
        &self.size
    }

    /// Number of guides eligible for regression and median calculations.
    #[inline]
    pub fn eligible_guides_count(&self) -> usize {
        self.eligible_guides_count
    }

    /// Regression coefficients if a valid OLS or pooled WLS fit was found.
    #[inline]
    pub fn regression_coefficients(&self) -> Option<&RegressionCoefficients> {
        self.regression_coefficients.as_ref()
    }

    /// Median dispersion used as fallback.
    #[inline]
    pub fn median_dispersion(&self) -> f64 {
        self.median_dispersion
    }
}

/// Fits a Negative Binomial mean-variance model using only the specified control samples.
pub fn fit_meanvar(
    counts: &CountMatrix,
    control_indices: &[usize],
) -> Result<NbModel, MageckError> {
    if counts.is_empty() {
        return Err(MageckError::Math("CountMatrix is empty".into()));
    }
    if counts.n_guides() == 0 || counts.n_samples() == 0 {
        return Err(MageckError::Math(
            "CountMatrix dimensions cannot be zero".into(),
        ));
    }
    if control_indices.is_empty() {
        return Err(MageckError::Math("Control list is empty".into()));
    }
    if control_indices.len() < 2 {
        return Err(MageckError::Math(
            "At least two control samples are required for unbiased variance".into(),
        ));
    }

    let mut sorted_indices = control_indices.to_vec();
    sorted_indices.sort_unstable();
    for i in 1..sorted_indices.len() {
        if sorted_indices[i] == sorted_indices[i - 1] {
            return Err(MageckError::Math(
                "Duplicate control indices are not allowed".into(),
            ));
        }
    }
    if *sorted_indices.last().unwrap() >= counts.n_samples() {
        return Err(MageckError::Math("Control index out of bounds".into()));
    }

    let n_guides = counts.n_guides();
    let n_control = control_indices.len();

    let mut mean = Vec::with_capacity(n_guides);
    let mut variance = Vec::with_capacity(n_guides);
    let mut raw_dispersion = Vec::with_capacity(n_guides);

    let mut eligible_raw_disps = Vec::new();
    let mut eligible_xs = Vec::new();
    let mut eligible_ys = Vec::new();

    for g in 0..n_guides {
        let mut sum = 0.0;
        for &s in control_indices {
            let val = counts.entry(s, g);
            if !val.is_finite() {
                return Err(MageckError::Math(
                    "Non-finite selected control count encountered".into(),
                ));
            }
            if val < 0.0 {
                return Err(MageckError::Math(
                    "Negative selected control count encountered".into(),
                ));
            }
            sum += val;
        }
        let m = sum / (n_control as f64);
        if !m.is_finite() {
            return Err(MageckError::Math("Non-finite mean calculated".into()));
        }

        let mut sq_sum = 0.0;
        for &s in control_indices {
            let val = counts.entry(s, g);
            let diff = val - m;
            sq_sum += diff * diff;
        }
        let v = sq_sum / ((n_control - 1) as f64);
        if !v.is_finite() {
            return Err(MageckError::Math("Non-finite variance calculated".into()));
        }

        mean.push(m);
        variance.push(v);

        let raw = if m > 0.0 && v > m {
            let raw_val = (v - m) / (m * m);
            if raw_val.is_finite() && raw_val > 0.0 {
                Some(raw_val)
            } else {
                None
            }
        } else {
            None
        };
        raw_dispersion.push(raw);

        if let Some(r) = raw {
            eligible_raw_disps.push(r);
            eligible_xs.push(m.ln());
            eligible_ys.push(r.ln());
        }
    }

    let median_dispersion = if eligible_raw_disps.is_empty() {
        POISSON_DISPERSION
    } else {
        eligible_raw_disps.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let len = eligible_raw_disps.len();
        let med = if len % 2 == 1 {
            eligible_raw_disps[len / 2]
        } else {
            (eligible_raw_disps[len / 2 - 1] + eligible_raw_disps[len / 2]) / 2.0
        };
        if med.is_finite() && med > 0.0 {
            med
        } else {
            return Err(MageckError::Math(
                "Invalid median dispersion calculated".into(),
            ));
        }
    };

    let eligible_guides_count = eligible_xs.len();
    let mut regression_coefficients = None;

    if eligible_guides_count >= 3 {
        let n_eligible = eligible_guides_count as f64;
        let sum_x: f64 = eligible_xs.iter().sum();
        let sum_y: f64 = eligible_ys.iter().sum();

        if sum_x.is_finite() && sum_y.is_finite() {
            let mean_x = sum_x / n_eligible;
            let mean_y = sum_y / n_eligible;

            let mut ss_x = 0.0;
            let mut cov_xy = 0.0;

            for i in 0..eligible_guides_count {
                let x_diff = eligible_xs[i] - mean_x;
                let y_diff = eligible_ys[i] - mean_y;
                ss_x += x_diff * x_diff;
                cov_xy += x_diff * y_diff;
            }

            if ss_x.is_finite() && cov_xy.is_finite() && ss_x > 0.0 {
                let slope = cov_xy / ss_x;
                let intercept = mean_y - slope * mean_x;

                if slope.is_finite() && intercept.is_finite() {
                    regression_coefficients = Some(RegressionCoefficients { intercept, slope });
                }
            }
        }
    }

    let mut fitted_dispersion = Vec::with_capacity(n_guides);
    let mut final_dispersion = Vec::with_capacity(n_guides);
    let mut size = Vec::with_capacity(n_guides);

    for g in 0..n_guides {
        let m = mean[g];

        if m == 0.0 {
            fitted_dispersion.push(POISSON_DISPERSION);
            final_dispersion.push(POISSON_DISPERSION);
            size.push(POISSON_SIZE);
            continue;
        }

        let mut fitted = median_dispersion;
        if let Some(coeff) = &regression_coefficients {
            let pred = (coeff.intercept + coeff.slope * m.ln()).exp();
            if pred.is_finite() && pred > 0.0 {
                fitted = pred;
            }
        }
        fitted_dispersion.push(fitted);

        let mut final_disp = fitted;
        if let Some(raw) = raw_dispersion[g] {
            if raw > fitted {
                final_disp = raw;
            }
        }

        if !final_disp.is_finite() || final_disp <= 0.0 {
            return Err(MageckError::Math(
                "Invalid final dispersion calculated".into(),
            ));
        }

        final_dispersion.push(final_disp);
        let s = 1.0 / final_disp;
        if !s.is_finite() || s <= 0.0 {
            return Err(MageckError::Math(
                "Invalid size parameter calculated".into(),
            ));
        }
        size.push(s);
    }

    Ok(NbModel {
        mean,
        variance,
        raw_dispersion,
        fitted_dispersion,
        final_dispersion,
        size,
        eligible_guides_count,
        regression_coefficients,
        median_dispersion,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn varying_means() -> CountMatrix {
        CountMatrix::from_vec(
            vec![
                100.0, 20.0, 300.0, 8.0, 5.0, 40.0, 100.0, 8.0, 3.0, 10.0, 60.0, 8.0, 4.0, 60.0,
                40.0, 8.0,
            ],
            4,
            4,
        )
        .unwrap()
    }

    #[test]
    fn p4f24_f_sh3_wls_matches_centered_python_and_differs_from_ols() {
        let counts = varying_means();
        let pooled = NbModel::fit_pooled_dispersion(&counts, 0, &[0, 1, 2, 3]).unwrap();
        let ordinary = fit_meanvar(&counts, &[0, 1, 2, 3]).unwrap();
        let wls = pooled.regression_coefficients().unwrap();
        let ols = ordinary.regression_coefficients().unwrap();
        // F_SH3: independently centered Python WLS on the identical fixture.
        // Absolute tolerance checks arithmetic, not agreement with another model.
        assert!((wls.slope - -0.15593805487328874).abs() <= 1e-9);
        assert!((wls.intercept - 0.6270541639895113).abs() <= 1e-9);
        assert!((wls.slope - ols.slope).abs() > 0.1);
        assert!((wls.intercept - ols.intercept).abs() > 0.1);
        assert_eq!(pooled.eligible_guides_count(), 3);
        assert_eq!(pooled.raw_dispersion()[3], None);
        assert_eq!(pooled.variance(), ordinary.variance());
        assert_eq!(pooled.raw_dispersion(), ordinary.raw_dispersion());
        assert_eq!(pooled.median_dispersion(), ordinary.median_dispersion());
    }

    #[test]
    fn p4f24_shrinkage_discards_large_raw_dispersion() {
        let counts = varying_means();
        let model = NbModel::fit_pooled_dispersion(&counts, 0, &[0, 1, 2, 3]).unwrap();
        assert!(model.raw_dispersion()[0].unwrap() > 3.0 * model.fitted_dispersion()[0]);
        assert_eq!(model.final_dispersion(), model.fitted_dispersion());
        for (&dispersion, &size) in model.final_dispersion().iter().zip(model.size()) {
            assert_eq!(size, 1.0 / dispersion);
        }
    }

    #[test]
    fn p4f24_fitted_trend_uses_control_mean_including_ineligible_guide() {
        let counts = varying_means();
        let pooled = fit_meanvar(&counts, &[0, 1, 2, 3]).unwrap();
        let model = NbModel::fit_pooled_dispersion(&counts, 0, &[0, 1, 2, 3]).unwrap();
        let coeff = model.regression_coefficients().unwrap();
        assert_eq!(model.mean(), counts.sample_slice(0));
        for g in 0..counts.n_guides() {
            let at_control = (coeff.intercept + coeff.slope * counts.entry(0, g).ln()).exp();
            assert_eq!(model.fitted_dispersion()[g], at_control);
        }
        let at_pooled = (coeff.intercept + coeff.slope * pooled.mean()[0].ln()).exp();
        assert!((model.fitted_dispersion()[0] - at_pooled).abs() > 0.1);
        assert!((model.final_dispersion()[0] - 0.9129564337871559).abs() < 1e-12);
    }

    #[test]
    fn p4f24_multi_control_preserves_conservative_maximum() {
        let counts = varying_means();
        let model = fit_meanvar(&counts, &[0, 1, 2, 3]).unwrap();
        assert!(model.raw_dispersion()[0].unwrap() > model.fitted_dispersion()[0]);
        for g in 0..counts.n_guides() {
            let expected = model.raw_dispersion()[g].map_or(model.fitted_dispersion()[g], |raw| {
                raw.max(model.fitted_dispersion()[g])
            });
            assert_eq!(model.final_dispersion()[g], expected);
            assert_eq!(model.size()[g], 1.0 / expected);
        }
    }

    #[test]
    fn p4f24_equal_pooled_means_fall_back_to_median() {
        let counts = CountMatrix::from_vec(
            vec![100.0; 4]
                .into_iter()
                .chain([5.0; 4])
                .chain([3.0; 4])
                .chain([4.0; 4])
                .collect(),
            4,
            4,
        )
        .unwrap();
        let model = NbModel::fit_pooled_dispersion(&counts, 0, &[0, 1, 2, 3]).unwrap();
        assert!(model.regression_coefficients().is_none());
        assert_eq!(model.final_dispersion(), &[model.median_dispersion(); 4]);
    }

    #[test]
    fn p4f24_fewer_than_three_eligible_guides_shrink_to_median() {
        let counts = CountMatrix::from_vec(vec![100.0, 20.0, 0.0, 10.0], 2, 2).unwrap();
        let model = NbModel::fit_pooled_dispersion(&counts, 0, &[0, 1]).unwrap();
        assert_eq!(model.eligible_guides_count(), 2);
        assert!(model.regression_coefficients().is_none());
        assert!(model.raw_dispersion()[0].unwrap() > model.median_dispersion());
        assert_eq!(model.final_dispersion(), &[model.median_dispersion(); 2]);
    }

    #[test]
    fn p4f24_zero_control_keeps_poisson_sentinel_with_wls() {
        let counts = CountMatrix::from_vec(
            vec![
                100.0, 20.0, 300.0, 0.0, 5.0, 40.0, 100.0, 8.0, 3.0, 10.0, 60.0, 8.0, 4.0, 60.0,
                40.0, 8.0,
            ],
            4,
            4,
        )
        .unwrap();
        let model = NbModel::fit_pooled_dispersion(&counts, 0, &[0, 1, 2, 3]).unwrap();
        assert!(model.regression_coefficients().is_some());
        assert_eq!(model.mean()[3], 0.0);
        assert_eq!(model.fitted_dispersion()[3], POISSON_DISPERSION);
        assert_eq!(model.final_dispersion()[3], POISSON_DISPERSION);
        assert_eq!(model.size()[3], POISSON_SIZE);
    }

    #[test]
    fn p4f24_no_eligible_guides_use_poisson_fallback() {
        let counts = CountMatrix::from_vec(vec![10.0, 0.0, 10.0, 0.0], 2, 2).unwrap();
        let model = NbModel::fit_pooled_dispersion(&counts, 0, &[0, 1]).unwrap();
        assert_eq!(model.eligible_guides_count(), 0);
        assert!(model.regression_coefficients().is_none());
        assert_eq!(model.final_dispersion(), &[POISSON_DISPERSION; 2]);
        assert_eq!(model.size(), &[POISSON_SIZE; 2]);
    }
}
