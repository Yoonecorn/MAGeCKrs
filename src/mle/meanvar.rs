//! Bounded translation of MAGeCK2's classical MeanVarModel LM route.
#![allow(dead_code)] // Registered privately before the later MLE orchestration consumes it.

use crate::error::MageckError;

const RESIDUAL_STABILIZER: f64 = 0.01;
const VARIANCE_FLOOR: f64 = 0.01;
const ALPHA_FLOOR: f64 = 0.01;
const FILTER_THRESHOLD: f64 = -1.0;
const FILTERED_MINIMUM: usize = 20;
const ALPHA_REFERENCE_MEAN: f64 = 10.0;

/// Finite, order-preserving observations accepted by the classical LM route.
#[derive(Debug, Clone)]
pub(super) struct MeanVarObservations {
    pub(super) means: Vec<f64>,
    pub(super) residuals: Vec<f64>,
}

impl MeanVarObservations {
    pub(super) fn new(means: Vec<f64>, residuals: Vec<f64>) -> Result<Self, MageckError> {
        if means.is_empty() {
            return Err(MageckError::Math(
                "MeanVarModel observations must not be empty".to_string(),
            ));
        }
        if means.len() != residuals.len() {
            return Err(MageckError::Math(
                "MeanVarModel means and residuals must have equal length".to_string(),
            ));
        }
        if means.len() < 2 {
            return Err(MageckError::Math(
                "MeanVarModel requires at least two observations".to_string(),
            ));
        }

        for (index, (&mean, &residual)) in means.iter().zip(&residuals).enumerate() {
            checked_mean(mean, index)?;
            if !residual.is_finite() {
                return Err(MageckError::Math(format!(
                    "MeanVarModel residual at observation {index} must be finite"
                )));
            }
            let predictor = mean.ln();
            let response = (residual * residual + RESIDUAL_STABILIZER).ln();
            if !predictor.is_finite() || !response.is_finite() {
                return Err(MageckError::Math(format!(
                    "MeanVarModel transformed observation {index} must be finite"
                )));
            }
        }

        Ok(Self { means, residuals })
    }
}

/// Fitted finite classical mean-variance model. This is private until the MLE
/// orchestration supplies validated, flattened EM observations.
#[derive(Debug, Clone, Copy)]
pub(super) struct MeanVarModel {
    pub(super) slope: f64,
    pub(super) intercept: f64,
}

impl MeanVarModel {
    pub(super) fn fit_from_observations(
        observations: MeanVarObservations,
    ) -> Result<Self, MageckError> {
        let mut all = Vec::with_capacity(observations.means.len());
        let mut filtered = Vec::with_capacity(observations.means.len());

        for (&mean, &residual) in observations.means.iter().zip(&observations.residuals) {
            let point = (mean.ln(), (residual * residual + RESIDUAL_STABILIZER).ln());
            if !point.0.is_finite() || !point.1.is_finite() {
                return Err(MageckError::Math(
                    "MeanVarModel transformed observation must be finite".to_string(),
                ));
            }
            all.push(point);
            if point.1 > FILTER_THRESHOLD {
                filtered.push(point);
            }
        }

        let points = if filtered.len() > FILTERED_MINIMUM {
            &filtered
        } else {
            &all
        };
        let (slope, intercept) = ordinary_least_squares(points)?;
        Ok(Self { slope, intercept })
    }

    fn variance(&self, mean: f64) -> Result<f64, MageckError> {
        checked_mean(mean, 0)?;
        let log_variance = mean.ln() * self.slope + self.intercept;
        if !log_variance.is_finite() {
            return Err(MageckError::Math(
                "MeanVarModel fitted log variance must be finite".to_string(),
            ));
        }
        let raw = log_variance.exp() - RESIDUAL_STABILIZER;
        if !raw.is_finite() {
            return Err(MageckError::Math(
                "MeanVarModel fitted variance must be finite".to_string(),
            ));
        }
        let variance = raw.max(VARIANCE_FLOOR);
        if !variance.is_finite() {
            return Err(MageckError::Math(
                "MeanVarModel bounded variance must be finite".to_string(),
            ));
        }
        Ok(variance)
    }

    fn alpha_floor(&self) -> Result<f64, MageckError> {
        let variance_at_ten = self.variance(ALPHA_REFERENCE_MEAN)?;
        let floor = ((variance_at_ten - ALPHA_REFERENCE_MEAN)
            / (ALPHA_REFERENCE_MEAN * ALPHA_REFERENCE_MEAN))
            .max(ALPHA_FLOOR);
        if !floor.is_finite() {
            return Err(MageckError::Math(
                "MeanVarModel alpha floor must be finite".to_string(),
            ));
        }
        Ok(floor)
    }

    pub(super) fn alpha(&self, mean: f64) -> Result<f64, MageckError> {
        let variance = self.variance(mean)?;
        let floor = self.alpha_floor()?;
        let raw = (variance - mean) / (mean * mean);
        if !raw.is_finite() {
            return Err(MageckError::Math(
                "MeanVarModel fitted alpha must be finite".to_string(),
            ));
        }
        let alpha = raw.max(floor);
        if !alpha.is_finite() {
            return Err(MageckError::Math(
                "MeanVarModel bounded alpha must be finite".to_string(),
            ));
        }
        Ok(alpha)
    }
}

fn checked_mean(mean: f64, index: usize) -> Result<(), MageckError> {
    if !mean.is_finite() || mean <= 0.0 {
        return Err(MageckError::Math(format!(
            "MeanVarModel mean at observation {index} must be finite and strictly greater than zero"
        )));
    }
    Ok(())
}

fn ordinary_least_squares(points: &[(f64, f64)]) -> Result<(f64, f64), MageckError> {
    let count = points.len() as f64;
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    for &(x, y) in points {
        sum_x += x;
        sum_y += y;
    }
    let mean_x = sum_x / count;
    let mean_y = sum_y / count;
    if !mean_x.is_finite() || !mean_y.is_finite() {
        return Err(MageckError::Math(
            "MeanVarModel OLS means must be finite".to_string(),
        ));
    }

    let mut predictor_variance = 0.0;
    let mut covariance = 0.0;
    for &(x, y) in points {
        predictor_variance += (x - mean_x) * (x - mean_x);
        covariance += (x - mean_x) * (y - mean_y);
    }
    predictor_variance /= count;
    covariance /= count;
    if !predictor_variance.is_finite() || predictor_variance == 0.0 {
        return Err(MageckError::Math(
            "MeanVarModel predictor variance must be finite and non-zero".to_string(),
        ));
    }
    let slope = covariance / predictor_variance;
    let intercept = mean_y - slope * mean_x;
    if !slope.is_finite() || !intercept.is_finite() {
        return Err(MageckError::Math(
            "MeanVarModel OLS coefficients must be finite".to_string(),
        ));
    }
    Ok((slope, intercept))
}

#[cfg(test)]
mod tests {
    use super::*;

    const ABS_TOLERANCE: f64 = 2e-12;
    const REL_TOLERANCE: f64 = 2e-12;

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = ABS_TOLERANCE + REL_TOLERANCE * expected.abs();
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:.17e}"
        );
    }

    struct OracleFixture {
        means: Vec<f64>,
        residuals: Vec<f64>,
        slope: f64,
        intercept: f64,
        queries: Vec<f64>,
        variances: Vec<f64>,
        alpha_floor: f64,
        alphas: Vec<f64>,
    }

    fn assert_fixture(fixture: OracleFixture) {
        let model = MeanVarModel::fit_from_observations(
            MeanVarObservations::new(fixture.means, fixture.residuals).unwrap(),
        )
        .unwrap();
        assert_close(model.slope, fixture.slope);
        assert_close(model.intercept, fixture.intercept);
        assert_close(model.alpha_floor().unwrap(), fixture.alpha_floor);
        for ((&query, &variance), &alpha) in fixture
            .queries
            .iter()
            .zip(&fixture.variances)
            .zip(&fixture.alphas)
        {
            assert_close(model.variance(query).unwrap(), variance);
            assert_close(model.alpha(query).unwrap(), alpha);
        }
    }

    #[test]
    fn canonical_minimal_oracle() {
        assert_fixture(OracleFixture {
            means: vec![10.0, 20.0, 25.0, 40.0],
            residuals: vec![1.0, 2.0, 3.0, 4.0],
            slope: 2.054_113_127_179_705,
            intercept: -4.675_590_970_423_031,
            queries: vec![1.0, 10.0, 20.0, 100.0],
            variances: vec![
                0.01,
                1.045_673_862_807_654_8,
                4.374_089_735_060_906,
                119.565_690_182_357_1,
            ],
            alpha_floor: 0.01,
            alphas: vec![0.01, 0.01, 0.01, 0.01],
        });
    }

    #[test]
    fn realistic_multi_gene_oracle_preserves_order() {
        let means = vec![8.0, 16.0, 11.0, 22.0, 15.0, 30.0, 6.0, 12.0, 24.0];
        let observations = MeanVarObservations::new(
            means.clone(),
            vec![0.2, 0.8, 0.5, 1.4, 0.9, 2.1, 0.1, 0.7, 1.8],
        )
        .unwrap();
        assert_eq!(observations.means, means);
        assert_fixture(OracleFixture {
            means: observations.means,
            residuals: observations.residuals,
            slope: 3.410_360_939_988_364,
            intercept: -9.750_982_233_890_11,
            queries: vec![2.0, 6.0, 12.0, 30.0, 60.0],
            variances: vec![
                0.01,
                0.016_240_907_623_861_18,
                0.268_997_165_716_172_3,
                6.339_197_746_067_745,
                67.495_598_553_064_04,
            ],
            alpha_floor: 0.01,
            alphas: vec![0.01, 0.01, 0.01, 0.01, 0.01],
        });
    }

    #[test]
    fn filtered_gt_twenty_oracle_and_alpha_floor_boundary() {
        let model = MeanVarModel::fit_from_observations(
            MeanVarObservations::new((1..=21).map(f64::from).collect(), vec![1.0; 21]).unwrap(),
        )
        .unwrap();
        assert_close(model.slope, 0.0);
        assert_close(model.intercept, 0.009_950_330_853_168_083);
        assert_close(model.variance(0.5).unwrap(), 1.0);
        assert_close(model.alpha_floor().unwrap(), 0.01);
        assert_close(model.alpha(0.5).unwrap(), 2.0);
        assert_eq!(model.alpha(1.0).unwrap(), model.alpha_floor().unwrap());
    }

    #[test]
    fn exactly_twenty_filtered_points_uses_all_observations() {
        let means: Vec<f64> = (1..=21).map(f64::from).collect();
        let mut residuals = vec![1.0; 20];
        residuals.push(0.0);
        let model = MeanVarModel::fit_from_observations(
            MeanVarObservations::new(means, residuals).unwrap(),
        )
        .unwrap();
        assert!(
            model.slope.abs() > 1e-3,
            "all-observations route must include y <= -1"
        );
    }

    #[test]
    fn truncation_and_order_oracle_observations() {
        assert_fixture(OracleFixture {
            means: vec![3.0, 9.0, 4.0, 16.0],
            residuals: vec![0.5, 1.5, 0.25, 2.5],
            slope: 2.411_986_450_257_406,
            intercept: -4.825_595_218_755_365,
            queries: vec![3.0, 4.0, 9.0, 16.0],
            variances: vec![
                0.103_522_154_826_298,
                0.217_211_948_414_373_64,
                1.596_536_613_160_072_5,
                6.425_639_480_137_809,
            ],
            alpha_floor: 0.01,
            alphas: vec![0.01, 0.01, 0.01, 0.01],
        });
    }

    #[test]
    fn repeated_calls_are_deterministic_and_variance_floor_is_exact() {
        let model = MeanVarModel::fit_from_observations(
            MeanVarObservations::new(vec![10.0, 20.0, 25.0, 40.0], vec![1.0, 2.0, 3.0, 4.0])
                .unwrap(),
        )
        .unwrap();
        assert_eq!(model.variance(1.0).unwrap(), 0.01);
        assert_eq!(model.variance(20.0).unwrap(), model.variance(20.0).unwrap());
        assert_eq!(model.alpha(20.0).unwrap(), model.alpha(20.0).unwrap());
    }

    #[test]
    fn rejects_malformed_lengths_and_empty_or_single_training_inputs() {
        assert!(MeanVarObservations::new(vec![], vec![]).is_err());
        assert!(MeanVarObservations::new(vec![1.0], vec![1.0, 2.0]).is_err());
        assert!(MeanVarObservations::new(vec![1.0], vec![1.0]).is_err());
    }

    #[test]
    fn rejects_non_positive_non_finite_and_degenerate_inputs() {
        assert!(MeanVarObservations::new(vec![0.0, 1.0], vec![1.0, 1.0]).is_err());
        assert!(MeanVarObservations::new(vec![-1.0, 1.0], vec![1.0, 1.0]).is_err());
        assert!(MeanVarObservations::new(vec![f64::NAN, 1.0], vec![1.0, 1.0]).is_err());
        assert!(MeanVarObservations::new(vec![1.0, 2.0], vec![f64::INFINITY, 1.0]).is_err());
        let observations = MeanVarObservations::new(vec![2.0, 2.0], vec![1.0, 2.0]).unwrap();
        assert!(MeanVarModel::fit_from_observations(observations).is_err());
    }

    #[test]
    fn rejects_invalid_evaluation_means() {
        let model = MeanVarModel::fit_from_observations(
            MeanVarObservations::new(vec![1.0, 2.0], vec![1.0, 2.0]).unwrap(),
        )
        .unwrap();
        assert!(model.variance(0.0).is_err());
        assert!(model.alpha(-1.0).is_err());
        assert!(model.variance(f64::INFINITY).is_err());
    }
}
