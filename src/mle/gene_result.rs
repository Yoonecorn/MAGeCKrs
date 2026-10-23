//! Private, owned one-gene fitted-result construction.

use super::{
    fit::FitResult,
    wald::{project_wald, WaldInput},
};
use crate::error::MageckError;

#[derive(Clone, Debug)]
pub(super) struct GeneResult {
    pub(super) gene_id: String,
    pub(super) condition_order: Vec<String>,
    pub(super) condition_beta: Vec<f64>,
    pub(super) condition_covariance: Vec<f64>,
    pub(super) standard_error: Vec<f64>,
    pub(super) z_score: Vec<f64>,
    pub(super) two_sided_p_value: Vec<f64>,
    pub(super) condition_fdr: Vec<f64>,
    pub(super) status: super::irls::IrlsStatus,
    pub(super) attempts: usize,
    pub(super) efficiencies: Option<Vec<f64>>,
    pub(super) means: Vec<f64>,
    pub(super) residuals: Vec<f64>,
}

pub(super) fn construct_gene_result(
    gene_id: &str,
    fit: &FitResult,
    condition_order: &[String],
) -> Result<GeneResult, MageckError> {
    let n_conditions = condition_order.len();
    let n_guides =
        fit.beta.len().checked_sub(n_conditions).ok_or_else(|| {
            MageckError::Math("gene-result beta shorter than condition order".into())
        })?;
    let wald = project_wald(WaldInput {
        full_beta: &fit.beta,
        n_guides,
        n_conditions,
        condition_covariance: &fit.condition_se,
        status: Some(fit.status),
        attempts: fit.attempts,
        means: &fit.means,
        residuals: &fit.residuals,
        efficiencies: fit.efficiencies.as_deref(),
        condition_order,
    })?;
    Ok(GeneResult {
        gene_id: gene_id.to_owned(),
        condition_order: wald.condition_order,
        condition_beta: wald.condition_beta,
        condition_covariance: wald.condition_covariance,
        standard_error: wald.standard_error,
        z_score: wald.z_score,
        two_sided_p_value: wald.two_sided_p_value,
        condition_fdr: Vec::new(),
        status: wald.status,
        attempts: wald.attempts,
        efficiencies: wald.efficiencies,
        means: wald.means,
        residuals: wald.residuals,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mle::{
        fit::FitResult,
        irls::IrlsStatus,
        wald::{project_wald, WaldInput},
    };

    fn names(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_owned()).collect()
    }

    fn fit(status: IrlsStatus, attempts: usize, efficiencies: Option<Vec<f64>>) -> FitResult {
        FitResult::gene_result_fixture(
            status,
            attempts,
            vec![1.5, 1.6, -1.2, 0.05, 2.0],
            efficiencies,
            vec![0.04, 0.01, -0.005, 0.01, 0.25, 0.03, -0.005, 0.03, 0.01],
            vec![5.0, 5.4, 4.8, 5.3, 6.1, 6.4, 3.2, 3.5],
            vec![0.2, -0.1, 0.3, -0.2, 0.1, -0.1, 0.4, -0.3],
        )
    }

    #[test]
    fn converged_construction_preserves_complete_covariance() {
        let order = names(&["negative_first", "small_middle", "positive_last"]);
        let result =
            construct_gene_result("GENE_A", &fit(IrlsStatus::Converged, 1, None), &order).unwrap();
        assert_eq!(result.gene_id, "GENE_A");
        assert_eq!(result.condition_order, order);
        assert_eq!(result.condition_beta, vec![-1.2, 0.05, 2.0]);
        assert_eq!(result.condition_covariance.len(), 9);
        assert_eq!(result.condition_covariance[2], -0.005);
        for index in 0..3 {
            assert_eq!(
                result.standard_error[index],
                result.condition_covariance[index * 3 + index].sqrt()
            );
        }
        assert_eq!(result.efficiencies, None);
        assert_eq!(result.means, vec![5.0, 5.4, 4.8, 5.3, 6.1, 6.4, 3.2, 3.5]);
        assert_eq!(
            result.residuals,
            vec![0.2, -0.1, 0.3, -0.2, 0.1, -0.1, 0.4, -0.3]
        );
    }

    #[test]
    fn retry_exhaustion_and_enabled_efficiency_are_retained() {
        let order = names(&["first", "second", "third"]);
        assert!(construct_gene_result(
            "GENE_B",
            &fit(IrlsStatus::IterationLimitReached, 3, Some(vec![0.7, 0.8])),
            &order
        )
        .is_err());
    }

    #[test]
    fn projection_matches_direct_wald_exactly() {
        let order = names(&["first", "second", "third"]);
        let fit = fit(IrlsStatus::Converged, 1, None);
        let direct = project_wald(WaldInput {
            full_beta: &fit.beta,
            n_guides: 2,
            n_conditions: 3,
            condition_covariance: &fit.condition_se,
            status: Some(fit.status),
            attempts: fit.attempts,
            means: &fit.means,
            residuals: &fit.residuals,
            efficiencies: fit.efficiencies.as_deref(),
            condition_order: &order,
        })
        .unwrap();
        let result = construct_gene_result("GENE_C", &fit, &order).unwrap();
        assert_eq!(result.condition_beta, direct.condition_beta);
        assert_eq!(result.condition_covariance, direct.condition_covariance);
        assert_eq!(result.standard_error, direct.standard_error);
        assert_eq!(result.z_score, direct.z_score);
        assert_eq!(result.two_sided_p_value, direct.two_sided_p_value);
    }

    #[test]
    fn malformed_and_non_finite_propagated_inputs_are_rejected() {
        let order = names(&["first", "second", "third"]);
        let mut malformed = fit(IrlsStatus::Converged, 1, None);
        malformed.beta.truncate(2);
        assert!(construct_gene_result("GENE", &malformed, &order).is_err());

        let mut non_finite = fit(IrlsStatus::Converged, 1, None);
        non_finite.condition_se[0] = f64::NAN;
        assert!(construct_gene_result("GENE", &non_finite, &order).is_err());

        let mut non_positive = fit(IrlsStatus::Converged, 1, None);
        non_positive.condition_se[4] = 0.0;
        assert!(construct_gene_result("GENE", &non_positive, &order).is_err());
    }
}
