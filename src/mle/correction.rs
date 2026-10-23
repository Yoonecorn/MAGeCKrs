//! Deterministic collection-level multiple-testing correction for MLE results.

use super::gene_result::GeneResult;
use super::serial::SerialMleOutput;
use crate::error::MageckError;
use crate::stats::bh_correction;

const FDR_EPSILON: f64 = 1e-12;

/// Adjusts two-sided Wald p-values in-place across all genes per condition using BH FDR.
pub(crate) fn adjust_mle_fdr(results: &mut [GeneResult]) -> Result<(), MageckError> {
    if results.is_empty() {
        return Ok(());
    }

    let n_conditions = results[0].condition_order.len();
    if n_conditions == 0 {
        return Err(MageckError::Math(
            "serial MLE FDR requires at least one condition".into(),
        ));
    }

    let reference_order = &results[0].condition_order;
    for gene in results.iter() {
        if !gene.status.converged() {
            return Err(MageckError::Math(
                "BH requires converged eligible targeting genes".into(),
            ));
        }
        if gene.condition_order.len() != n_conditions
            || gene.two_sided_p_value.len() != n_conditions
            || &gene.condition_order != reference_order
        {
            return Err(MageckError::Math(
                "serial MLE FDR condition alignment mismatch across genes".into(),
            ));
        }
        for &p_value in &gene.two_sided_p_value {
            if !p_value.is_finite() || !(0.0..=1.0).contains(&p_value) {
                return Err(MageckError::Math(format!(
                    "serial MLE FDR requires finite p-values in [0, 1], got {p_value} for gene {}",
                    gene.gene_id
                )));
            }
        }
    }

    let n_genes = results.len();
    let mut condition_fdrs = Vec::with_capacity(n_conditions);
    for condition in 0..n_conditions {
        let mut order: Vec<_> = (0..n_genes).collect();
        order.sort_by(|&a, &b| {
            results[a].two_sided_p_value[condition]
                .total_cmp(&results[b].two_sided_p_value[condition])
                .then_with(|| {
                    results[a]
                        .gene_id
                        .as_bytes()
                        .cmp(results[b].gene_id.as_bytes())
                })
        });
        let p_values: Vec<_> = order
            .iter()
            .map(|&i| results[i].two_sided_p_value[condition])
            .collect();
        let adjusted = bh_correction(&p_values);
        let mut fdr_values = vec![0.0; n_genes];
        for (&i, q) in order.iter().zip(adjusted) {
            fdr_values[i] = q;
        }
        if fdr_values.len() != n_genes {
            return Err(MageckError::Math(
                "BH correction output length mismatch".into(),
            ));
        }
        condition_fdrs.push(fdr_values);
    }

    for (gene_index, gene) in results.iter().enumerate() {
        for (condition, fdr_by_gene) in condition_fdrs.iter().enumerate() {
            let fdr = fdr_by_gene[gene_index];
            let p_value = gene.two_sided_p_value[condition];
            if !fdr.is_finite() || !(0.0..=1.0).contains(&fdr) || fdr + FDR_EPSILON < p_value {
                return Err(MageckError::Math(
                    "computed FDR violates probability post-conditions".into(),
                ));
            }
        }
    }

    for (gene_index, gene) in results.iter_mut().enumerate() {
        let mut fdr_values = Vec::with_capacity(n_conditions);
        for fdr_by_gene in &condition_fdrs {
            fdr_values.push(fdr_by_gene[gene_index]);
        }
        gene.condition_fdr = fdr_values;
    }

    Ok(())
}

/// Convenience wrapper for adjusting serial MLE output in-place.
pub(crate) fn adjust_mle_output_fdr(output: &mut SerialMleOutput) -> Result<(), MageckError> {
    adjust_mle_fdr(&mut output.results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mle::irls::IrlsStatus;

    fn gene(id: &str, conditions: &[&str], p_values: &[f64]) -> GeneResult {
        GeneResult {
            gene_id: id.to_owned(),
            condition_order: conditions
                .iter()
                .map(|condition| (*condition).to_owned())
                .collect(),
            condition_beta: Vec::new(),
            condition_covariance: Vec::new(),
            standard_error: Vec::new(),
            z_score: Vec::new(),
            two_sided_p_value: p_values.to_vec(),
            condition_fdr: vec![0.99],
            status: IrlsStatus::Converged,
            attempts: 1,
            efficiencies: None,
            means: Vec::new(),
            residuals: Vec::new(),
        }
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() <= 1e-12, "{actual} != {expected}");
    }

    #[test]
    fn fixture_one_matches_condition_wise_direct_oracle() {
        let mut results = vec![
            gene(
                "GENE_ALPHA",
                &["beta_1", "beta_2"],
                &[0.11158253268494193, 0.35880285802687006],
            ),
            gene(
                "GENE_BETA",
                &["beta_1", "beta_2"],
                &[0.6163517079559528, 0.14808795825040513],
            ),
            gene(
                "GENE_GAMMA",
                &["beta_1", "beta_2"],
                &[0.05603914247267051, 0.023114238157811065],
            ),
        ];

        adjust_mle_fdr(&mut results).unwrap();

        let expected = [
            [0.1673737990274129, 0.35880285802687006],
            [0.6163517079559528, 0.2221319373756077],
            [0.1673737990274129, 0.0693427144734332],
        ];
        for (result, expected_fdr) in results.iter().zip(expected) {
            assert_eq!(result.condition_fdr.len(), 2);
            for (actual, expected_value) in result.condition_fdr.iter().zip(expected_fdr) {
                assert_close(*actual, expected_value);
            }
        }
        assert_eq!(
            results
                .iter()
                .map(|result| &result.gene_id)
                .collect::<Vec<_>>(),
            ["GENE_ALPHA", "GENE_BETA", "GENE_GAMMA"]
        );
    }

    #[test]
    fn fixture_two_matches_condition_wise_direct_oracle() {
        let mut results = vec![
            gene(
                "GENE_ALPHA",
                &["beta_1", "beta_2"],
                &[0.14179143907100336, 0.16846300350455112],
            ),
            gene(
                "GENE_BETA",
                &["beta_1", "beta_2"],
                &[0.07192937001547361, 0.08977077095386152],
            ),
        ];

        adjust_mle_fdr(&mut results).unwrap();

        assert_close(results[0].condition_fdr[0], 0.14179143907100336);
        assert_close(results[1].condition_fdr[0], 0.14179143907100336);
        assert_close(results[0].condition_fdr[1], 0.16846300350455112);
        assert_close(results[1].condition_fdr[1], 0.16846300350455112);
    }

    #[test]
    fn retained_nonconverged_fixture_is_corrected() {
        let mut results = vec![
            gene("GENE_NORMAL", &["beta_1"], &[0.08715633426895365]),
            gene("GENE_EXHAUSTED", &["beta_1"], &[0.9998912604311312]),
        ];
        results[1].status = IrlsStatus::IterationLimitReached;
        assert!(adjust_mle_fdr(&mut results).is_err());
        assert_eq!(results[0].condition_fdr, vec![0.99]);
    }

    #[test]
    fn empty_and_single_gene_collections_have_defined_semantics() {
        let mut empty: Vec<GeneResult> = Vec::new();
        assert!(adjust_mle_fdr(&mut empty).is_ok());

        let mut single = vec![gene("GENE", &["beta_1"], &[0.042])];
        adjust_mle_fdr(&mut single).unwrap();
        assert_eq!(single[0].condition_fdr, vec![0.042]);
    }

    #[test]
    fn exact_ties_receive_identical_fdr_and_replace_prior_values() {
        let mut results = vec![
            gene("GENE_A", &["beta_1"], &[0.2]),
            gene("GENE_B", &["beta_1"], &[0.2]),
            gene("GENE_C", &["beta_1"], &[0.8]),
        ];
        results[0].condition_fdr = vec![0.1, 0.2];

        adjust_mle_fdr(&mut results).unwrap();

        assert_eq!(results[0].condition_fdr.len(), 1);
        assert_eq!(results[0].condition_fdr, results[1].condition_fdr);
        assert_close(results[0].condition_fdr[0], 0.3);
    }

    #[test]
    fn invalid_alignment_and_probability_domains_are_rejected() {
        let mut label_mismatch = vec![
            gene("GENE_A", &["beta_1", "beta_2"], &[0.1, 0.2]),
            gene("GENE_B", &["beta_1", "other"], &[0.3, 0.4]),
        ];
        assert!(matches!(
            adjust_mle_fdr(&mut label_mismatch),
            Err(MageckError::Math(_))
        ));

        let mut nan = vec![gene("GENE", &["beta_1"], &[f64::NAN])];
        assert!(matches!(
            adjust_mle_fdr(&mut nan),
            Err(MageckError::Math(_))
        ));

        let mut out_of_range = vec![gene("GENE", &["beta_1"], &[1.05])];
        assert!(matches!(
            adjust_mle_fdr(&mut out_of_range),
            Err(MageckError::Math(_))
        ));
    }
}
