//! Deterministic, serial multi-gene classical MLE orchestration.

use super::{
    fit::{setup_and_fit, FitDispersion, FitInput},
    gene_result::{construct_gene_result, GeneResult},
    meanvar::{MeanVarModel, MeanVarObservations},
};
use crate::error::MageckError;

#[derive(Clone, Debug)]
pub(crate) struct SerialGeneInput {
    pub(crate) gene_id: String,
    pub(crate) sgrna_ids: Vec<String>,
    pub(crate) counts: Vec<f64>,
}

#[derive(Clone, Debug)]
pub(crate) struct SerialMleInput<'a> {
    pub(crate) genes: &'a [SerialGeneInput],
    pub(crate) n_samples: usize,
    pub(crate) design_matrix: &'a [f64],
    pub(crate) n_conditions: usize,
    pub(crate) condition_labels: &'a [String],
    pub(crate) size_factors: &'a [f64],
    pub(crate) genes_varmodeling: usize,
    pub(crate) alpha_val: f64,
    pub(crate) estimate_efficiency: bool,
    pub(crate) update_efficiency: bool,
}

#[derive(Clone, Debug)]
pub(crate) struct SerialMleOutput {
    pub(crate) records: Vec<super::diagnostics::Record>,
    pub(crate) failure: Option<String>,
    pub(crate) results: Vec<GeneResult>,
    pub(crate) meanvar_model: Option<MeanVarModel>,
}

fn math(message: impl Into<String>) -> MageckError {
    MageckError::Math(message.into())
}

fn checked_product(left: usize, right: usize, name: &str) -> Result<usize, MageckError> {
    left.checked_mul(right)
        .ok_or_else(|| math(format!("serial MLE {name} dimension overflow")))
}

fn validate(
    input: &SerialMleInput<'_>,
    training_genes: &[SerialGeneInput],
) -> Result<(), MageckError> {
    if input.n_samples < 2 || input.n_conditions == 0 {
        return Err(math(
            "serial MLE requires at least two samples and one condition",
        ));
    }
    if input.condition_labels.len() != input.n_conditions {
        return Err(math("serial MLE condition label count mismatch"));
    }
    if input.size_factors.len() != input.n_samples {
        return Err(math("serial MLE size factor count mismatch"));
    }
    if input
        .size_factors
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(math("serial MLE size factors must be finite and positive"));
    }
    let columns = input
        .n_conditions
        .checked_add(1)
        .ok_or_else(|| math("serial MLE design dimension overflow"))?;
    if input.design_matrix.len() != checked_product(input.n_samples, columns, "design")? {
        return Err(math("serial MLE design shape mismatch"));
    }
    if input.design_matrix.iter().any(|value| !value.is_finite()) {
        return Err(math("serial MLE design entries must be finite"));
    }
    if (0..input.n_samples).any(|row| input.design_matrix[row * columns] != 1.0)
        || input.design_matrix[..columns]
            != std::iter::once(1.0)
                .chain(std::iter::repeat(0.0).take(input.n_conditions))
                .collect::<Vec<_>>()
    {
        return Err(math("serial MLE design requires literal row zero baseline"));
    }
    if !input.alpha_val.is_finite() || input.alpha_val <= 0.0 {
        return Err(math("serial MLE alpha must be finite and positive"));
    }
    for gene in input
        .genes
        .iter()
        .chain(training_genes.iter().take(input.genes_varmodeling))
    {
        if gene.gene_id.is_empty() || gene.sgrna_ids.is_empty() {
            return Err(math(
                "serial MLE genes require non-empty identifiers and guides",
            ));
        }
        let expected = checked_product(input.n_samples, gene.sgrna_ids.len(), "count")?;
        if gene.counts.len() != expected {
            return Err(math("serial MLE gene count shape mismatch"));
        }
        if gene
            .counts
            .iter()
            .any(|count| !count.is_finite() || *count <= 0.0 || count.fract() != 0.0)
        {
            return Err(math("fit counts must be positive finite integers"));
        }
    }
    if input.genes_varmodeling > training_genes.len() {
        return Err(math("genes_varmodeling exceeds available gene count"));
    }
    Ok(())
}

fn fit_input(
    input: &SerialMleInput<'_>,
    gene: &SerialGeneInput,
    dispersion: FitDispersion,
    estimate_efficiency: bool,
    update_efficiency: bool,
) -> FitInput {
    FitInput {
        counts: gene.counts.clone(),
        n_samples: input.n_samples,
        n_guides: gene.sgrna_ids.len(),
        design: input.design_matrix.to_vec(),
        n_conditions: input.n_conditions,
        size_factors: input.size_factors.to_vec(),
        dispersion,
        alpha_val: input.alpha_val,
        estimate_efficiency,
        restart_beta: None,
        restart_efficiencies: None,
        update_efficiency,
    }
}

fn vector_dispersion(
    model: MeanVarModel,
    gene: &SerialGeneInput,
) -> Result<FitDispersion, MageckError> {
    let guides = gene.sgrna_ids.len();
    let mut counts = gene.counts.clone();
    counts.extend_from_slice(&gene.counts[guides..]);
    counts
        .into_iter()
        .map(|count| model.alpha(count))
        .collect::<Result<Vec<_>, _>>()
        .map(FitDispersion::Vector)
}

#[cfg(test)]
pub(crate) fn run_mle_serial(
    input: &SerialMleInput<'_>,
    training_genes: &[SerialGeneInput],
) -> Result<SerialMleOutput, MageckError> {
    run_mle_with_execution(
        input,
        training_genes,
        &crate::execution::Execution::serial(),
    )
}

pub(crate) fn run_mle_with_execution(
    input: &SerialMleInput<'_>,
    training_genes: &[SerialGeneInput],
    execution: &crate::execution::Execution,
) -> Result<SerialMleOutput, MageckError> {
    validate(input, training_genes)?;
    let mut records = Vec::new();
    let mut failure = None;
    let meanvar_model = if input.genes_varmodeling == 0 {
        None
    } else {
        let max_guides = training_genes
            .iter()
            .take(input.genes_varmodeling)
            .map(|g| g.sgrna_ids.len())
            .max()
            .unwrap_or(0);
        let outcomes = execution.map_stage(
            "mle.training",
            crate::execution::mle_worker_bytes(input.n_samples, max_guides, input.n_conditions),
            training_genes
                .iter()
                .take(input.genes_varmodeling)
                .collect(),
            |gene| {
                Ok((
                    gene,
                    setup_and_fit(&fit_input(
                        input,
                        gene,
                        FitDispersion::Scalar(input.alpha_val),
                        false,
                        false,
                    )),
                ))
            },
        )?;
        let mut means = Vec::new();
        let mut residuals = Vec::new();
        for (gene, outcome) in outcomes {
            match outcome {
                Ok(fit) => {
                    let error = if !fit.status.converged() {
                        Some(format!("training fit {}", fit.status.name()))
                    } else {
                        None
                    };
                    if failure.is_none() {
                        failure = error
                            .as_ref()
                            .map(|e| format!("gene '{}': {e}", gene.gene_id));
                    }
                    let physical =
                        checked_product(input.n_samples, gene.sgrna_ids.len(), "observation")?;
                    if error.is_none() {
                        means.extend_from_slice(&gene.counts);
                        residuals.extend_from_slice(&fit.residuals[..physical]);
                    }
                    records.push(super::diagnostics::Record::new(
                        gene,
                        "training",
                        Some(fit),
                        error,
                    ));
                }
                Err(e) => {
                    let error = e.to_string();
                    if failure.is_none() {
                        failure = Some(format!("gene '{}': {error}", gene.gene_id));
                    }
                    records.push(super::diagnostics::Record::new(
                        gene,
                        "training",
                        None,
                        Some(error),
                    ));
                }
            }
        }
        if failure.is_some() {
            None
        } else {
            match MeanVarModel::fit_from_observations(MeanVarObservations::new(means, residuals)?) {
                Ok(m) => Some(m),
                Err(e) => {
                    failure = Some(e.to_string());
                    None
                }
            }
        }
    };
    if let Some(error) = failure {
        for gene in input.genes {
            records.push(super::diagnostics::Record::new(
                gene,
                "fit",
                None,
                Some("NOT_RUN: training failed".into()),
            ));
        }
        return Ok(SerialMleOutput {
            results: Vec::new(),
            records,
            failure: Some(error),
            meanvar_model,
        });
    }
    let max_guides = input
        .genes
        .iter()
        .map(|g| g.sgrna_ids.len())
        .max()
        .unwrap_or(0);
    let outcomes = execution.map_stage(
        "mle.fit",
        crate::execution::mle_worker_bytes(input.n_samples, max_guides, input.n_conditions),
        input.genes.iter().collect(),
        |gene| {
            let outcome = (|| {
                let dispersion = match meanvar_model {
                    Some(model) => vector_dispersion(model, gene)?,
                    None => FitDispersion::Scalar(input.alpha_val),
                };
                setup_and_fit(&fit_input(
                    input,
                    gene,
                    dispersion,
                    input.estimate_efficiency,
                    input.update_efficiency,
                ))
            })();
            Ok((gene, outcome))
        },
    )?;
    let mut results = Vec::new();
    for (gene, outcome) in outcomes {
        match outcome {
            Ok(fit) => {
                let projected = construct_gene_result(&gene.gene_id, &fit, input.condition_labels);
                let error = match projected {
                    Ok(result) => {
                        results.push(result);
                        None
                    }
                    Err(e) => Some(e.to_string()),
                };
                if failure.is_none() {
                    failure = error
                        .as_ref()
                        .map(|e| format!("gene '{}': {e}", gene.gene_id));
                }
                records.push(super::diagnostics::Record::new(
                    gene,
                    "fit",
                    Some(fit),
                    error,
                ));
            }
            Err(e) => {
                let error = e.to_string();
                if failure.is_none() {
                    failure = Some(format!("gene '{}': {error}", gene.gene_id));
                }
                records.push(super::diagnostics::Record::new(
                    gene,
                    "fit",
                    None,
                    Some(error),
                ));
            }
        }
    }
    Ok(SerialMleOutput {
        results,
        meanvar_model,
        records,
        failure,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn labels() -> Vec<String> {
        vec!["beta_1".to_owned(), "beta_2".to_owned()]
    }

    fn assert_close(actual: &[f64], expected: &[f64], tolerance: f64) {
        assert_eq!(actual.len(), expected.len());
        for (&actual, &expected) in actual.iter().zip(expected) {
            assert!(
                (actual - expected).abs() < tolerance,
                "{actual} != {expected}"
            );
        }
    }

    fn genes() -> Vec<SerialGeneInput> {
        vec![
            SerialGeneInput {
                gene_id: "GENE_ALPHA".to_owned(),
                sgrna_ids: vec!["a1".to_owned(), "a2".to_owned()],
                counts: vec![100., 120., 150., 180., 80., 90.],
            },
            SerialGeneInput {
                gene_id: "GENE_BETA".to_owned(),
                sgrna_ids: vec!["b1".to_owned(), "b2".to_owned(), "b3".to_owned()],
                counts: vec![50., 60., 70., 40., 50., 60., 100., 120., 140.],
            },
            SerialGeneInput {
                gene_id: "GENE_GAMMA".to_owned(),
                sgrna_ids: vec![
                    "c1".to_owned(),
                    "c2".to_owned(),
                    "c3".to_owned(),
                    "c4".to_owned(),
                ],
                counts: vec![
                    200., 250., 300., 350., 180., 220., 260., 310., 250., 300., 350., 400.,
                ],
            },
        ]
    }

    fn input<'a>(
        genes: &'a [SerialGeneInput],
        labels: &'a [String],
        k: usize,
    ) -> SerialMleInput<'a> {
        SerialMleInput {
            genes,
            n_samples: 3,
            design_matrix: &[1., 0., 0., 1., 1., 0., 1., 0., 1.],
            n_conditions: 2,
            condition_labels: labels,
            size_factors: &[1., 1.1, 0.9],
            genes_varmodeling: k,
            alpha_val: 0.01,
            estimate_efficiency: true,
            update_efficiency: true,
        }
    }

    #[test]
    fn fixture_01_two_pass_standard() {
        let genes = genes();
        let labels = labels();
        let output = run_mle_serial(&input(&genes, &labels, 3), &genes).unwrap();
        let model = output.meanvar_model.unwrap();
        assert!((model.slope - -0.06422946392139549).abs() < 1e-10);
        assert!((model.intercept - 2.9127244718126177).abs() < 1e-10);
        assert!((model.alpha(10.).unwrap() - 0.05866351351442944).abs() < 1e-12);
        assert!(output.failure.is_none());
        assert_close(
            &output.results[0].condition_beta,
            &[0.2924446626792171, -0.16865964333242214],
            1e-10,
        );
        assert_close(
            &output.results[0].standard_error,
            &[0.18379938904226006, 0.18381274971812891],
            1e-09,
        );
        assert_close(
            &output.results[0].z_score,
            &[1.5911079150104073, -0.917562266986683],
            1e-09,
        );
        assert_close(
            &output.results[0].two_sided_p_value,
            &[0.11158529155971224, 0.35884808000688284],
            1e-12,
        );
        assert_close(
            output.results[0].efficiencies.as_deref().unwrap(),
            &[0.6599673192387402, 0.6825893990024281],
            1e-10,
        );
        assert_close(
            &output.results[1].condition_beta,
            &[-0.27767741210008945, 0.8013685910825371],
            1e-10,
        );
        assert_close(
            &output.results[1].standard_error,
            &[0.5542135439444013, 0.5540763233621174],
            1e-09,
        );
        assert_close(
            &output.results[1].z_score,
            &[-0.5010296394487719, 1.4463144467532165],
            1e-09,
        );
        assert_close(
            &output.results[1].two_sided_p_value,
            &[0.6163502634723724, 0.14808902553045053],
            1e-12,
        );
        assert_close(
            output.results[1].efficiencies.as_deref().unwrap(),
            &[0.9969250353682094, 0.9970273784062338, 0.9970965526494261],
            1e-10,
        );
        assert_close(
            &output.results[2].condition_beta,
            &[-0.22846339795726658, 0.27162141483445185],
            1e-10,
        );
        assert_close(
            &output.results[2].standard_error,
            &[0.11957238229586488, 0.11957724623598802],
            1e-09,
        );
        assert_close(
            &output.results[2].z_score,
            &[-1.9106702866550433, 2.2715142168301963],
            1e-09,
        );
        assert_close(
            &output.results[2].two_sided_p_value,
            &[0.05604696556845347, 0.02311586732664111],
            1e-12,
        );
        assert_close(
            output.results[2].efficiencies.as_deref().unwrap(),
            &[
                0.6873767491652414,
                0.6630141279386322,
                0.645007159299658,
                0.62162815603346,
            ],
            1e-10,
        );
    }

    #[test]
    fn disabled_efficiency_and_exact_k_boundaries() {
        let genes = genes();
        let labels = labels();
        let mut disabled = input(&genes, &labels, 2);
        disabled.estimate_efficiency = false;
        disabled.update_efficiency = false;
        let disabled = run_mle_serial(&disabled, &genes).unwrap();
        assert!(disabled.results.iter().all(|r| r.efficiencies.is_none()));
        assert_close(
            &disabled.results[0].condition_beta,
            &[0.23753001817494712, -0.22265768016018284],
            1e-10,
        );
        assert_close(
            &disabled.results[0].standard_error,
            &[0.1616777794341553, 0.1616787457112203],
            1e-09,
        );
        assert_close(
            &disabled.results[0].z_score,
            &[1.4691568563488548, -1.377161105380413],
            1e-09,
        );
        assert_close(
            &disabled.results[0].two_sided_p_value,
            &[0.14179024928574668, 0.16846244400868793],
            1e-12,
        );
        assert_close(
            &disabled.results[1].condition_beta,
            &[-0.555524805605611, 0.523596631907591],
            1e-10,
        );
        assert_close(
            &disabled.results[1].standard_error,
            &[0.30870063559287675, 0.3086151653166655],
            1e-09,
        );
        assert_close(
            &disabled.results[1].z_score,
            &[-1.79955834732473, 1.6966004615176191],
            1e-09,
        );
        assert_close(
            &disabled.results[1].two_sided_p_value,
            &[0.07193040304780123, 0.08977222251102361],
            1e-12,
        );
        let zero = run_mle_serial(&input(&genes, &labels, 0), &genes).unwrap();
        assert!(zero.meanvar_model.is_none());
        assert_eq!(zero.results.len(), 3);
        assert!(run_mle_serial(&input(&genes, &labels, 4), &genes).is_err());
        let empty: Vec<SerialGeneInput> = Vec::new();
        assert!(run_mle_serial(&input(&empty, &labels, 0), &empty)
            .unwrap()
            .results
            .is_empty());
        assert!(run_mle_serial(&input(&empty, &labels, 1), &empty).is_err());
    }

    #[test]
    fn training_includes_skipped_genes_before_window_boundary() {
        let all_genes = genes();
        let labels = labels();
        // The middle, three-guide gene is training-only. Its position must
        // consume the second window slot, excluding the final gene.
        let fitted_genes = vec![all_genes[0].clone(), all_genes[2].clone()];
        let expected = run_mle_serial(&input(&all_genes, &labels, 2), &all_genes).unwrap();
        let actual = run_mle_serial(&input(&fitted_genes, &labels, 2), &all_genes).unwrap();
        let without_skipped =
            run_mle_serial(&input(&fitted_genes, &labels, 2), &fitted_genes).unwrap();
        let expected_model = expected.meanvar_model.unwrap();
        let actual_model = actual.meanvar_model.unwrap();
        assert_eq!(actual_model.slope.to_bits(), expected_model.slope.to_bits());
        assert_eq!(
            actual_model.intercept.to_bits(),
            expected_model.intercept.to_bits()
        );
        assert_ne!(
            actual_model.slope.to_bits(),
            without_skipped.meanvar_model.unwrap().slope.to_bits()
        );
        assert_eq!(actual.results.len(), 2);
        assert_eq!(actual.results[0].gene_id, "GENE_ALPHA");
        assert_eq!(actual.results[1].gene_id, "GENE_GAMMA");
        assert_close(
            &actual.results[0].condition_beta,
            &expected.results[0].condition_beta,
            1e-12,
        );

        // Training may exceed the number of final fits, including no final fits.
        let empty = Vec::new();
        let mut training_only = input(&empty, &labels, 3);
        assert!(run_mle_serial(&training_only, &all_genes)
            .unwrap()
            .results
            .is_empty());
        training_only.genes_varmodeling = 2;
        let output = run_mle_serial(&training_only, &all_genes).unwrap();
        assert_eq!(
            output.meanvar_model.unwrap().slope.to_bits(),
            expected_model.slope.to_bits()
        );
    }

    #[test]
    fn strict_validation_errors() {
        let mut genes = genes();
        let labels = labels();
        genes[0].counts[0] = 0.;
        assert!(matches!(
            run_mle_serial(&input(&genes, &labels, 0), &genes),
            Err(MageckError::Math(_))
        ));
        genes[0].counts[0] = 1.5;
        assert!(run_mle_serial(&input(&genes, &labels, 0), &genes).is_err());
        let mut bad = input(&genes, &labels, 0);
        bad.size_factors = &[1., -1., 0.9];
        assert!(run_mle_serial(&bad, &genes).is_err());
        let mut bad_design = input(&genes, &labels, 0);
        bad_design.design_matrix = &[1., 1., 0., 1., 1., 0., 1., 0., 1.];
        assert!(run_mle_serial(&bad_design, &genes).is_err());
    }

    #[test]
    fn retained_nonconverged_scalar_pass_is_preserved() {
        let labels = vec!["beta_1".to_owned()];
        let genes = vec![
            SerialGeneInput {
                gene_id: "GENE_NORMAL".to_owned(),
                sgrna_ids: vec!["n1".to_owned(), "n2".to_owned()],
                counts: vec![100., 150., 120., 170., 110., 160., 130., 180.],
            },
            SerialGeneInput {
                gene_id: "GENE_EXHAUSTED".to_owned(),
                sgrna_ids: vec!["e1".to_owned(), "e2".to_owned()],
                counts: vec![1., 100000., 100000., 1., 1., 1000000., 1000000., 1.],
            },
        ];
        let input = SerialMleInput {
            genes: &genes,
            n_samples: 4,
            design_matrix: &[1., 0., 1., 1., 1., 0., 1., 1.],
            n_conditions: 1,
            condition_labels: &labels,
            size_factors: &[1.; 4],
            genes_varmodeling: 0,
            alpha_val: 0.01,
            estimate_efficiency: false,
            update_efficiency: false,
        };
        let output = run_mle_serial(&input, &genes).unwrap();
        assert_eq!(output.results.len(), 1);
        assert!(output.failure.is_some());
        assert_eq!(output.records[1].fit.as_ref().unwrap().attempts, 3);
        assert!(!output.records[1].fit.as_ref().unwrap().status.converged());
    }
    #[test]
    fn threaded_fits_and_training_assembly_are_bit_exact() {
        let genes = genes();
        let labels = labels();
        fn bits(values: &[f64]) -> Vec<u64> {
            values.iter().map(|v| v.to_bits()).collect()
        }
        for training in [0, 2, 3] {
            let input = input(&genes, &labels, training);
            let expected = run_mle_serial(&input, &genes).unwrap();
            for threads in [1, 2, 4, 8] {
                let execution =
                    crate::execution::Execution::new(Some(threads), genes.len()).unwrap();
                let actual = run_mle_with_execution(&input, &genes, &execution).unwrap();
                match (actual.meanvar_model, expected.meanvar_model) {
                    (Some(a), Some(b)) => {
                        assert_eq!(a.slope.to_bits(), b.slope.to_bits());
                        assert_eq!(a.intercept.to_bits(), b.intercept.to_bits());
                    }
                    (None, None) => {}
                    _ => panic!("model mismatch"),
                }
                for (a, b) in actual.results.iter().zip(&expected.results) {
                    assert_eq!(a.gene_id, b.gene_id);
                    assert_eq!(a.condition_order, b.condition_order);
                    assert_eq!(a.status, b.status);
                    assert_eq!(a.attempts, b.attempts);
                    for (av, bv) in [
                        (&a.condition_beta, &b.condition_beta),
                        (&a.condition_covariance, &b.condition_covariance),
                        (&a.standard_error, &b.standard_error),
                        (&a.z_score, &b.z_score),
                        (&a.two_sided_p_value, &b.two_sided_p_value),
                        (&a.condition_fdr, &b.condition_fdr),
                        (&a.means, &b.means),
                        (&a.residuals, &b.residuals),
                    ] {
                        assert_eq!(bits(av), bits(bv));
                    }
                    assert_eq!(
                        a.efficiencies.as_deref().map(bits),
                        b.efficiencies.as_deref().map(bits)
                    );
                }
            }
        }
    }
}
