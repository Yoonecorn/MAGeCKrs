//! Normalization size factor estimation for classical MAGeCK MLE.

use std::collections::HashSet;

use super::design_matrix::DesignMatrix;
use super::ingest::{ingest_mle_input, MleIngestBundle};
use super::serial::SerialGeneInput;
use crate::cli::NormMethod;
use crate::data::ParsedCountTable;
use crate::error::MageckError;

pub(crate) fn compute_mle_size_factors(
    count_table: &ParsedCountTable,
    design: &DesignMatrix,
    method: NormMethod,
    control_sgrnas: Option<&[String]>,
) -> Result<Vec<f64>, MageckError> {
    let bundle = ingest_mle_input(count_table, design)?;
    compute_mle_size_factors_for_bundle(&bundle, method, control_sgrnas)
}

pub(crate) fn compute_mle_size_factors_for_bundle(
    bundle: &MleIngestBundle,
    method: NormMethod,
    control_sgrnas: Option<&[String]>,
) -> Result<Vec<f64>, MageckError> {
    let factors = match method {
        NormMethod::Median => normalize_median(bundle)?,
        NormMethod::Total => normalize_total(bundle)?,
        NormMethod::Control => normalize_control(
            bundle,
            control_sgrnas.ok_or_else(|| {
                MageckError::Config("Control normalization requires a control sgRNA list".into())
            })?,
        )?,
        NormMethod::None => normalize_none(&bundle.design)?,
    };
    validate_factors(&factors)?;
    Ok(factors)
}

pub(crate) fn normalize_median(bundle: &MleIngestBundle) -> Result<Vec<f64>, MageckError> {
    normalize_median_genes(&bundle.genes, &bundle.design)
}

pub(crate) fn normalize_total(bundle: &MleIngestBundle) -> Result<Vec<f64>, MageckError> {
    let n_samples = bundle.design.n_samples();
    let mut totals = vec![0.0; n_samples];
    for gene in &bundle.genes {
        let n_guides = gene.sgrna_ids.len();
        for (sample, total) in totals.iter_mut().enumerate() {
            let start = sample * n_guides;
            *total += gene.counts[start..start + n_guides].iter().sum::<f64>();
        }
    }
    for (sample, total) in totals.iter().enumerate() {
        if *total == 0.0 {
            return Err(MageckError::Math(format!(
                "Sample '{}' has zero total count; cannot compute total-count size factor",
                bundle.design.sample_labels()[sample]
            )));
        }
    }
    let average = totals.iter().sum::<f64>() / n_samples as f64;
    Ok(totals.into_iter().map(|total| total / average).collect())
}

pub(crate) fn normalize_control(
    bundle: &MleIngestBundle,
    control_sgrnas: &[String],
) -> Result<Vec<f64>, MageckError> {
    let controls: HashSet<&str> = control_sgrnas.iter().map(String::as_str).collect();
    let mut selected = Vec::new();
    for gene in &bundle.genes {
        let n_guides = gene.sgrna_ids.len();
        for (guide_index, guide) in gene.sgrna_ids.iter().enumerate() {
            if controls.contains(guide.as_str()) {
                selected.push(SerialGeneInput {
                    gene_id: gene.gene_id.clone(),
                    sgrna_ids: vec![guide.clone()],
                    counts: (0..bundle.design.n_samples())
                        .map(|sample| gene.counts[sample * n_guides + guide_index])
                        .collect(),
                });
            }
        }
    }
    if selected.is_empty() {
        return Err(MageckError::Config(
            "No control sgRNAs found in count table".into(),
        ));
    }
    normalize_median_genes(&selected, &bundle.design)
}

pub(crate) fn normalize_none(design: &DesignMatrix) -> Result<Vec<f64>, MageckError> {
    Ok(vec![1.0; design.n_samples()])
}

fn normalize_median_genes(
    genes: &[SerialGeneInput],
    design: &DesignMatrix,
) -> Result<Vec<f64>, MageckError> {
    let n_samples = design.n_samples();
    let mut ratios = vec![Vec::new(); n_samples];
    for gene in genes {
        let n_guides = gene.sgrna_ids.len();
        for guide in 0..n_guides {
            let values: Vec<f64> = (0..n_samples)
                .map(|sample| gene.counts[sample * n_guides + guide])
                .collect();
            if values.iter().sum::<f64>() == 0.0 {
                continue;
            }
            let geometric_mean = (values.iter().map(|value| (value + 1.0).ln()).sum::<f64>()
                / n_samples as f64)
                .exp();
            for (sample, value) in values.into_iter().enumerate() {
                ratios[sample].push(value / geometric_mean);
            }
        }
    }
    if ratios.iter().any(Vec::is_empty) {
        return Err(MageckError::Math(
            "Every guide has zero total count across selected samples; cannot compute median-ratio size factors"
                .into(),
        ));
    }
    let mut factors = Vec::with_capacity(n_samples);
    let mut requires_total_fallback = false;
    for sample_ratios in &mut ratios {
        sample_ratios.sort_by(f64::total_cmp);
        let factor = sample_ratios[sample_ratios.len() / 2];
        requires_total_fallback |= factor <= 0.0;
        factors.push(factor);
    }
    if requires_total_fallback {
        let fallback = MleIngestBundle {
            genes: genes.to_vec(),
            design: design.clone(),
            diagnostics: super::ingest::IngestDiagnostics {
                selected_samples: design.sample_labels().to_vec(),
                excluded_samples: Vec::new(),
                sample_index_mapping: Vec::new(),
            },
        };
        return normalize_total(&fallback);
    }
    Ok(factors)
}

fn validate_factors(factors: &[f64]) -> Result<(), MageckError> {
    if factors
        .iter()
        .any(|factor| !factor.is_finite() || *factor <= 0.0)
    {
        return Err(MageckError::Math(
            "Computed normalization size factor is not positive finite".into(),
        ));
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::excessive_precision)]
mod tests {
    use super::*;
    use crate::error::MageckError;

    fn table(input: &str) -> ParsedCountTable {
        ParsedCountTable::from_reader(input.as_bytes()).unwrap()
    }

    fn design(input: &str) -> DesignMatrix {
        DesignMatrix::from_str(input).unwrap()
    }

    fn assert_close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.iter().zip(expected) {
            assert!((actual - expected).abs() <= 1e-12 * expected.abs().max(1.0));
        }
    }

    fn assert_error<T: std::fmt::Debug>(result: Result<T, MageckError>, expected: MageckError) {
        match (result, expected) {
            (Err(MageckError::Config(actual)), MageckError::Config(expected))
            | (Err(MageckError::Math(actual)), MageckError::Math(expected)) => {
                assert_eq!(actual, expected);
            }
            (other, _) => panic!("unexpected result: {other:?}"),
        }
    }

    #[test]
    fn test_fixture_01_median_normalization_reordered_samples() {
        let factors = compute_mle_size_factors(&table("sgRNA\tGene\tExtra_Sample\tTreat_2\tDay0_Ctrl\tTreat_1\ng1\tGENE_A\t999\t200\t100\t150\ng2\tGENE_A\t999\t120\t50\t80\ng3\tGENE_B\t999\t400\t200\t300\ng4\tGENE_B\t999\t160\t80\t120\n"), &design("sample\tbaseline\tbeta_1\tbeta_2\nDay0_Ctrl\t1\t0\t0\nTreat_1\t1\t1\t0\nTreat_2\t1\t0\t1\n"), NormMethod::Median, None).unwrap();
        assert_close(
            &factors,
            &[0.6903293380996394, 1.0320765351786687, 1.3802240294504677],
        );
    }

    #[test]
    fn test_fixture_02_total_normalization_excluded_samples() {
        let factors = compute_mle_size_factors(&table("sgRNA\tGene\tS_Extra\tS2_High\tS0_Base\tS1_Med\ng1\tG1\t9999\t300\t120\t180\ng2\tG2\t9999\t600\t240\t360\ng3\tG3\t9999\t900\t360\t540\n"), &design("sample\tbaseline\tcond_1\tcond_2\nS0_Base\t1\t0\t0\nS1_Med\t1\t1\t0\nS2_High\t1\t0\t1\n"), NormMethod::Total, None).unwrap();
        assert_close(
            &factors,
            &[0.6009975062344139, 0.9002493765586035, 1.4987531172069826],
        );
    }

    #[test]
    fn test_fixture_03_control_normalization_deduplication() {
        let controls = vec![
            "ctrl_1".into(),
            "ctrl_2".into(),
            "ctrl_1".into(),
            "ghost".into(),
        ];
        let factors = compute_mle_size_factors(&table("sgRNA\tGene\tDay0\tTreat1\tTreat2\ntarget_1\tTARGET_A\t100\t200\t400\ntarget_2\tTARGET_B\t150\t300\t600\nctrl_1\tCONTROL\t100\t200\t300\nctrl_2\tCONTROL\t80\t160\t240\nctrl_3\tCONTROL\t120\t240\t360\n"), &design("sample\tbaseline\tbeta_1\tbeta_2\nDay0\t1\t0\t0\nTreat1\t1\t1\t0\nTreat2\t1\t0\t1\n"), NormMethod::Control, Some(&controls)).unwrap();
        assert_close(
            &factors,
            &[0.5491216233532654, 1.0928063989505579, 1.6364911745478505],
        );
    }

    #[test]
    fn test_fixture_04_none_normalization() {
        let factors = compute_mle_size_factors(
            &table("sgRNA\tGene\tS0\tS1\tS2\ng\tG\t1\t2\t3\n"),
            &design("sample\tbaseline\tb1\tb2\nS0\t1\t0\t0\nS1\t1\t1\t0\nS2\t1\t0\t1\n"),
            NormMethod::None,
            None,
        )
        .unwrap();
        assert_eq!(factors, vec![1.0; 3]);
    }

    #[test]
    fn test_fixture_05_zero_containing_guides() {
        let factors = compute_mle_size_factors(&table("sgRNA\tGene\tDay0\tTreat1\tTreat2\ng1\tG1\t0\t50\t100\ng2\tG2\t30\t0\t90\ng3\tG3\t40\t80\t0\ng4\tG4\t60\t120\t180\ng5\tG_ZERO\t0\t0\t0\n"), &design("sample\tbaseline\tbeta_1\tbeta_2\nDay0\t1\t0\t0\nTreat1\t1\t1\t0\nTreat2\t1\t0\t1\n"), NormMethod::Median, None).unwrap();
        assert_close(
            &factors,
            &[0.5483470763870151, 1.0877048564398166, 1.6270626364926184],
        );
    }

    #[test]
    fn test_fixture_06_control_errors() {
        let count = table("sgRNA\tGene\tS0\tS1\ng\tG\t1\t2\n");
        let matrix = design("sample\tbaseline\tbeta\nS0\t1\t0\nS1\t1\t1\n");
        assert_error(
            compute_mle_size_factors(&count, &matrix, NormMethod::Control, None),
            MageckError::Config("Control normalization requires a control sgRNA list".into()),
        );
        assert_error(
            compute_mle_size_factors(
                &count,
                &matrix,
                NormMethod::Control,
                Some(&["ghost".into()]),
            ),
            MageckError::Config("No control sgRNAs found in count table".into()),
        );
    }

    #[test]
    fn test_fixture_07_degenerate_sample_states() {
        let matrix = design("sample\tbaseline\tbeta\nS0\t1\t0\nS1\t1\t1\n");
        let empty = MleIngestBundle {
            genes: Vec::new(),
            design: matrix.clone(),
            diagnostics: super::super::ingest::IngestDiagnostics {
                selected_samples: vec!["S0".into(), "S1".into()],
                excluded_samples: Vec::new(),
                sample_index_mapping: Vec::new(),
            },
        };
        assert_error(
            normalize_total(&empty),
            MageckError::Math(
                "Sample 'S0' has zero total count; cannot compute total-count size factor".into(),
            ),
        );
        assert_error(normalize_median(&empty), MageckError::Math("Every guide has zero total count across selected samples; cannot compute median-ratio size factors".into()));
        assert_error(compute_mle_size_factors(&table("sgRNA\tGene\tS0\tS1\ng\tG\t1\t2\n"), &design("sample\tbaseline\tbeta_1\tbeta_2\nS0\t1\t0\t0\nS1\t1\t1\t0\nMissingSample\t1\t0\t1\n"), NormMethod::Control, None), MageckError::Config("Sample 'MissingSample' specified in design matrix cannot be found in count table".into()));
    }

    #[test]
    fn output_factor_validation_and_determinism() {
        assert_error(
            validate_factors(&[0.0]),
            MageckError::Math("Computed normalization size factor is not positive finite".into()),
        );
        assert_error(
            validate_factors(&[f64::NAN]),
            MageckError::Math("Computed normalization size factor is not positive finite".into()),
        );
        let count = table("sgRNA\tGene\tS0\tS1\ng\tG\t1\t2\n");
        let matrix = design("sample\tbaseline\tbeta\nS0\t1\t0\nS1\t1\t1\n");
        let first = compute_mle_size_factors(&count, &matrix, NormMethod::Median, None).unwrap();
        let second = compute_mle_size_factors(&count, &matrix, NormMethod::Median, None).unwrap();
        assert_eq!(first, second);
    }
}
