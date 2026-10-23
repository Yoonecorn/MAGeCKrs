//! Private deterministic bridge from neutral count tables to serial MLE inputs.

use std::collections::HashMap;

use super::design_matrix::DesignMatrix;
use super::serial::{SerialGeneInput, SerialMleInput};
use crate::data::ParsedCountTable;
use crate::error::MageckError;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct IngestDiagnostics {
    pub(crate) selected_samples: Vec<String>,
    pub(crate) excluded_samples: Vec<String>,
    pub(crate) sample_index_mapping: Vec<usize>,
}

#[derive(Clone, Debug)]
pub(crate) struct MleIngestBundle {
    pub(crate) genes: Vec<SerialGeneInput>,
    pub(crate) design: DesignMatrix,
    pub(crate) diagnostics: IngestDiagnostics,
}

impl MleIngestBundle {
    pub(crate) fn as_serial_mle_input<'a>(
        &'a self,
        size_factors: &'a [f64],
        genes_varmodeling: usize,
        alpha_val: f64,
        estimate_efficiency: bool,
        update_efficiency: bool,
    ) -> Result<SerialMleInput<'a>, MageckError> {
        if size_factors.len() != self.design.n_samples() {
            return Err(MageckError::Math(
                "serial MLE size factor count mismatch".into(),
            ));
        }
        if size_factors
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(MageckError::Math(
                "serial MLE size factors must be finite and positive".into(),
            ));
        }
        if !alpha_val.is_finite() || alpha_val <= 0.0 {
            return Err(MageckError::Math(
                "serial MLE alpha must be finite and positive".into(),
            ));
        }
        if genes_varmodeling > self.genes.len() {
            return Err(MageckError::Math(
                "genes_varmodeling exceeds available gene count".into(),
            ));
        }
        Ok(SerialMleInput {
            genes: &self.genes,
            n_samples: self.design.n_samples(),
            design_matrix: self.design.matrix(),
            n_conditions: self.design.n_conditions(),
            condition_labels: self.design.condition_labels(),
            size_factors,
            genes_varmodeling,
            alpha_val,
            estimate_efficiency,
            update_efficiency,
        })
    }
}

fn validate_and_transform_analytical_count(raw_count: f64) -> Result<f64, MageckError> {
    if !raw_count.is_finite() {
        return Err(MageckError::Math("raw counts must be finite".into()));
    }
    if raw_count < 0.0 {
        return Err(MageckError::Math("raw counts must be non-negative".into()));
    }
    if raw_count.fract() != 0.0 {
        return Err(MageckError::Math("raw counts must be integers".into()));
    }
    let transformed = raw_count + 1.0;
    if !transformed.is_finite() || transformed <= 0.0 || transformed.fract() != 0.0 {
        return Err(MageckError::Math(
            "analytical counts must be positive finite integers".into(),
        ));
    }
    Ok(transformed)
}

pub(crate) fn bridge_count_table_and_design(
    count_table: &ParsedCountTable,
    design: &DesignMatrix,
) -> Result<(Vec<SerialGeneInput>, IngestDiagnostics), MageckError> {
    let source_indices: HashMap<&str, usize> = count_table
        .sample_labels()
        .iter()
        .enumerate()
        .map(|(index, label)| (label.as_str(), index))
        .collect();
    let sample_index_mapping = design
        .sample_labels()
        .iter()
        .map(|sample| {
            source_indices.get(sample.as_str()).copied().ok_or_else(|| {
                MageckError::Config(format!(
                    "Sample '{sample}' specified in design matrix cannot be found in count table"
                ))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let excluded_samples = count_table
        .sample_labels()
        .iter()
        .enumerate()
        .filter(|(index, _)| !sample_index_mapping.contains(index))
        .map(|(_, label)| label.clone())
        .collect();
    let diagnostics = IngestDiagnostics {
        selected_samples: design.sample_labels().to_vec(),
        excluded_samples,
        sample_index_mapping,
    };

    let mut genes = Vec::<SerialGeneInput>::new();
    let mut gene_indices = HashMap::<&str, usize>::new();
    let mut grouped_guide_indices = Vec::<Vec<usize>>::new();
    for guide_index in 0..count_table.n_guides() {
        let gene_id = &count_table.genes()[guide_index];
        let gene_index = if let Some(index) = gene_indices.get(gene_id.as_str()) {
            *index
        } else {
            let index = genes.len();
            gene_indices.insert(gene_id, index);
            grouped_guide_indices.push(Vec::new());
            genes.push(SerialGeneInput {
                gene_id: gene_id.clone(),
                sgrna_ids: Vec::new(),
                counts: Vec::new(),
            });
            index
        };
        grouped_guide_indices[gene_index].push(guide_index);
        genes[gene_index]
            .sgrna_ids
            .push(count_table.sgrna_ids()[guide_index].clone());
    }
    // Preserve first-appearance gene and guide order without rescanning the
    // complete count table for each gene. Validation stays gene/sample/guide.
    for (gene, guide_indices) in genes.iter_mut().zip(&grouped_guide_indices) {
        let capacity = diagnostics
            .sample_index_mapping
            .len()
            .checked_mul(guide_indices.len())
            .ok_or_else(|| MageckError::Math("MLE ingest count dimension overflow".into()))?;
        gene.counts = Vec::with_capacity(capacity);
        for &sample_index in &diagnostics.sample_index_mapping {
            for &guide_index in guide_indices {
                let raw_count = count_table.counts().entry(sample_index, guide_index);
                let count = validate_and_transform_analytical_count(raw_count)?;
                gene.counts.push(count);
            }
        }
    }
    Ok((genes, diagnostics))
}

pub(crate) fn ingest_mle_input(
    count_table: &ParsedCountTable,
    design: &DesignMatrix,
) -> Result<MleIngestBundle, MageckError> {
    let (genes, diagnostics) = bridge_count_table_and_design(count_table, design)?;
    Ok(MleIngestBundle {
        genes,
        design: design.clone(),
        diagnostics,
    })
}

#[cfg(test)]
mod tests {
    use super::{ingest_mle_input, validate_and_transform_analytical_count, DesignMatrix};
    use crate::data::ParsedCountTable;
    use crate::error::MageckError;

    fn table(input: &str) -> ParsedCountTable {
        ParsedCountTable::from_reader(input.as_bytes()).unwrap()
    }
    fn design(input: &str) -> DesignMatrix {
        DesignMatrix::from_str(input).unwrap()
    }

    fn assert_math_error<T: std::fmt::Debug>(result: Result<T, MageckError>, message: &str) {
        match result {
            Err(MageckError::Math(actual)) => assert_eq!(actual, message),
            other => panic!("expected math error, got {other:?}"),
        }
    }

    #[test]
    fn fixture_01_canonical_sample_reordering() {
        let bundle = ingest_mle_input(&table("sgRNA\tGene\tTreat1\tTreat2\tDay0\na1\tG1\t30\t40\t10\na2\tG1\t31\t41\t11\n"), &design("sample\tbaseline\tbeta_1\tbeta_2\nDay0\t1\t0\t0\nTreat1\t1\t1\t0\nTreat2\t1\t0\t1\n")).unwrap();
        assert_eq!(bundle.diagnostics.sample_index_mapping, vec![2, 0, 1]);
        assert_eq!(
            bundle.genes[0].counts,
            vec![11.0, 12.0, 31.0, 32.0, 41.0, 42.0]
        );
    }

    #[test]
    fn fixture_02_extra_sample_exclusion_and_diagnostics() {
        let bundle = ingest_mle_input(
            &table("sgRNA\tGene\tExtra1\tTreat\tExtra2\tDay0\na1\tG1\t9\t20\t8\t10\n"),
            &design("sample\tbaseline\tbeta\nDay0\t1\t0\nTreat\t1\t1\n"),
        )
        .unwrap();
        assert_eq!(bundle.diagnostics.selected_samples, vec!["Day0", "Treat"]);
        assert_eq!(
            bundle.diagnostics.excluded_samples,
            vec!["Extra1", "Extra2"]
        );
        assert_eq!(bundle.genes[0].counts, vec![11.0, 21.0]);
    }

    #[test]
    fn fixture_03_gene_encounter_and_interleaved_guides() {
        let bundle = ingest_mle_input(
            &table("sgRNA\tGene\tDay0\tTreat\na1\tA\t10\t20\nb1\tB\t11\t21\na2\tA\t12\t22\n"),
            &design("sample\tbaseline\tbeta\nDay0\t1\t0\nTreat\t1\t1\n"),
        )
        .unwrap();
        assert_eq!(bundle.genes[0].gene_id, "A");
        assert_eq!(bundle.genes[0].sgrna_ids, vec!["a1", "a2"]);
        assert_eq!(bundle.genes[1].gene_id, "B");
        assert_eq!(bundle.genes[1].sgrna_ids, vec!["b1"]);
    }

    #[test]
    fn fixture_05a_missing_required_sample_error() {
        assert!(ingest_mle_input(
            &table("sgRNA\tGene\tDay0\na1\tG1\t10\n"),
            &design("sample\tbaseline\tbeta\nDay0\t1\t0\nTreat\t1\t1\n")
        )
        .is_err());
    }

    #[test]
    fn fixture_05c_case_sensitivity_mismatch_error() {
        assert!(ingest_mle_input(
            &table("sgRNA\tGene\tday0\tTreat\na1\tG1\t10\t20\n"),
            &design("sample\tbaseline\tbeta\nDay0\t1\t0\nTreat\t1\t1\n")
        )
        .is_err());
    }

    #[test]
    fn fixture_05d_zero_count_transformed_under_od_mle_5() {
        let bundle = ingest_mle_input(
            &table("sgRNA\tGene\tDay0\tTreat\na1\tG1\t0\t20\n"),
            &design("sample\tbaseline\tbeta\nDay0\t1\t0\nTreat\t1\t1\n"),
        )
        .unwrap();
        assert_eq!(bundle.genes[0].counts, vec![1.0, 21.0]);
    }

    #[test]
    fn fixture_05e_fractional_count_rejected_under_fit_boundary() {
        assert_math_error(
            ingest_mle_input(
                &table("sgRNA\tGene\tDay0\tTreat\na1\tG1\t10.5\t20\n"),
                &design("sample\tbaseline\tbeta\nDay0\t1\t0\nTreat\t1\t1\n"),
            ),
            "raw counts must be integers",
        );
    }

    #[test]
    fn defensive_analytical_count_validation() {
        assert_math_error(
            validate_and_transform_analytical_count(-1.0),
            "raw counts must be non-negative",
        );
        assert_math_error(
            validate_and_transform_analytical_count(f64::NAN),
            "raw counts must be finite",
        );
        assert_eq!(validate_and_transform_analytical_count(0.0).unwrap(), 1.0);
        assert_eq!(validate_and_transform_analytical_count(42.0).unwrap(), 43.0);
    }

    #[test]
    fn fixture_05f_size_factor_validation_errors() {
        let bundle = ingest_mle_input(
            &table("sgRNA\tGene\tDay0\tTreat\na1\tG1\t10\t20\n"),
            &design("sample\tbaseline\tbeta\nDay0\t1\t0\nTreat\t1\t1\n"),
        )
        .unwrap();
        assert!(bundle
            .as_serial_mle_input(&[1.0], 0, 0.01, false, false)
            .is_err());
        assert!(bundle
            .as_serial_mle_input(&[1.0, 0.0], 0, 0.01, false, false)
            .is_err());
    }

    #[test]
    fn fixture_06_end_to_end_bridge_to_serial_mle_input() {
        let bundle = ingest_mle_input(&table("sgRNA\tGene\tExtra_Control\tTreat_2\tDay0_Ctrl\tTreat_1\na1\tGENE_ALPHA\t999\t80\t100\t150\na2\tGENE_ALPHA\t999\t90\t120\t180\nb1\tGENE_BETA\t999\t100\t50\t40\nb2\tGENE_BETA\t999\t120\t60\t50\nb3\tGENE_BETA\t999\t140\t70\t60\nc1\tGENE_GAMMA\t999\t250\t200\t180\nc2\tGENE_GAMMA\t999\t300\t250\t220\nc3\tGENE_GAMMA\t999\t350\t300\t260\nc4\tGENE_GAMMA\t999\t400\t350\t310\n"), &design("sample\tbaseline\tbeta_1\tbeta_2\nDay0_Ctrl\t1\t0\t0\nTreat_1\t1\t1\t0\nTreat_2\t1\t0\t1\n")).unwrap();
        assert_eq!(
            bundle.diagnostics.selected_samples,
            vec!["Day0_Ctrl", "Treat_1", "Treat_2"]
        );
        assert_eq!(bundle.diagnostics.excluded_samples, vec!["Extra_Control"]);
        assert_eq!(bundle.diagnostics.sample_index_mapping, vec![2, 3, 1]);
        assert_eq!(bundle.design.entry(0, 0), 1.0);
        assert_eq!(bundle.design.entry(0, 1), 0.0);
        assert_eq!(bundle.design.condition_labels(), &["beta_1", "beta_2"]);
        let input = bundle
            .as_serial_mle_input(&[1.0, 1.1, 0.9], 3, 0.01, true, true)
            .unwrap();
        assert_eq!(
            input.design_matrix,
            &[1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0]
        );
        assert_eq!(input.genes[0].sgrna_ids, vec!["a1", "a2"]);
        assert_eq!(
            input.genes[0].counts,
            vec![101.0, 121.0, 151.0, 181.0, 81.0, 91.0]
        );
        assert_eq!(input.genes[1].sgrna_ids, vec!["b1", "b2", "b3"]);
        assert_eq!(
            input.genes[1].counts,
            vec![51.0, 61.0, 71.0, 41.0, 51.0, 61.0, 101.0, 121.0, 141.0]
        );
        assert_eq!(input.genes[2].sgrna_ids, vec!["c1", "c2", "c3", "c4"]);
        assert_eq!(
            input.genes[2].counts,
            vec![
                201.0, 251.0, 301.0, 351.0, 181.0, 221.0, 261.0, 311.0, 251.0, 301.0, 351.0, 401.0
            ]
        );
    }
}
