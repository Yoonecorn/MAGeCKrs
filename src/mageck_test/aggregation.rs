//! In-memory gene LFC and result aggregation for one selected comparison.
#![allow(dead_code)] // Consumed by M3 serialization and M4 orchestration.

use std::collections::{HashMap, HashSet};

use crate::data::ParsedCountTable;
use crate::error::MageckError;
use crate::mageck_test::rra::GeneFdrResult;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GeneLfcSummary {
    pub(crate) gene_id: String,
    pub(crate) neg_lfc: f64,
    pub(crate) pos_lfc: f64,
}

/// Complete gene-summary payload; ranking and serialization belong to M3.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AggregatedGeneResult {
    pub(crate) gene_id: String,
    pub(crate) num_sgrnas: usize,
    pub(crate) num_valid_sgrnas: usize,
    pub(crate) neg_score: f64,
    pub(crate) neg_pvalue: f64,
    pub(crate) neg_fdr: f64,
    pub(crate) neg_lfc: f64,
    pub(crate) neg_goodsgrna: usize,
    pub(crate) pos_score: f64,
    pub(crate) pos_pvalue: f64,
    pub(crate) pos_fdr: f64,
    pub(crate) pos_lfc: f64,
    pub(crate) pos_goodsgrna: usize,
}

fn finite(value: f64, context: &str) -> Result<f64, MageckError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MageckError::Math(format!(
            "LFC aggregation {context}: non-finite value {value}"
        )))
    }
}

fn validate_columns(
    indices: &[usize],
    n_samples: usize,
    direction: &str,
) -> Result<(), MageckError> {
    if indices.is_empty() {
        return Err(MageckError::Math(format!(
            "LFC aggregation: empty {direction} sample selection"
        )));
    }
    let mut seen = HashSet::new();
    for &index in indices {
        if index >= n_samples || !seen.insert(index) {
            return Err(MageckError::Math(format!(
                "LFC aggregation: invalid or duplicate {direction} sample index {index} (samples {n_samples})"
            )));
        }
    }
    Ok(())
}

fn condition_mean(
    table: &ParsedCountTable,
    columns: &[usize],
    guide: usize,
    direction: &str,
) -> Result<f64, MageckError> {
    let mut sum = 0.0;
    for &sample in columns {
        let context = format!(
            "gene {:?}, guide {guide}, {direction} sample {sample}",
            table.genes()[guide]
        );
        let raw = finite(
            table.counts().get(sample, guide).ok_or_else(|| {
                MageckError::Math(format!("LFC aggregation {context}: missing count"))
            })?,
            &context,
        )?;
        if raw < 0.0 {
            return Err(MageckError::Math(format!(
                "LFC aggregation {context}: negative raw count {raw}"
            )));
        }
        sum = finite(sum + raw, &format!("{context}, replicate sum"))?;
    }
    finite(sum / columns.len() as f64, "replicate mean")
}

/// Uses all physical guides, independently of the RRA validity mask.
/// Replicate means precede the +1 pseudocount. Both directions retain the
/// signed mean of guide log2(treatment/control) ratios. Call once per comparison.
pub(crate) fn aggregate_lfc(
    table: &ParsedCountTable,
    control_columns: &[usize],
    treatment_columns: &[usize],
) -> Result<Vec<GeneLfcSummary>, MageckError> {
    validate_columns(control_columns, table.n_samples(), "control")?;
    validate_columns(treatment_columns, table.n_samples(), "treatment")?;
    if control_columns
        .iter()
        .any(|c| treatment_columns.contains(c))
    {
        return Err(MageckError::Math(
            "LFC aggregation: overlapping control and treatment sample selections".into(),
        ));
    }
    let mut positions = HashMap::new();
    let mut genes: Vec<(String, f64, usize)> = Vec::new();
    for (guide, gene_id) in table.genes().iter().enumerate() {
        let context = format!("gene {gene_id:?}, guide {guide}");
        let control = finite(
            condition_mean(table, control_columns, guide, "control")? + 1.0,
            &context,
        )?;
        let treatment = finite(
            condition_mean(table, treatment_columns, guide, "treatment")? + 1.0,
            &context,
        )?;
        let ratio = finite(treatment / control, &format!("{context}, ratio"))?;
        let lfc = finite(ratio.log2(), &format!("{context}, log2 ratio"))?;
        let next = genes.len();
        let position = *positions.entry(gene_id).or_insert_with(|| {
            genes.push((gene_id.clone(), 0.0, 0));
            next
        });
        let (_, sum, count) = &mut genes[position];
        *sum = finite(*sum + lfc, &format!("{context}, gene sum"))?;
        *count += 1;
    }
    genes
        .into_iter()
        .map(|(gene_id, sum, count)| {
            let lfc = finite(sum / count as f64, &format!("gene {gene_id:?}, gene mean"))?;
            Ok(GeneLfcSummary {
                gene_id,
                neg_lfc: lfc,
                pos_lfc: lfc,
            })
        })
        .collect()
}

/// Match exact opaque identifiers, retaining every FDR field and source order.
pub(crate) fn aggregate_results(
    fdr: &[GeneFdrResult],
    lfc: &[GeneLfcSummary],
) -> Result<Vec<AggregatedGeneResult>, MageckError> {
    let mut by_gene = HashMap::new();
    for (index, summary) in lfc.iter().enumerate() {
        if by_gene.insert(summary.gene_id.as_str(), summary).is_some() {
            return Err(MageckError::Math(format!(
                "Result aggregation: duplicate LFC gene {:?} at index {index}",
                summary.gene_id
            )));
        }
        finite(
            summary.neg_lfc,
            &format!("gene {:?}, negative summary", summary.gene_id),
        )?;
        finite(
            summary.pos_lfc,
            &format!("gene {:?}, positive summary", summary.gene_id),
        )?;
    }
    let mut seen = HashSet::new();
    let mut output = Vec::with_capacity(fdr.len());
    for (index, result) in fdr.iter().enumerate() {
        let gene = &result.empirical;
        if !seen.insert(gene.gene_id.as_str()) {
            return Err(MageckError::Math(format!(
                "Result aggregation: duplicate FDR gene {:?} at index {index}",
                gene.gene_id
            )));
        }
        let summary = by_gene.get(gene.gene_id.as_str()).ok_or_else(|| MageckError::Math(format!(
            "Result aggregation gene-set mismatch: FDR gene {:?} at index {index} missing from LFC table", gene.gene_id
        )))?;
        output.push(AggregatedGeneResult {
            gene_id: gene.gene_id.clone(),
            num_sgrnas: gene.num_sgrnas,
            num_valid_sgrnas: gene.num_valid_sgrnas,
            neg_score: gene.neg_score,
            neg_pvalue: gene.neg_pvalue,
            neg_fdr: result.neg_fdr,
            neg_lfc: summary.neg_lfc,
            neg_goodsgrna: gene.neg_goodsgrna,
            pos_score: gene.pos_score,
            pos_pvalue: gene.pos_pvalue,
            pos_fdr: result.pos_fdr,
            pos_lfc: summary.pos_lfc,
            pos_goodsgrna: gene.pos_goodsgrna,
        });
    }
    for (index, summary) in lfc.iter().enumerate() {
        if !seen.contains(summary.gene_id.as_str()) {
            return Err(MageckError::Math(format!(
                "Result aggregation gene-set mismatch: LFC gene {:?} at index {index} missing from FDR results", summary.gene_id
            )));
        }
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mageck_test::rra::GeneEmpiricalResult;

    fn table(rows: &str) -> ParsedCountTable {
        ParsedCountTable::from_reader(format!("sgRNA\tGene\tc1\tc2\tt1\tt2\n{rows}").as_bytes())
            .unwrap()
    }

    fn fdr(id: &str) -> GeneFdrResult {
        GeneFdrResult {
            empirical: GeneEmpiricalResult {
                gene_id: id.into(),
                num_sgrnas: 3,
                num_valid_sgrnas: 2,
                neg_score: 0.12,
                pos_score: 0.34,
                neg_pvalue: 0.23,
                pos_pvalue: 0.45,
                neg_goodsgrna: 1,
                pos_goodsgrna: 2,
            },
            neg_fdr: 0.56,
            pos_fdr: 0.67,
        }
    }

    fn summary(id: &str, value: f64) -> GeneLfcSummary {
        GeneLfcSummary {
            gene_id: id.into(),
            neg_lfc: value,
            pos_lfc: value,
        }
    }

    fn math<T: std::fmt::Debug>(result: Result<T, MageckError>, needle: &str) {
        match result {
            Err(MageckError::Math(message)) => assert!(message.contains(needle), "{message}"),
            other => panic!("expected Math error containing {needle:?}: {other:?}"),
        }
    }

    #[test]
    fn replicate_mean_pseudocount_then_mean_of_guide_logs() {
        let counts = table("s1\tG\t0\t2\t2\t12\ns2\tG\t0\t0\t0\t0\n");
        let result = aggregate_lfc(&counts, &[0, 1], &[2, 3]).unwrap();
        // First guide: log2((7+1)/(1+1)) = 2. Second guide: 0.
        assert_eq!(result, vec![summary("G", 1.0)]);
    }

    #[test]
    fn first_appearance_order_and_opaque_identifiers() {
        let counts = table(
            "s1\tZ,A\t0\t0\t3\t3\ns2\tA\t3\t3\t0\t0\ns3\tZ,A\t3\t3\t0\t0\ns4\ta\t0\t0\t0\t0\n",
        );
        for _ in 0..10 {
            assert_eq!(
                aggregate_lfc(&counts, &[0, 1], &[2, 3]).unwrap(),
                vec![summary("Z,A", 0.0), summary("A", -2.0), summary("a", 0.0)]
            );
        }
    }

    #[test]
    fn comparison_columns_select_only_requested_samples() {
        let counts = table("s1\tG\t0\t3\t7\t15\n");
        assert_eq!(
            aggregate_lfc(&counts, &[0], &[2]).unwrap(),
            vec![summary("G", 3.0)]
        );
        assert_eq!(
            aggregate_lfc(&counts, &[1], &[3]).unwrap(),
            vec![summary("G", 2.0)]
        );
    }

    #[test]
    fn non_finite_replicate_sum_is_contextual_math_error() {
        let counts = table("s1\tG\t1e308\t1e308\t0\t0\n");
        math(
            aggregate_lfc(&counts, &[0, 1], &[2, 3]),
            "gene \"G\", guide 0, control sample 1, replicate sum",
        );
    }

    #[test]
    fn large_finite_counts_keep_finite_signed_lfc() {
        let counts = table("s1\tG\t1e308\t0\t0\t0\n");
        let result = aggregate_lfc(&counts, &[0], &[2]).unwrap();
        assert!(result[0].neg_lfc.is_finite());
        assert!((result[0].neg_lfc + 1e308_f64.log2()).abs() < 1e-12);
        assert_eq!(result[0].neg_lfc, result[0].pos_lfc);
    }

    #[test]
    fn invalid_sample_selections_are_rejected() {
        let counts = table("s1\tG\t0\t0\t0\t0\n");
        math(aggregate_lfc(&counts, &[], &[2]), "empty control");
        math(aggregate_lfc(&counts, &[0], &[]), "empty treatment");
        math(aggregate_lfc(&counts, &[4], &[2]), "control sample index 4");
        math(aggregate_lfc(&counts, &[0], &[2, 2]), "duplicate treatment");
        math(aggregate_lfc(&counts, &[0], &[0]), "overlapping");
    }

    #[test]
    fn every_field_merges_and_fdr_order_wins() {
        let output = aggregate_results(
            &[fdr("B"), fdr("A")],
            &[summary("A", -2.0), summary("B", 3.0)],
        )
        .unwrap();
        assert_eq!(
            output[0],
            AggregatedGeneResult {
                gene_id: "B".into(),
                num_sgrnas: 3,
                num_valid_sgrnas: 2,
                neg_score: 0.12,
                neg_pvalue: 0.23,
                neg_fdr: 0.56,
                neg_lfc: 3.0,
                neg_goodsgrna: 1,
                pos_score: 0.34,
                pos_pvalue: 0.45,
                pos_fdr: 0.67,
                pos_lfc: 3.0,
                pos_goodsgrna: 2,
            }
        );
        assert_eq!(output[1].gene_id, "A");
        assert_eq!(output[1].neg_lfc, -2.0);
        assert_eq!(output[1].pos_lfc, -2.0);
    }

    #[test]
    fn missing_and_extra_genes_are_errors_in_source_order() {
        math(
            aggregate_results(&[fdr("A")], &[summary("a", 0.0)]),
            "FDR gene \"A\"",
        );
        math(
            aggregate_results(&[fdr("A")], &[summary("A", 0.0), summary("extra", 0.0)]),
            "LFC gene \"extra\"",
        );
        math(aggregate_results(&[fdr("A")], &[]), "missing from LFC");
        math(
            aggregate_results(&[], &[summary("A", 0.0)]),
            "missing from FDR",
        );
    }

    #[test]
    fn duplicate_gene_summaries_are_errors() {
        math(
            aggregate_results(&[fdr("A"), fdr("A")], &[summary("A", 0.0)]),
            "duplicate FDR",
        );
        math(
            aggregate_results(&[fdr("A")], &[summary("A", 0.0), summary("A", 0.0)]),
            "duplicate LFC",
        );
    }

    #[test]
    fn non_finite_summary_is_rejected_in_both_directions() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            math(
                aggregate_results(&[fdr("A")], &[summary("A", value)]),
                "negative summary",
            );
            let mut lfc = summary("A", 0.0);
            lfc.pos_lfc = value;
            math(aggregate_results(&[fdr("A")], &[lfc]), "positive summary");
        }
    }

    #[test]
    fn zero_valid_gene_retains_all_physical_guide_lfc() {
        let counts = table("s1\tG\t3\t3\t0\t0\ns2\tG\t0\t0\t0\t0\n");
        let lfc = aggregate_lfc(&counts, &[0, 1], &[2, 3]).unwrap();
        let mut gene = fdr("G");
        gene.empirical.num_sgrnas = 2;
        gene.empirical.num_valid_sgrnas = 0;
        gene.empirical.neg_score = 1.0;
        gene.empirical.pos_score = 1.0;
        gene.empirical.neg_pvalue = 1.0;
        gene.empirical.pos_pvalue = 1.0;
        gene.empirical.neg_goodsgrna = 0;
        gene.empirical.pos_goodsgrna = 0;
        gene.neg_fdr = 1.0;
        gene.pos_fdr = 1.0;
        let output = aggregate_results(&[gene], &lfc).unwrap();
        assert_eq!(output[0].num_valid_sgrnas, 0);
        assert_eq!(output[0].num_sgrnas, 2);
        assert_eq!(output[0].neg_fdr, 1.0);
        assert_eq!(output[0].pos_fdr, 1.0);
        assert_eq!(output[0].neg_goodsgrna, 0);
        assert_eq!(output[0].pos_goodsgrna, 0);
        assert_eq!(output[0].neg_score, 1.0);
        assert_eq!(output[0].pos_score, 1.0);
        assert_eq!(output[0].neg_pvalue, 1.0);
        assert_eq!(output[0].pos_pvalue, 1.0);
        assert_eq!(output[0].neg_lfc, -1.0);
        assert_eq!(output[0].pos_lfc, -1.0);
    }

    #[test]
    fn empty_inputs_merge_to_empty_output() {
        assert!(aggregate_results(&[], &[]).unwrap().is_empty());
    }
}
