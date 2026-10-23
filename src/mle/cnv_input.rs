//! Owner-approved M11A CNV matrix adapter. Scientific calculations stay in stats::cnv.

use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::io::BufRead;

use super::{gene_result::GeneResult, serial::SerialGeneInput};
use crate::error::MageckError;
use crate::stats::cnv::{fit_cnv_model, transform_log2_cn, CnvError, CnvModel};

fn schema(message: impl Into<String>) -> MageckError {
    MageckError::Config(format!("CNV matrix: {}", message.into()))
}

fn identifier(value: &str, label: &str) -> Result<(), MageckError> {
    if value.trim().is_empty() || value.contains(['\t', '\r', '\n']) {
        return Err(schema(format!("blank or invalid {label}: {value:?}")));
    }
    Ok(())
}

struct Matrix {
    profiles: Vec<String>,
    genes: Vec<String>,
    rows: HashMap<String, usize>,
    values: Vec<Vec<f64>>,
}

impl Matrix {
    fn parse(reader: impl BufRead) -> Result<Self, MageckError> {
        let mut lines = reader.lines();
        let header = lines.next().ok_or_else(|| schema("missing header"))??;
        let fields: Vec<_> = header.split('\t').collect();
        if fields.len() < 2 || fields[0] != "SYMBOL" {
            return Err(schema(
                "header requires SYMBOL followed by named log2 profiles",
            ));
        }
        let mut profiles = Vec::new();
        let mut seen = HashSet::new();
        for profile in &fields[1..] {
            identifier(profile, "profile header")?;
            if !seen.insert(*profile) {
                return Err(schema(format!("duplicate profile header {profile:?}")));
            }
            profiles.push((*profile).to_owned());
        }
        let mut genes = Vec::new();
        let mut rows = HashMap::new();
        let mut values = vec![Vec::new(); profiles.len()];
        for (index, line) in lines.enumerate() {
            let line = line?;
            let row: Vec<_> = line.split('\t').collect();
            let line_number = index + 2;
            if row.len() != fields.len() {
                return Err(schema(format!(
                    "line {line_number}: expected {} fields, got {}",
                    fields.len(),
                    row.len()
                )));
            }
            identifier(row[0], &format!("gene key on line {line_number}"))?;
            if rows.insert(row[0].to_owned(), genes.len()).is_some() {
                return Err(schema(format!(
                    "line {line_number}: duplicate gene key {:?}",
                    row[0]
                )));
            }
            for (column, profile) in profiles.iter().enumerate() {
                let value = row[column + 1].parse::<f64>().map_err(|_| schema(format!("line {line_number}, gene {:?}, profile {profile:?}: non-numeric log2 value {:?}", row[0], row[column + 1])))?;
                if !value.is_finite() {
                    return Err(schema(format!("line {line_number}, gene {:?}, profile {profile:?}: non-finite log2 value {:?}", row[0], row[column + 1])));
                }
                values[column].push(value);
            }
            genes.push(row[0].to_owned());
        }
        Ok(Self {
            profiles,
            genes,
            rows,
            values,
        })
    }
}

fn contextual(condition: &str, phase: &str, error: CnvError) -> MageckError {
    MageckError::Math(format!(
        "CNV condition {condition:?}, profile {condition:?}, {phase}: {error}"
    ))
}

fn finite(value: f64, message: &str) -> Result<f64, CnvError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CnvError::NonFiniteCalculation(message.to_owned()))
    }
}

// CN-only guards copied in the exact order/arithmetic of fit_cnv_model's
// population preflight. No regression or optimizer is run before ordinary MLE.
// Tests compare these guards directly with the unchanged kernel.
fn validate_population(cn: &[f64]) -> Result<(), CnvError> {
    if cn.is_empty() {
        return Err(CnvError::EmptyPopulation);
    }
    if cn.len() < 2 {
        return Err(CnvError::InsufficientPopulation { count: cn.len() });
    }
    let mean = finite(
        cn.iter().sum::<f64>() / cn.len() as f64,
        "Non-finite population mean",
    )?;
    let variance = cn.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>();
    if variance == 0.0 {
        return Err(CnvError::RankDeficientPopulation { variance: 0.0 });
    }
    finite(variance, "Non-finite population variance")?;
    let mut count = 0;
    let mut sum = 0.0;
    for &x in cn {
        if x <= 1.0 {
            count += 1;
            sum += x;
        }
    }
    if count < 2 {
        return Err(CnvError::InsufficientLowSideObservations { bp: 1.0, count });
    }
    let mean = finite(sum / count as f64, "Non-finite low-side mean")?;
    let mut variance = 0.0;
    for &x in cn {
        if x <= 1.0 {
            variance += (x - mean) * (x - mean);
        }
    }
    if variance == 0.0 {
        return Err(CnvError::RankDeficientLowSide {
            bp: 1.0,
            variance: 0.0,
        });
    }
    finite(variance, "Non-finite low-side variance")?;
    let upper = crate::stats::type7_percentile_999(cn).map_err(|e| match e {
        MageckError::Math(message) => CnvError::PercentileCalculationFailure(message),
        other => CnvError::PercentileCalculationFailure(other.to_string()),
    })?;
    if upper < 1.0 {
        return Err(CnvError::InvalidBounds { lower: 1.0, upper });
    }
    Ok(())
}

pub(super) struct PreparedCnv {
    conditions: Vec<String>,
    genes: Vec<String>,
    matched: Vec<usize>,
    unmatched: Vec<usize>,
    copy_numbers: Vec<Vec<f64>>,
    extra_genes: Vec<String>,
    unused_profiles: Vec<String>,
}

pub(super) fn prepare(
    reader: impl BufRead,
    conditions: &[String],
    genes: &[SerialGeneInput],
) -> Result<PreparedCnv, MageckError> {
    // Fully parse all rows/profiles, including unused data, before mapping.
    let matrix = Matrix::parse(reader)?;
    let mut seen = HashSet::new();
    if conditions.is_empty() {
        return Err(schema("no fitted conditions"));
    }
    for condition in conditions {
        identifier(condition, "design condition")?;
        if !seen.insert(condition) {
            return Err(schema(format!("duplicate design condition {condition:?}")));
        }
    }
    let mut gene_ids = HashSet::new();
    for gene in genes {
        identifier(&gene.gene_id, "fitted gene")?;
        if !gene_ids.insert(gene.gene_id.as_str()) {
            return Err(schema(format!("duplicate fitted gene {:?}", gene.gene_id)));
        }
    }
    let matched: Vec<_> = genes
        .iter()
        .enumerate()
        .filter_map(|(i, g)| matrix.rows.contains_key(&g.gene_id).then_some(i))
        .collect();
    let unmatched = genes
        .iter()
        .enumerate()
        .filter_map(|(i, g)| (!matrix.rows.contains_key(&g.gene_id)).then_some(i))
        .collect();
    let mut copy_numbers = Vec::new();
    for condition in conditions {
        let column = matrix.profiles.iter().position(|p| p == condition).ok_or_else(|| schema(format!("condition {condition:?} requires an exact profile header; no broadcast is supported")))?;
        // Transform the whole active column, including extra CNV genes.
        let transformed = transform_log2_cn(&matrix.values[column])
            .map_err(|e| contextual(condition, "log2 transformation (CNV row index)", e))?;
        let cn: Vec<_> = matched
            .iter()
            .map(|&i| transformed[matrix.rows[&genes[i].gene_id]])
            .collect();
        validate_population(&cn).map_err(|e| contextual(condition, "population validation", e))?;
        copy_numbers.push(cn);
    }
    Ok(PreparedCnv {
        conditions: conditions.to_vec(),
        genes: genes.iter().map(|g| g.gene_id.clone()).collect(),
        matched,
        unmatched,
        copy_numbers,
        extra_genes: matrix
            .genes
            .into_iter()
            .filter(|g| !gene_ids.contains(g.as_str()))
            .collect(),
        unused_profiles: matrix
            .profiles
            .into_iter()
            .filter(|p| !conditions.contains(p))
            .collect(),
    })
}

/// Serialize deterministic preparation facts before fitting; validate every
/// identifier that will enter the sidecar, including skipped/excluded entries.
pub(super) fn run_info(
    bundle: &super::ingest::MleIngestBundle,
    all_genes: &[SerialGeneInput],
    factors: &[f64],
    config: &crate::config::MleConfig,
    training_window: usize,
) -> Result<String, MageckError> {
    let mut out = String::new();
    row(
        &mut out,
        "",
        "",
        "analytical-count-transform",
        "raw_count + 1.0",
    );
    row(&mut out, "", "", "normalization", config.norm_method);
    row(
        &mut out,
        "",
        "",
        "genes-varmodeling-requested",
        config.genes_varmodeling,
    );
    row(
        &mut out,
        "",
        "",
        "variance-training-window",
        training_window,
    );
    row(
        &mut out,
        "",
        "",
        "max-sgrnas-per-gene",
        config.max_sgrnas_per_gene,
    );
    row(
        &mut out,
        "",
        "",
        "initial-scalar-alpha",
        super::run::DEFAULT_MLE_ALPHA,
    );
    for ((sample, source), factor) in bundle
        .diagnostics
        .selected_samples
        .iter()
        .zip(&bundle.diagnostics.sample_index_mapping)
        .zip(factors)
    {
        identifier(sample, "selected sample")?;
        row(&mut out, "", sample, "selected-sample", sample);
        row(&mut out, "", sample, "source-sample-index", source);
        row(
            &mut out,
            "",
            sample,
            "size-factor",
            format_args!("{factor:.17e}"),
        );
    }
    for sample in &bundle.diagnostics.excluded_samples {
        identifier(sample, "excluded sample")?;
        row(&mut out, "", sample, "excluded-sample", sample);
    }
    let mut skipped = 0;
    for gene in all_genes {
        identifier(&gene.gene_id, "input gene")?;
        for guide in &gene.sgrna_ids {
            identifier(guide, "sgRNA")?;
        }
        if gene.sgrna_ids.len() > config.max_sgrnas_per_gene {
            skipped += 1;
            row(&mut out, "", "", "skipped-gene", &gene.gene_id);
        }
    }
    row(&mut out, "", "", "skipped-genes", skipped);
    Ok(out)
}

fn row(
    out: &mut String,
    condition: &str,
    profile: &str,
    metric: &str,
    value: impl std::fmt::Display,
) {
    writeln!(out, "{condition}\t{profile}\t{metric}\t{value}").unwrap();
}

impl PreparedCnv {
    pub(super) fn correct(&self, results: &mut [GeneResult]) -> Result<String, MageckError> {
        if results.len() != self.genes.len() {
            return Err(schema("fitted result population length mismatch"));
        }
        for (gene, result) in self.genes.iter().zip(results.iter()) {
            if result.gene_id != *gene
                || result.condition_order != self.conditions
                || result.condition_beta.len() != self.conditions.len()
            {
                return Err(schema(format!(
                    "fitted result alignment mismatch for gene {gene:?}"
                )));
            }
        }
        let mut pending = Vec::new();
        for (column, condition) in self.conditions.iter().enumerate() {
            let scores: Vec<_> = self
                .matched
                .iter()
                .map(|&i| results[i].condition_beta[column])
                .collect();
            let model = fit_cnv_model(&scores, &self.copy_numbers[column])
                .map_err(|e| contextual(condition, "model fitting (matched gene index)", e))?;
            let corrected = model
                .correct_slice(&scores, &self.copy_numbers[column])
                .map_err(|e| contextual(condition, "score correction (matched gene index)", e))?;
            pending.push((model, corrected));
        }
        // Build all model diagnostics before the infallible commit of betas.
        let mut info = String::from("Condition\tProfile\tMetric\tValue\n");
        row(&mut info, "", "", "adapter", "M11A-CNV-ADAPTER-1");
        row(&mut info, "", "", "inference", "Only matched betas are CNV-corrected; covariance, standard errors, z, p-values, BH FDR and sgRNA efficiencies are retained from ordinary MLE and are not recomputed.");
        row(
            &mut info,
            "",
            "",
            "cnv-transform",
            "exp2(log2-CNV); diploid reference=1",
        );
        row(&mut info, "", "", "fitted-genes", self.genes.len());
        row(&mut info, "", "", "extra-cnv-genes", self.extra_genes.len());
        for gene in &self.extra_genes {
            row(&mut info, "", "", "extra-cnv-gene", gene);
        }
        row(
            &mut info,
            "",
            "",
            "unused-profiles",
            self.unused_profiles.len(),
        );
        for profile in &self.unused_profiles {
            row(&mut info, "", profile, "unused-profile", profile);
        }
        for (column, (model, corrected)) in pending.iter().enumerate() {
            let condition = &self.conditions[column];
            row(
                &mut info,
                condition,
                condition,
                "matched-genes",
                self.matched.len(),
            );
            row(
                &mut info,
                condition,
                condition,
                "unmatched-genes",
                self.unmatched.len(),
            );
            for &i in &self.matched {
                row(
                    &mut info,
                    condition,
                    condition,
                    "matched-gene",
                    &self.genes[i],
                );
            }
            for &i in &self.unmatched {
                row(
                    &mut info,
                    condition,
                    condition,
                    "unmatched-gene",
                    &self.genes[i],
                );
            }
            model_info(&mut info, condition, model);
            debug_assert_eq!(corrected.len(), self.matched.len());
        }
        for (column, (_, corrected)) in pending.into_iter().enumerate() {
            for (&index, value) in self.matched.iter().zip(corrected) {
                results[index].condition_beta[column] = value;
            }
        }
        Ok(info)
    }
}

fn model_info(out: &mut String, condition: &str, model: &CnvModel) {
    let d = &model.diagnostics;
    for (key, value) in [
        ("breakpoint", model.breakpoint),
        ("slope", model.slope),
        ("intercept", model.intercept),
        ("sse", model.sse),
        ("mu-ref", model.mu_ref),
        ("upper-bound", model.upper_bound),
        ("minimizer", d.minimizer),
        ("objective-value", d.objective_value),
        ("final-bracket-lower", d.final_bracket.0),
        ("final-bracket-upper", d.final_bracket.1),
    ] {
        row(out, condition, condition, key, format_args!("{value:.17e}"));
    }
    row(out, condition, condition, "iterations", d.iterations);
    row(
        out,
        condition,
        condition,
        "function-evaluations",
        d.func_evals,
    );
    row(out, condition, condition, "converged", d.converged);
    row(
        out,
        condition,
        condition,
        "stop-reason",
        format_args!("{:?}", d.stop_reason),
    );
}

#[cfg(test)]
mod tests {
    use super::super::irls::IrlsStatus;
    use super::*;
    use std::io::Cursor;

    fn genes(count: usize) -> Vec<SerialGeneInput> {
        (0..count)
            .map(|i| SerialGeneInput {
                gene_id: format!("G{i}"),
                sgrna_ids: vec![format!("s{i}")],
                counts: vec![],
            })
            .collect()
    }
    fn matrix() -> String {
        let mut text = String::from("SYMBOL\tA\tB\tunused\n");
        let a = [-1., 0., 0.5, 0.8, 1., 1.2, 1.5, 1.8, 2., 2.2];
        let b = [-2., -1., 0., 0.4, 0.7, 1., 1.1, 1.3, 1.6, 2.];
        for i in (0..10).rev() {
            writeln!(text, "G{i}\t{}\t{}\t1024", a[i], b[i]).unwrap();
        }
        text.push_str("EXTRA\t0\t0\t1024\n");
        text
    }
    fn results() -> Vec<GeneResult> {
        let scores = [-1.5, -0.8, -0.3, 0.2, 0.5, 0.9, 1.4, 1.8, 2., 2., 0.123];
        scores
            .iter()
            .enumerate()
            .map(|(i, &score)| GeneResult {
                gene_id: format!("G{i}"),
                condition_order: vec!["A".into(), "B".into()],
                condition_beta: vec![score, -score],
                condition_covariance: vec![0.04, 0.01, 0.01, 0.09],
                standard_error: vec![0.2, 0.3],
                z_score: vec![2., -1.],
                two_sided_p_value: vec![0.04, 0.2],
                condition_fdr: vec![0.08, 0.3],
                status: IrlsStatus::IterationLimitReached,
                attempts: 11,
                efficiencies: Some(vec![0.8]),
                means: vec![4.],
                residuals: vec![2.],
            })
            .collect()
    }
    #[test]
    fn full_schema_errors_precede_mapping_and_population() {
        for (input, message) in [
            ("", "missing header"),
            ("Gene\tA\nG0\t0\n", "SYMBOL"),
            ("SYMBOL\t\n", "profile header"),
            ("SYMBOL\tA\tA\n", "duplicate profile"),
            ("SYMBOL\tA\nG0\n", "expected 2 fields"),
            ("SYMBOL\tA\n \t0\n", "gene key"),
            ("SYMBOL\tA\nG0\t0\nG0\t1\n", "duplicate gene"),
            ("SYMBOL\tA\nG0\tbad\n", "non-numeric"),
            ("SYMBOL\tA\nG0\tNaN\n", "non-finite"),
            ("SYMBOL\tA\tunused\nG0\t0\tinf\n", "non-finite"),
            ("SYMBOL\tA\nG0\t0\t1\n", "expected 2 fields"),
        ] {
            let error = prepare(Cursor::new(input), &["missing".into()], &genes(3))
                .err()
                .unwrap();
            assert!(matches!(error, MageckError::Config(_)));
            assert!(error.to_string().contains(message), "{input:?}: {error}");
        }
        let bytes = b"SYMBOL\tA\nG0\t\xff\n";
        assert!(matches!(
            prepare(Cursor::new(bytes), &["A".into()], &genes(3)),
            Err(MageckError::Io(_))
        ));
    }
    #[test]
    fn exact_keys_mapping_and_full_active_transformation() {
        let order = vec!["A".into(), "B".into()];
        let prepared = prepare(Cursor::new(matrix()), &order, &genes(11)).unwrap();
        assert_eq!(prepared.matched, (0..10).collect::<Vec<_>>());
        assert_eq!(prepared.unmatched, [10]);
        assert_eq!(prepared.extra_genes, ["EXTRA"]);
        assert_eq!(prepared.copy_numbers[0][0], 0.5);
        assert_eq!(prepared.copy_numbers[1][0], 0.25);
        assert_eq!(prepared.unused_profiles, ["unused"]);
        for label in ["a", "baseline"] {
            assert!(prepare(Cursor::new(matrix()), &[label.into()], &genes(11))
                .err()
                .unwrap()
                .to_string()
                .contains("exact profile"));
        }
        let overflow = matrix().replace("EXTRA\t0\t0", "EXTRA\t0\t1024");
        let error = prepare(Cursor::new(overflow), &order, &genes(11))
            .err()
            .unwrap()
            .to_string();
        assert!(error.contains("condition \"B\", profile \"B\""));
        assert!(error.contains("overflow"));
        let renamed = matrix().replace("G5\t", "g5\t");
        assert_eq!(
            prepare(Cursor::new(renamed), &order, &genes(11))
                .unwrap()
                .unmatched,
            [5, 10]
        );
        let mut invalid = genes(3);
        invalid[1].gene_id = "G0".into();
        assert!(prepare(Cursor::new(matrix()), &order, &invalid).is_err());
    }
    #[test]
    fn cn_only_preflight_agrees_with_accepted_kernel() {
        for cn in [
            vec![],
            vec![1.],
            vec![1., 1.],
            vec![0.5, 2., 3.],
            vec![1., 1., 2.],
            vec![0.2, 0.5],
            vec![0.5, 1., 1.],
            vec![0.5, 1., 2.],
            vec![0.5, 1., 1e308],
            vec![0., f64::from_bits(1), 2.],
        ] {
            let adapter = validate_population(&cn);
            let kernel = fit_cnv_model(&vec![0.; cn.len()], &cn).map(|_| ());
            assert_eq!(adapter, kernel, "{cn:?}");
        }
    }
    #[test]
    fn accepted_numerical_fixture_correction_and_inference_preservation() {
        let prepared =
            prepare(Cursor::new(matrix()), &["A".into(), "B".into()], &genes(11)).unwrap();
        let expected = [
            -0.866460192,
            -0.8,
            -0.8248415615,
            -0.7390341308,
            -0.7670796159,
            -0.743904925,
            -0.916762739,
            -0.9416309282,
            -0.7416309282,
            -0.7416309282,
        ];
        let original = results();
        let scores_b: Vec<_> = original[..10].iter().map(|r| r.condition_beta[1]).collect();
        let model_b = fit_cnv_model(&scores_b, &prepared.copy_numbers[1]).unwrap();
        let expected_b = model_b
            .correct_slice(&scores_b, &prepared.copy_numbers[1])
            .unwrap();
        let mut corrected = original.clone();
        let info = prepared.correct(&mut corrected).unwrap();
        for (i, (before, after)) in original.iter().zip(&corrected).enumerate() {
            if i < 10 {
                assert!((after.condition_beta[0] - expected[i]).abs() <= 1e-5);
                assert_eq!(after.condition_beta[1].to_bits(), expected_b[i].to_bits());
            } else {
                assert_eq!(after.condition_beta, before.condition_beta);
            }
            let mut restored = after.clone();
            restored.condition_beta = before.condition_beta.clone();
            assert_eq!(format!("{restored:?}"), format!("{before:?}"));
        }
        assert!(info.contains("not recomputed"));
        assert!(info.contains("IntervalToleranceMet"));
        assert!(info.contains("A\tA\titerations\t7\n"));
        assert!(info.contains("A\tA\tfunction-evaluations\t8\n"));
        for _ in 0..100 {
            let mut repeat = original.clone();
            assert_eq!(prepared.correct(&mut repeat).unwrap(), info);
            for (a, b) in repeat.iter().zip(&corrected) {
                assert_eq!(
                    a.condition_beta
                        .iter()
                        .map(|v| v.to_bits())
                        .collect::<Vec<_>>(),
                    b.condition_beta
                        .iter()
                        .map(|v| v.to_bits())
                        .collect::<Vec<_>>()
                );
            }
        }
    }
    #[test]
    fn late_condition_failure_or_alignment_error_never_mutates_betas() {
        let prepared =
            prepare(Cursor::new(matrix()), &["A".into(), "B".into()], &genes(11)).unwrap();
        for case in 0..4 {
            let mut candidate = results();
            match case {
                0 => candidate[8].condition_beta[1] = f64::NAN,
                1 => candidate[8].condition_beta[1] = f64::MAX,
                2 => candidate[8].gene_id = "wrong".into(),
                _ => {
                    candidate[8].condition_order.swap(0, 1);
                }
            }
            let before = format!("{candidate:?}");
            let error = prepared.correct(&mut candidate).unwrap_err();
            if case < 2 {
                assert!(error.to_string().contains("condition \"B\", profile \"B\""));
            }
            assert_eq!(format!("{candidate:?}"), before);
        }
        assert!(prepared.correct(&mut results()[..10]).is_err());
        let mut malformed = results();
        malformed[0].condition_beta.pop();
        assert!(prepared.correct(&mut malformed).is_err());
    }
    #[test]
    fn degenerate_interval_telemetry_is_truthful() {
        let prepared = prepare(
            Cursor::new("SYMBOL\tA\tB\nG0\t-1\t-1\nG1\t0\t0\nG2\t0\t0\n"),
            &["A".into(), "B".into()],
            &genes(3),
        )
        .unwrap();
        let info = prepared.correct(&mut results()[..3]).unwrap();
        assert!(info.contains("A\tA\titerations\t0\n"));
        assert!(info.contains("A\tA\tfunction-evaluations\t1\n"));
        assert!(info.contains("DegenerateIntervalBypass"));
    }
}
