//! Observed directional gene scores from masked rank fractions.

#[path = "control_pool.rs"]
pub(crate) mod control_pool;

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use crate::error::MageckError;
use crate::mageck_test::count_table::CountTable;
use crate::prng::Xoshiro256PlusPlus;
use crate::stats::{self, DirectionalAlpha};

pub(crate) const DEFAULT_RRA_PERMUTATION_SEED: u64 = 0x0123_4567_89ab_cdef;

/// Observed gene-level RRA lo-value scores and good-sgRNA counts.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneLoValueResult {
    /// Opaque gene identifier.
    pub gene_id: String,
    /// Number of physical sgRNA observations.
    pub num_sgrnas: usize,
    /// Number of valid sgRNA observations.
    pub num_valid_sgrnas: usize,
    /// Negative observed RRA lo-value.
    pub neg_score: f64,
    /// Positive observed RRA lo-value.
    pub pos_score: f64,
    /// Valid negative ranks at or below the negative alpha.
    pub neg_goodsgrna: usize,
    /// Valid positive ranks at or below the positive alpha.
    pub pos_goodsgrna: usize,
}

#[derive(Debug, Clone, PartialEq)]
struct GeneGroup {
    gene_id: String,
    guide_indices: Vec<usize>,
}

/// Group every physical observation in first-appearance gene order.
/// Identifiers have already been validated by the caller.
fn group_sgrnas_by_gene(genes: &[String]) -> Vec<GeneGroup> {
    let mut group_index: HashMap<&str, usize> = HashMap::new();
    let mut groups: Vec<GeneGroup> = Vec::new();
    for (index, gene) in genes.iter().enumerate() {
        if let Some(&position) = group_index.get(gene.as_str()) {
            groups[position].guide_indices.push(index);
        } else {
            group_index.insert(gene.as_str(), groups.len());
            groups.push(GeneGroup {
                gene_id: gene.clone(),
                guide_indices: vec![index],
            });
        }
    }
    groups
}

fn validate_directional_alpha(alpha: DirectionalAlpha) -> Result<(), MageckError> {
    if !alpha.negative.is_finite() || !alpha.positive.is_finite() {
        return Err(MageckError::Math(format!(
            "Directional alpha values must be finite: neg={}, pos={}",
            alpha.negative, alpha.positive
        )));
    }
    if !(0.05..=0.5).contains(&alpha.negative) {
        return Err(MageckError::Math(format!(
            "Negative directional alpha {} out of bounds [0.05, 0.5]",
            alpha.negative
        )));
    }
    if !(0.01..=0.5).contains(&alpha.positive) {
        return Err(MageckError::Math(format!(
            "Positive directional alpha {} out of bounds [0.01, 0.5]",
            alpha.positive
        )));
    }

    Ok(())
}

/// Compute observed directional RRA lo-values and good-sgRNA counts.
///
/// The three input slices are aligned by physical guide row. None excludes a
/// guide in both directions. n_valid is the global valid-guide count.
pub fn compute_gene_lo_values(
    genes: &[String],
    neg_ranks: &[Option<f64>],
    pos_ranks: &[Option<f64>],
    alpha: DirectionalAlpha,
    n_valid: usize,
) -> Result<Vec<GeneLoValueResult>, MageckError> {
    // Validation order is part of the accepted P4F5 contract.
    if genes.len() != neg_ranks.len() || genes.len() != pos_ranks.len() {
        return Err(MageckError::Math(format!(
            "Length mismatch: genes ({}), neg_ranks ({}), pos_ranks ({})",
            genes.len(),
            neg_ranks.len(),
            pos_ranks.len()
        )));
    }
    if genes.is_empty() {
        return Err(MageckError::Math(
            "Input guide slices must not be empty".into(),
        ));
    }
    for (index, gene) in genes.iter().enumerate() {
        if gene.is_empty() {
            return Err(MageckError::Math(format!(
                "Row {index}: gene identifier is empty"
            )));
        }
    }
    if n_valid == 0 {
        return Err(MageckError::Math(
            "Global valid guide universe N_valid cannot be 0".into(),
        ));
    }

    let mut observed_valid_count = 0;
    for (index, (&neg, &pos)) in neg_ranks.iter().zip(pos_ranks).enumerate() {
        match (neg, pos) {
            (None, None) => {}
            (Some(neg), Some(pos)) => {
                if !neg.is_finite() || !pos.is_finite() {
                    return Err(MageckError::Math(format!(
                        "Row {index}: non-finite rank fraction (neg={neg}, pos={pos})"
                    )));
                }
                if !(0.0 < neg && neg <= 1.0 && 0.0 < pos && pos <= 1.0) {
                    return Err(MageckError::Math(format!(
                        "Row {index}: rank fraction out of bounds (0.0, 1.0] (neg={neg}, pos={pos})"
                    )));
                }
                observed_valid_count += 1;
            }
            _ => {
                return Err(MageckError::Math(format!(
                    "Row {index}: directional mask mismatch between negative and positive ranks"
                )));
            }
        }
    }
    if observed_valid_count != n_valid {
        return Err(MageckError::Math(format!(
            "Observed valid count ({observed_valid_count}) does not match n_valid ({n_valid})"
        )));
    }
    validate_directional_alpha(alpha)?;

    let mut results = Vec::new();
    for group in group_sgrnas_by_gene(genes) {
        let mut valid_neg: Vec<f64> = group
            .guide_indices
            .iter()
            .filter_map(|&index| neg_ranks[index])
            .collect();
        let mut valid_pos: Vec<f64> = group
            .guide_indices
            .iter()
            .filter_map(|&index| pos_ranks[index])
            .collect();
        let num_valid_sgrnas = valid_neg.len();
        let (neg_score, pos_score, neg_goodsgrna, pos_goodsgrna) = if num_valid_sgrnas == 0 {
            (1.0, 1.0, 0, 0)
        } else {
            valid_neg.sort_unstable_by(f64::total_cmp);
            valid_pos.sort_unstable_by(f64::total_cmp);
            let smallest_neg = valid_neg[0];
            let smallest_pos = valid_pos[0];
            valid_neg.retain(|&rank| rank <= alpha.negative);
            valid_pos.retain(|&rank| rank <= alpha.positive);
            let neg_score = if valid_neg.is_empty() {
                stats::lo_value(&[smallest_neg], num_valid_sgrnas)
            } else {
                stats::lo_value(&valid_neg, num_valid_sgrnas)
            };
            let pos_score = if valid_pos.is_empty() {
                stats::lo_value(&[smallest_pos], num_valid_sgrnas)
            } else {
                stats::lo_value(&valid_pos, num_valid_sgrnas)
            };
            if !neg_score.is_finite()
                || neg_score < 0.0
                || neg_score > 1.0
                || !pos_score.is_finite()
                || pos_score < 0.0
                || pos_score > 1.0
            {
                return Err(MageckError::Math(format!(
                        "Failed to compute lo-value for gene '{}': output is NaN or out of bounds (neg={neg_score}, pos={pos_score})",
                        group.gene_id
                    )));
            }
            let neg_goodsgrna = valid_neg.len();
            let pos_goodsgrna = valid_pos.len();
            (neg_score, pos_score, neg_goodsgrna, pos_goodsgrna)
        };
        results.push(GeneLoValueResult {
            gene_id: group.gene_id,
            num_sgrnas: group.guide_indices.len(),
            num_valid_sgrnas,
            neg_score,
            pos_score,
            neg_goodsgrna,
            pos_goodsgrna,
        });
    }
    Ok(results)
}

/// Compute observed gene scores from a count table's aligned gene rows.
pub fn compute_gene_lo_values_from_count_table(
    table: &CountTable,
    neg_ranks: &[Option<f64>],
    pos_ranks: &[Option<f64>],
    alpha: DirectionalAlpha,
    n_valid: usize,
) -> Result<Vec<GeneLoValueResult>, MageckError> {
    compute_gene_lo_values(table.genes(), neg_ranks, pos_ranks, alpha, n_valid)
}

/// Finalized upstream control-resolution counts. These are informational to
/// the numerical P4F6 engine and never change its draw sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ControlResolutionDiagnostics {
    pub total_supplied: usize,
    pub uniquely_resolved: usize,
    pub unresolved: usize,
    pub duplicated: usize,
    pub validity_excluded: usize,
    pub final_valid: usize,
}

#[derive(Debug, Clone)]
#[allow(dead_code)] // Consumed by the later in-crate RRA pipeline.
pub(crate) enum NullModelConfig<'a> {
    Ordinary,
    Control {
        neg_control_percentiles: &'a [f64],
        pos_control_percentiles: &'a [f64],
        diagnostics: &'a ControlResolutionDiagnostics,
    },
}

/// P4F5 observed scores together with P4F6 empirical significance.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneEmpiricalResult {
    pub gene_id: String,
    pub num_sgrnas: usize,
    pub num_valid_sgrnas: usize,
    pub neg_score: f64,
    pub pos_score: f64,
    pub neg_pvalue: f64,
    pub pos_pvalue: f64,
    pub neg_goodsgrna: usize,
    pub pos_goodsgrna: usize,
}

/// A complete P4F6 gene result with separate directional BH adjustments.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)] // Consumed by the later gene-level pipeline.
pub(crate) struct GeneFdrResult {
    pub(crate) empirical: GeneEmpiricalResult,
    pub(crate) neg_fdr: f64,
    pub(crate) pos_fdr: f64,
}

/// Correct eligible genes in two independent families without changing input order.
#[allow(dead_code)] // Consumed by the later gene-level pipeline.
pub(crate) fn apply_directional_bh(
    results: &[GeneEmpiricalResult],
) -> Result<Vec<GeneFdrResult>, MageckError> {
    // Validate each tier for the whole slice before starting the next tier.
    for (index, gene) in results.iter().enumerate() {
        for (direction, score) in [("negative", gene.neg_score), ("positive", gene.pos_score)] {
            if !score.is_finite() || !(0.0..=1.0).contains(&score) {
                return Err(MageckError::Math(format!(
                    "P4F7 structural validation: gene '{}' at source index {index}, {direction} observed score {score} is outside [0, 1] or non-finite",
                    gene.gene_id
                )));
            }
        }
    }

    let mut seen = HashSet::with_capacity(results.len());
    for (index, gene) in results.iter().enumerate() {
        if gene.gene_id.is_empty() || !seen.insert(gene.gene_id.as_str()) {
            return Err(MageckError::Math(format!(
                "P4F7 gene identity validation: gene '{}' at source index {index} is empty or duplicated",
                gene.gene_id
            )));
        }
    }

    for (index, gene) in results.iter().enumerate() {
        if gene.num_valid_sgrnas > gene.num_sgrnas
            || gene.neg_goodsgrna > gene.num_valid_sgrnas
            || gene.pos_goodsgrna > gene.num_valid_sgrnas
            || (gene.num_valid_sgrnas == 0
                && (gene.neg_score != 1.0
                    || gene.pos_score != 1.0
                    || gene.neg_goodsgrna != 0
                    || gene.pos_goodsgrna != 0))
        {
            return Err(MageckError::Math(format!(
                "P4F7 guide-count validation: gene '{}' at source index {index} has physical={}, valid={}, negative good={}, positive good={}, negative score={}, positive score={}",
                gene.gene_id,
                gene.num_sgrnas,
                gene.num_valid_sgrnas,
                gene.neg_goodsgrna,
                gene.pos_goodsgrna,
                gene.neg_score,
                gene.pos_score
            )));
        }
    }

    for (direction, get_pvalue) in [
        (
            "negative",
            (|gene: &GeneEmpiricalResult| gene.neg_pvalue) as fn(&GeneEmpiricalResult) -> f64,
        ),
        (
            "positive",
            (|gene: &GeneEmpiricalResult| gene.pos_pvalue) as fn(&GeneEmpiricalResult) -> f64,
        ),
    ] {
        for (index, gene) in results.iter().enumerate() {
            let value = get_pvalue(gene);
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(MageckError::Math(format!(
                    "P4F7 empirical p-value validation: gene '{}' at source index {index}, {direction} value {value} is outside [0, 1] or non-finite",
                    gene.gene_id
                )));
            }
        }
    }

    let mut eligible_indices = Vec::new();
    let mut negative = Vec::new();
    let mut positive = Vec::new();
    for (index, gene) in results.iter().enumerate() {
        if gene.num_valid_sgrnas > 0 {
            eligible_indices.push(index);
            negative.push(gene.neg_pvalue);
            positive.push(gene.pos_pvalue);
        }
    }

    let neg_adjusted = stats::bh_correction(&negative);
    let pos_adjusted = stats::bh_correction(&positive);
    for (direction, adjusted) in [("negative", &neg_adjusted), ("positive", &pos_adjusted)] {
        if adjusted.len() != eligible_indices.len() {
            return Err(MageckError::Math(format!(
                "P4F7 {direction} BH output length validation: expected {}, got {}",
                eligible_indices.len(),
                adjusted.len()
            )));
        }
        for (family_index, &value) in adjusted.iter().enumerate() {
            let source_index = eligible_indices[family_index];
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(MageckError::Math(format!(
                    "P4F7 {direction} BH output validation: gene '{}' at source index {source_index}, adjusted value {value} is outside [0, 1] or non-finite",
                    results[source_index].gene_id
                )));
            }
        }
    }

    let mut fdr_by_source = vec![None; results.len()];
    for (family_index, &source_index) in eligible_indices.iter().enumerate() {
        fdr_by_source[source_index] =
            Some((neg_adjusted[family_index], pos_adjusted[family_index]));
    }
    Ok(results
        .iter()
        .cloned()
        .zip(fdr_by_source)
        .map(|(empirical, fdr)| {
            let (neg_fdr, pos_fdr) = fdr.unwrap_or((1.0, 1.0));
            GeneFdrResult {
                empirical,
                neg_fdr,
                pos_fdr,
            }
        })
        .collect())
}

#[cfg(test)]
mod p4f7_tests {
    use super::{apply_directional_bh, GeneEmpiricalResult};

    fn gene(id: &str, valid: usize, neg: f64, pos: f64) -> GeneEmpiricalResult {
        GeneEmpiricalResult {
            gene_id: id.into(),
            num_sgrnas: 1,
            num_valid_sgrnas: valid,
            neg_score: 1.0,
            pos_score: 1.0,
            neg_pvalue: neg,
            pos_pvalue: pos,
            neg_goodsgrna: 0,
            pos_goodsgrna: 0,
        }
    }

    #[test]
    fn p4f7_separate_families_and_zero_valid_exclusion() {
        let mut input = [
            gene("a", 1, 0.01, 0.2),
            gene("zero", 0, 1.0, 1.0),
            gene("b", 1, 0.05, 0.1),
        ];
        input[0].neg_score = 0.0;
        input[2].pos_score = 0.0;
        let output = apply_directional_bh(&input).unwrap();
        assert_eq!((output[0].neg_fdr, output[0].pos_fdr), (0.02, 0.2));
        assert_eq!((output[1].neg_fdr, output[1].pos_fdr), (1.0, 1.0));
        assert_eq!((output[2].neg_fdr, output[2].pos_fdr), (0.05, 0.2));
        assert!(output
            .iter()
            .zip(input)
            .all(|(result, gene)| result.empirical == gene));
    }

    #[test]
    fn p4f7_rejects_invalid_zero_valid_input_before_bh() {
        let mut input = [gene("zero", 0, 1.0, 1.0)];
        input[0].pos_pvalue = f64::NAN;
        let error = apply_directional_bh(&input).unwrap_err().to_string();
        assert!(error.contains("empirical p-value validation"));
        assert!(error.contains("zero"));
        assert!(error.contains("positive"));
    }
}

#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)] // Kept as the deterministic internal null-pool representation.
pub(crate) struct StratifiedNullPool {
    pub(crate) strata: BTreeMap<usize, Vec<f64>>,
    pub(crate) permutation_count: usize,
}

fn validate_control_pool(pool: &[f64]) -> Result<(), MageckError> {
    for (index, &value) in pool.iter().enumerate() {
        if !value.is_finite() {
            return Err(MageckError::Math(format!(
                "Control pool contains non-finite percentile at index {index}: {value}"
            )));
        }
        if !(0.0 < value && value <= 1.0) {
            return Err(MageckError::Math(format!(
                "Control pool percentile out of domain at index {index}: {value}"
            )));
        }
    }
    Ok(())
}

fn validate_observed(gene: &GeneLoValueResult, n_valid: usize) -> Result<(), MageckError> {
    let k = gene.num_valid_sgrnas;
    if k > gene.num_sgrnas || k > n_valid {
        return Err(MageckError::Math(format!(
            "Inconsistent valid-guide cardinality for gene '{}': k={k}, physical={}, n_valid={n_valid}",
            gene.gene_id, gene.num_sgrnas
        )));
    }
    if !gene.neg_score.is_finite() || !gene.pos_score.is_finite() {
        return Err(MageckError::Math(format!(
            "Non-finite observed score for gene '{}' with valid guides",
            gene.gene_id
        )));
    }
    if !(0.0 <= gene.neg_score
        && gene.neg_score <= 1.0
        && 0.0 <= gene.pos_score
        && gene.pos_score <= 1.0)
    {
        return Err(MageckError::Math(format!(
            "Observed score out of domain for gene '{}'",
            gene.gene_id
        )));
    }
    if k == 0
        && (gene.neg_score != 1.0
            || gene.pos_score != 1.0
            || gene.neg_goodsgrna != 0
            || gene.pos_goodsgrna != 0)
    {
        return Err(MageckError::Math(format!(
            "Invalid zero-valid sentinel for gene '{}'",
            gene.gene_id
        )));
    }
    Ok(())
}

fn draw_percentiles(
    rng: &mut Xoshiro256PlusPlus,
    k: usize,
    control_pool: Option<&[f64]>,
) -> Result<Vec<f64>, MageckError> {
    let mut percentiles = Vec::with_capacity(k);
    for _ in 0..k {
        let value = match control_pool {
            Some(pool) => pool[rng.next_bounded(pool.len() as u64)? as usize],
            None => rng.next_f64_open(),
        };
        if !value.is_finite()
            || value <= 0.0
            || (control_pool.is_none() && value >= 1.0)
            || (control_pool.is_some() && value > 1.0)
        {
            return Err(MageckError::Math(format!(
                "Invalid percentile generated during permutation round: {value}"
            )));
        }
        percentiles.push(value);
    }
    Ok(percentiles)
}

fn score_percentiles(mut percentiles: Vec<f64>, k: usize, alpha: f64) -> Result<f64, MageckError> {
    percentiles.sort_unstable_by(f64::total_cmp);
    let smallest = percentiles[0];
    percentiles.retain(|&rank| rank <= alpha);
    let score = if percentiles.is_empty() {
        stats::lo_value(&[smallest], k)
    } else {
        stats::lo_value(&percentiles, k)
    };
    validate_null_score(score, k)?;
    Ok(score)
}

fn null_score(
    rng: &mut Xoshiro256PlusPlus,
    k: usize,
    control_pool: Option<&[f64]>,
    alpha: f64,
) -> Result<f64, MageckError> {
    score_percentiles(draw_percentiles(rng, k, control_pool)?, k, alpha)
}

fn validate_null_score(score: f64, k: usize) -> Result<(), MageckError> {
    if !score.is_finite() || score < 0.0 || score > 1.0 {
        return Err(MageckError::Math(format!(
            "Invalid null lo-value ({score}) generated in stratum {k}"
        )));
    }
    Ok(())
}

#[cfg(test)]
fn make_pool(
    rng: &mut Xoshiro256PlusPlus,
    required: &BTreeSet<usize>,
    permutation_count: usize,
    control_pool: Option<&[f64]>,
    alpha: f64,
) -> Result<StratifiedNullPool, MageckError> {
    make_pool_with_execution(
        rng,
        required,
        permutation_count,
        control_pool,
        alpha,
        &crate::execution::Execution::serial(),
    )
}

fn make_pool_with_execution(
    rng: &mut Xoshiro256PlusPlus,
    required: &BTreeSet<usize>,
    permutation_count: usize,
    control_pool: Option<&[f64]>,
    alpha: f64,
    execution: &crate::execution::Execution,
) -> Result<StratifiedNullPool, MageckError> {
    let mut profile = crate::execution::Profile::new("rra-pool");
    let mut strata = BTreeMap::new();
    for &k in required {
        let mut scores = Vec::with_capacity(permutation_count);
        if execution.is_parallel() {
            // At most 64 rounds and approximately 65536 percentiles per batch,
            // except one round when k itself exceeds that bound. No P*k buffer.
            let batch_size = (65536 / k).clamp(1, 64);
            while scores.len() < permutation_count {
                let rounds = batch_size.min(permutation_count - scores.len());
                let draws = (0..rounds)
                    .map(|_| draw_percentiles(rng, k, control_pool))
                    .collect();
                scores.extend(execution.map(draws, |draw| score_percentiles(draw?, k, alpha))?);
            }
        } else {
            for _ in 0..permutation_count {
                scores.push(null_score(rng, k, control_pool, alpha)?);
            }
        }
        if scores.len() != permutation_count {
            return Err(MageckError::Math(format!(
                "Stratum pool length mismatch for k={k}: expected {permutation_count}, got {}",
                scores.len()
            )));
        }
        scores.sort_unstable_by(f64::total_cmp);
        strata.insert(k, scores);
        profile.stage(&format!("stratum-{k}"));
    }
    profile.finish();
    Ok(StratifiedNullPool {
        strata,
        permutation_count,
    })
}

fn empirical_for_gene(
    gene: &GeneLoValueResult,
    pool: &StratifiedNullPool,
    observed: f64,
) -> Result<f64, MageckError> {
    let k = gene.num_valid_sgrnas;
    let scores = pool.strata.get(&k).ok_or_else(|| {
        MageckError::Math(format!(
            "Missing required stratum pool for guide cardinality k={k}"
        ))
    })?;
    stats::empirical_p_value(observed, scores, pool.permutation_count).map_err(|error| {
        MageckError::Math(format!(
            "Empirical p-value calculation failed for gene '{}': {error:?}",
            gene.gene_id
        ))
    })
}

/// Compute deterministic directional empirical significance using one PRNG stream.
/// Ordinary mode scores one shared percentile stream with each directional
/// alpha; control mode draws negative then positive pools with their alphas.
#[allow(dead_code)] // The public test pipeline has not yet been wired.
pub(crate) fn compute_empirical_significance(
    gene_results: &[GeneLoValueResult],
    null_config: &NullModelConfig<'_>,
    permutation_count: usize,
    base_seed: u64,
    n_valid: usize,
    alpha: DirectionalAlpha,
) -> Result<Vec<GeneEmpiricalResult>, MageckError> {
    compute_empirical_with_execution(
        gene_results,
        null_config,
        permutation_count,
        base_seed,
        n_valid,
        alpha,
        &crate::execution::Execution::serial(),
    )
}

pub(crate) fn compute_empirical_with_execution(
    gene_results: &[GeneLoValueResult],
    null_config: &NullModelConfig<'_>,
    permutation_count: usize,
    base_seed: u64,
    n_valid: usize,
    alpha: DirectionalAlpha,
    execution: &crate::execution::Execution,
) -> Result<Vec<GeneEmpiricalResult>, MageckError> {
    if permutation_count == 0 {
        return Err(MageckError::Config(
            "Permutation count must be greater than zero (P > 0)".into(),
        ));
    }
    if let NullModelConfig::Control {
        neg_control_percentiles,
        pos_control_percentiles,
        ..
    } = null_config
    {
        if neg_control_percentiles.is_empty() || pos_control_percentiles.is_empty() {
            return Err(MageckError::Config(
                "Requested control-null sampling but valid control pool is empty".into(),
            ));
        }
    }
    if n_valid == 0 {
        return Err(MageckError::Math(
            "Globally empty valid guide universe (n_valid == 0)".into(),
        ));
    }
    if gene_results.is_empty() {
        return Err(MageckError::Math(
            "Input gene results slice is empty".into(),
        ));
    }
    if let NullModelConfig::Control {
        neg_control_percentiles,
        pos_control_percentiles,
        diagnostics,
    } = null_config
    {
        let _ = diagnostics; // Upstream diagnostics do not alter numerical sampling.
        validate_control_pool(neg_control_percentiles)?;
        validate_control_pool(pos_control_percentiles)?;
    }

    validate_directional_alpha(alpha)?;

    let mut required = BTreeSet::new();
    for gene in gene_results {
        validate_observed(gene, n_valid)?;
        if gene.num_valid_sgrnas > 0 {
            required.insert(gene.num_valid_sgrnas);
        }
    }
    // Detect impossible draw counts before advancing the generator.
    for &k in &required {
        permutation_count
            .checked_mul(k)
            .ok_or_else(|| MageckError::Config("Permutation draw count overflow".into()))?;
    }

    let mut profile = crate::execution::Profile::new("rra-empirical");
    let mut rng = Xoshiro256PlusPlus::from_seed(base_seed);
    let (negative, positive) = match null_config {
        NullModelConfig::Ordinary => {
            // Replaying the cloned stream preserves common draws without
            // advancing the analysis stream twice when directional alphas differ.
            let mut positive_rng = rng.clone();
            let negative = make_pool_with_execution(
                &mut rng,
                &required,
                permutation_count,
                None,
                alpha.negative,
                execution,
            )?;
            let positive = if alpha.positive == alpha.negative {
                negative.clone()
            } else {
                make_pool_with_execution(
                    &mut positive_rng,
                    &required,
                    permutation_count,
                    None,
                    alpha.positive,
                    execution,
                )?
            };
            (negative, positive)
        }
        NullModelConfig::Control {
            neg_control_percentiles,
            pos_control_percentiles,
            ..
        } => {
            let negative = make_pool_with_execution(
                &mut rng,
                &required,
                permutation_count,
                Some(neg_control_percentiles),
                alpha.negative,
                execution,
            )?;
            let positive = make_pool_with_execution(
                &mut rng,
                &required,
                permutation_count,
                Some(pos_control_percentiles),
                alpha.positive,
                execution,
            )?;
            (negative, positive)
        }
    };

    profile.stage("null-pools");
    let result = execution.map(gene_results.iter().collect(), |gene| {
        let (neg_pvalue, pos_pvalue) = if gene.num_valid_sgrnas == 0 {
            (1.0, 1.0)
        } else {
            (
                empirical_for_gene(gene, &negative, gene.neg_score)?,
                empirical_for_gene(gene, &positive, gene.pos_score)?,
            )
        };
        Ok(GeneEmpiricalResult {
            gene_id: gene.gene_id.clone(),
            num_sgrnas: gene.num_sgrnas,
            num_valid_sgrnas: gene.num_valid_sgrnas,
            neg_score: gene.neg_score,
            pos_score: gene.pos_score,
            neg_pvalue,
            pos_pvalue,
            neg_goodsgrna: gene.neg_goodsgrna,
            pos_goodsgrna: gene.pos_goodsgrna,
        })
    });
    profile.stage("gene-lookup");
    profile.finish();
    result
}

#[allow(dead_code)] // Called when the later test pipeline omits an explicit seed.
pub(crate) fn compute_empirical_significance_default(
    gene_results: &[GeneLoValueResult],
    null_config: &NullModelConfig<'_>,
    permutation_count: usize,
    n_valid: usize,
    alpha: DirectionalAlpha,
) -> Result<Vec<GeneEmpiricalResult>, MageckError> {
    compute_empirical_significance(
        gene_results,
        null_config,
        permutation_count,
        DEFAULT_RRA_PERMUTATION_SEED,
        n_valid,
        alpha,
    )
}

#[cfg(test)]
mod tests {
    use super::group_sgrnas_by_gene;

    #[test]
    fn grouping_preserves_scattered_first_appearance() {
        let genes = ["A", "B", "A", "C", "B"].map(str::to_string);
        let groups = group_sgrnas_by_gene(&genes);
        assert_eq!(
            groups
                .iter()
                .map(|g| g.gene_id.as_str())
                .collect::<Vec<_>>(),
            ["A", "B", "C"]
        );
        assert_eq!(groups[0].guide_indices, [0, 2]);
        assert_eq!(groups[1].guide_indices, [1, 4]);
        assert_eq!(groups[2].guide_indices, [3]);
    }

    #[test]
    fn grouping_keeps_duplicate_physical_observations() {
        let genes = ["GENE", "GENE", "GENE"].map(str::to_string);
        let groups = group_sgrnas_by_gene(&genes);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].guide_indices, [0, 1, 2]);
    }

    #[test]
    fn grouping_treats_delimiters_as_content() {
        let genes = ["A,B", "A", "A;B", "A/B", "A,B"].map(str::to_string);
        let groups = group_sgrnas_by_gene(&genes);
        assert_eq!(
            groups
                .iter()
                .map(|g| g.gene_id.as_str())
                .collect::<Vec<_>>(),
            ["A,B", "A", "A;B", "A/B"]
        );
        assert_eq!(groups[0].guide_indices, [0, 4]);
    }

    #[test]
    fn grouping_preserves_long_first_appearance_order() {
        let genes = [
            "GAPDH", "TP53", "MYC", "BRCA1", "TP53", "EGFR", "MYC", "PTEN", "GAPDH",
        ]
        .map(str::to_string);
        let groups = group_sgrnas_by_gene(&genes);
        assert_eq!(
            groups
                .iter()
                .map(|g| g.gene_id.as_str())
                .collect::<Vec<_>>(),
            ["GAPDH", "TP53", "MYC", "BRCA1", "EGFR", "PTEN"]
        );
        assert_eq!(groups[0].guide_indices, [0, 8]);
    }
}

#[cfg(test)]
mod p4f6_tests {
    use super::{
        compute_empirical_significance, compute_empirical_significance_default, empirical_for_gene,
        make_pool, validate_null_score, ControlResolutionDiagnostics, GeneLoValueResult,
        NullModelConfig, StratifiedNullPool, DEFAULT_RRA_PERMUTATION_SEED,
    };
    use crate::error::MageckError;
    use crate::prng::Xoshiro256PlusPlus;
    use std::collections::{BTreeMap, BTreeSet};

    fn alpha() -> crate::stats::DirectionalAlpha {
        crate::stats::DirectionalAlpha {
            negative: 0.25,
            positive: 0.25,
        }
    }

    fn gene(id: &str, k: usize, score: f64) -> GeneLoValueResult {
        GeneLoValueResult {
            gene_id: id.into(),
            num_sgrnas: k,
            num_valid_sgrnas: k,
            neg_score: score,
            pos_score: score,
            neg_goodsgrna: 0,
            pos_goodsgrna: 0,
        }
    }

    fn diagnostics() -> ControlResolutionDiagnostics {
        ControlResolutionDiagnostics {
            total_supplied: 2,
            uniquely_resolved: 2,
            unresolved: 0,
            duplicated: 0,
            validity_excluded: 0,
            final_valid: 2,
        }
    }

    #[test]
    fn p4f6_ordinary_reuses_one_pool_for_both_directions_and_default_seed() {
        let genes = [gene("A", 2, 0.2), gene("B", 1, 0.7)];
        let defaulted = compute_empirical_significance_default(
            &genes,
            &NullModelConfig::Ordinary,
            30,
            5,
            alpha(),
        )
        .unwrap();
        let explicit = compute_empirical_significance(
            &genes,
            &NullModelConfig::Ordinary,
            30,
            DEFAULT_RRA_PERMUTATION_SEED,
            5,
            alpha(),
        )
        .unwrap();
        assert_eq!(defaulted, explicit);
        assert!(defaulted.iter().all(|g| g.neg_pvalue == g.pos_pvalue));
        let different =
            compute_empirical_significance(&genes, &NullModelConfig::Ordinary, 30, 1, 5, alpha())
                .unwrap();
        assert_ne!(defaulted, different);
    }

    #[test]
    fn p4f6_strata_are_ascending_sorted_exact_and_draw_once_per_ordinary_percentile() {
        let required = BTreeSet::from([1, 2, 4]);
        let mut rng = Xoshiro256PlusPlus::from_seed(0x5eed);
        let pool = make_pool(&mut rng, &required, 7, None, 0.25).unwrap();
        assert_eq!(pool.strata.keys().copied().collect::<Vec<_>>(), [1, 2, 4]);
        assert!(pool.strata.values().all(|v| v.len() == 7));
        assert!(pool
            .strata
            .values()
            .all(|v| v.windows(2).all(|w| w[0] <= w[1])));
        let mut reference = Xoshiro256PlusPlus::from_seed(0x5eed);
        for _ in 0..(7 * (1 + 2 + 4)) {
            reference.next_u64();
        }
        assert_eq!(rng, reference);
    }

    #[test]
    fn p4f6_control_accepts_one_preserves_replacement_and_direction_order() {
        let d = diagnostics();
        let config = NullModelConfig::Control {
            neg_control_percentiles: &[1.0],
            pos_control_percentiles: &[0.25, 0.75],
            diagnostics: &d,
        };
        let genes = [gene("A", 4, 1.0)];
        let observed =
            compute_empirical_significance(&genes, &config, 10, 0x5eed, 5, alpha()).unwrap();
        assert_eq!(observed[0].neg_score, 1.0);
        assert!(0.0 < observed[0].neg_pvalue && observed[0].neg_pvalue < 1.0);
        assert!(0.0 < observed[0].pos_pvalue && observed[0].pos_pvalue < 1.0);
        let required = BTreeSet::from([4]);
        let mut rng = Xoshiro256PlusPlus::from_seed(0x5eed);
        let negative = make_pool(&mut rng, &required, 10, Some(&[1.0]), 0.25).unwrap();
        assert_eq!(negative.strata[&4], vec![1.0; 10]);
        let positive = make_pool(&mut rng, &required, 10, Some(&[0.25, 0.75]), 0.25).unwrap();
        assert_eq!(
            observed[0].pos_pvalue,
            empirical_for_gene(&genes[0], &positive, 1.0).unwrap()
        );
        let mut reference = Xoshiro256PlusPlus::from_seed(0x5eed);
        for _ in 0..40 {
            reference.next_bounded(1).unwrap();
        }
        for _ in 0..40 {
            reference.next_bounded(2).unwrap();
        }
        assert_eq!(rng, reference);
    }

    #[test]
    fn p4f6_zero_valid_and_zero_physical_are_visible_without_draws() {
        let mut zero_valid = gene("masked", 0, 1.0);
        zero_valid.num_sgrnas = 3;
        let zero_physical = gene("absent", 0, 1.0);
        let genes = [zero_valid, zero_physical];
        let result =
            compute_empirical_significance(&genes, &NullModelConfig::Ordinary, 4, 9, 5, alpha())
                .unwrap();
        assert_eq!(result.len(), 2);
        assert!(result
            .iter()
            .all(|g| g.neg_pvalue == 1.0 && g.pos_pvalue == 1.0));
        let required = BTreeSet::new();
        let mut rng = Xoshiro256PlusPlus::from_seed(9);
        let before = rng.clone();
        let pools = make_pool(&mut rng, &required, 4, None, 0.25).unwrap();
        assert!(pools.strata.is_empty());
        assert_eq!(rng, before);
    }

    #[test]
    fn p4f6_validation_precedes_draws_and_preserves_error_types() {
        let valid = [gene("A", 1, 0.5)];
        assert!(matches!(
            compute_empirical_significance(&valid, &NullModelConfig::Ordinary, 0, 1, 1, alpha()),
            Err(MageckError::Config(_))
        ));
        assert!(matches!(
            compute_empirical_significance(&valid, &NullModelConfig::Ordinary, 2, 1, 0, alpha()),
            Err(MageckError::Math(_))
        ));
        let d = diagnostics();
        for value in [0.0, -0.1, 1.1, f64::NAN, f64::INFINITY] {
            let pool = [value];
            let config = NullModelConfig::Control {
                neg_control_percentiles: &pool,
                pos_control_percentiles: &[0.5],
                diagnostics: &d,
            };
            assert!(matches!(
                compute_empirical_significance(&valid, &config, 2, 1, 1, alpha()),
                Err(MageckError::Math(_))
            ));
        }
        let empty = NullModelConfig::Control {
            neg_control_percentiles: &[],
            pos_control_percentiles: &[0.5],
            diagnostics: &d,
        };
        assert!(matches!(
            compute_empirical_significance(&valid, &empty, 2, 1, 1, alpha()),
            Err(MageckError::Config(_))
        ));
        let mut bad = gene("bad", 1, f64::NAN);
        assert!(matches!(
            compute_empirical_significance(
                &[bad.clone()],
                &NullModelConfig::Ordinary,
                2,
                1,
                1,
                alpha()
            ),
            Err(MageckError::Math(_))
        ));
        bad.neg_score = 0.5;
        bad.pos_score = 1.1;
        assert!(matches!(
            compute_empirical_significance(&[bad], &NullModelConfig::Ordinary, 2, 1, 1, alpha()),
            Err(MageckError::Math(_))
        ));
    }

    #[test]
    fn p4f6_missing_and_incomplete_strata_fail_through_existing_kernel() {
        let gene = gene("A", 1, 0.5);
        let missing = StratifiedNullPool {
            strata: BTreeMap::new(),
            permutation_count: 2,
        };
        assert!(matches!(
            empirical_for_gene(&gene, &missing, gene.neg_score),
            Err(MageckError::Math(message)) if message.contains("Missing required stratum")
        ));
        let incomplete = StratifiedNullPool {
            strata: BTreeMap::from([(1, vec![0.5])]),
            permutation_count: 2,
        };
        assert!(matches!(
            empirical_for_gene(&gene, &incomplete, gene.neg_score),
            Err(MageckError::Math(message)) if message.contains("PoolLengthMismatch")
        ));
    }

    #[test]
    fn p4f6_null_score_domain_accepts_zero_and_one_and_rejects_invalid_values() {
        assert!(validate_null_score(0.0, 1).is_ok());
        assert!(validate_null_score(1.0, 1).is_ok());
        let output = compute_empirical_significance(
            &[gene("zero score", 1, 0.0)],
            &NullModelConfig::Ordinary,
            2,
            1,
            1,
            alpha(),
        )
        .unwrap();
        assert_eq!((output[0].neg_score, output[0].pos_score), (0.0, 0.0));
        for score in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                validate_null_score(score, 1),
                Err(MageckError::Math(_))
            ));
        }
    }
}

#[cfg(test)]
mod p4f15b_tests {
    use super::*;

    fn observe(ranks: &[Option<f64>], alpha: DirectionalAlpha) -> GeneLoValueResult {
        let genes = vec!["FIXTURE".to_string(); ranks.len()];
        let n_valid = ranks.iter().filter(|rank| rank.is_some()).count();
        compute_gene_lo_values(&genes, ranks, ranks, alpha, n_valid)
            .unwrap()
            .remove(0)
    }

    fn near(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < 1e-14, "{actual} != {expected}");
    }

    #[test]
    fn p4f15b_f1_both_directions_use_valid_gene_count_before_alpha() {
        let gene = observe(
            &[Some(0.90), Some(0.01), None, Some(0.02)],
            DirectionalAlpha {
                negative: 0.05,
                positive: 0.05,
            },
        );
        assert_eq!((gene.num_sgrnas, gene.num_valid_sgrnas), (4, 3));
        assert_eq!((gene.neg_goodsgrna, gene.pos_goodsgrna), (2, 2));
        near(gene.neg_score, 0.001184);
        near(gene.pos_score, 0.001184);
    }

    #[test]
    fn p4f15b_f3_both_fallbacks_keep_full_gene_count() {
        let gene = observe(
            &[Some(0.9), Some(0.8)],
            DirectionalAlpha {
                negative: 0.1,
                positive: 0.1,
            },
        );
        assert_eq!(gene.num_valid_sgrnas, 2);
        assert_eq!((gene.neg_goodsgrna, gene.pos_goodsgrna), (0, 0));
        near(gene.neg_score, 0.96);
        near(gene.pos_score, 0.96);
    }

    #[test]
    fn p4f15b_known_null_draws_apply_inclusive_alpha_and_full_n() {
        // Seed 3 draws control indices 0,1,2. Inject an unsorted known vector.
        let draws = [0.90, 0.01, 0.02];
        let mut indices = Xoshiro256PlusPlus::from_seed(3);
        for expected in 0..3 {
            assert_eq!(indices.next_bounded(3).unwrap(), expected);
        }
        let gene = observe(
            &[Some(0.90), Some(0.01), Some(0.02)],
            DirectionalAlpha {
                negative: 0.05,
                positive: 0.01,
            },
        );
        for (alpha, expected) in [(0.05, gene.neg_score), (0.01, gene.pos_score)] {
            let mut rng = Xoshiro256PlusPlus::from_seed(3);
            let actual = null_score(&mut rng, 3, Some(&draws), alpha).unwrap();
            near(actual, expected);
            assert_eq!(rng, indices);
        }
        near(gene.neg_score, 0.001184);
        // The unfiltered minimum is 0.001184; this value proves alpha is used.
        near(gene.pos_score, 0.029701);
    }

    #[test]
    fn p4f15b_known_null_empty_selection_uses_first_rank_full_n() {
        // Seed 3 draws control indices 0,1, preserving the exact two draws.
        let draws = [0.8, 0.9];
        let mut rng = Xoshiro256PlusPlus::from_seed(3);
        let mut reference = rng.clone();
        assert_eq!(reference.next_bounded(2).unwrap(), 0);
        assert_eq!(reference.next_bounded(2).unwrap(), 1);
        near(null_score(&mut rng, 2, Some(&draws), 0.1).unwrap(), 0.96);
        assert_eq!(rng, reference);
    }

    #[test]
    fn p4f15b_ordinary_threads_distinct_alphas_over_shared_draws() {
        let alpha = DirectionalAlpha {
            negative: 0.05,
            positive: 0.5,
        };
        let gene = observe(&[Some(0.1), Some(0.2), Some(0.3)], alpha);
        let actual = compute_empirical_significance(
            std::slice::from_ref(&gene),
            &NullModelConfig::Ordinary,
            31,
            3,
            3,
            alpha,
        )
        .unwrap()
        .remove(0);
        let required = BTreeSet::from([3]);
        let mut negative_rng = Xoshiro256PlusPlus::from_seed(3);
        let mut positive_rng = negative_rng.clone();
        let negative = make_pool(&mut negative_rng, &required, 31, None, alpha.negative).unwrap();
        let positive = make_pool(&mut positive_rng, &required, 31, None, alpha.positive).unwrap();
        assert_eq!(negative_rng, positive_rng);
        assert_ne!(negative, positive);
        assert_eq!(
            actual.neg_pvalue,
            empirical_for_gene(&gene, &negative, gene.neg_score).unwrap()
        );
        assert_eq!(
            actual.pos_pvalue,
            empirical_for_gene(&gene, &positive, gene.pos_score).unwrap()
        );
    }

    #[test]
    fn p4f15b_control_threads_distinct_alphas_in_original_draw_order() {
        let alpha = DirectionalAlpha {
            negative: 0.05,
            positive: 0.5,
        };
        let gene = observe(&[Some(0.1), Some(0.2), Some(0.3)], alpha);
        let diagnostics = ControlResolutionDiagnostics {
            total_supplied: 3,
            uniquely_resolved: 3,
            unresolved: 0,
            duplicated: 0,
            validity_excluded: 0,
            final_valid: 3,
        };
        let draws = [0.01, 0.2, 0.8];
        let config = NullModelConfig::Control {
            neg_control_percentiles: &draws,
            pos_control_percentiles: &draws,
            diagnostics: &diagnostics,
        };
        let actual =
            compute_empirical_significance(std::slice::from_ref(&gene), &config, 17, 3, 3, alpha)
                .unwrap()
                .remove(0);
        let required = BTreeSet::from([3]);
        let mut rng = Xoshiro256PlusPlus::from_seed(3);
        let negative = make_pool(&mut rng, &required, 17, Some(&draws), alpha.negative).unwrap();
        let positive = make_pool(&mut rng, &required, 17, Some(&draws), alpha.positive).unwrap();
        assert_eq!(
            actual.neg_pvalue,
            empirical_for_gene(&gene, &negative, gene.neg_score).unwrap()
        );
        assert_eq!(
            actual.pos_pvalue,
            empirical_for_gene(&gene, &positive, gene.pos_score).unwrap()
        );
    }

    #[test]
    fn p4f15b_empirical_alpha_validation_matches_observed_scoring() {
        let valid = observe(
            &[Some(0.2)],
            DirectionalAlpha {
                negative: 0.25,
                positive: 0.25,
            },
        );
        for alpha in [
            DirectionalAlpha {
                negative: f64::NAN,
                positive: 0.25,
            },
            DirectionalAlpha {
                negative: 1.0,
                positive: 0.25,
            },
            DirectionalAlpha {
                negative: 0.25,
                positive: 0.0,
            },
        ] {
            assert!(compute_empirical_significance(
                std::slice::from_ref(&valid),
                &NullModelConfig::Ordinary,
                3,
                3,
                1,
                alpha,
            )
            .is_err());
        }
    }
}

#[cfg(test)]
mod threading_tests {
    use super::*;
    #[test]
    fn bounded_batches_preserve_every_null_score_and_next_prng_draw() {
        let required = [1, 4, 999, 1000].into_iter().collect();
        let controls = [0.01, 0.1, 0.2, 0.5, 1.0];
        for control in [None, Some(controls.as_slice())] {
            for seed in [0, u64::MAX] {
                let mut serial_rng = Xoshiro256PlusPlus::from_seed(seed);
                let expected = make_pool(&mut serial_rng, &required, 129, control, 0.05).unwrap();
                for threads in [2, 4, 8] {
                    let execution = crate::execution::Execution::new(Some(threads), 129).unwrap();
                    let mut rng = Xoshiro256PlusPlus::from_seed(seed);
                    let actual = make_pool_with_execution(
                        &mut rng, &required, 129, control, 0.05, &execution,
                    )
                    .unwrap();
                    assert_eq!(actual, expected);
                    let mut expected_rng = serial_rng.clone();
                    for _ in 0..8 {
                        assert_eq!(rng.next_u64(), expected_rng.next_u64());
                    }
                }
            }
        }
    }
}
