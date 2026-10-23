use crate::data::CountMatrix;
use crate::mageck_test::nb_model::NbModel;
use crate::stats;
use crate::MageckError;
use statrs::distribution::{ContinuousCDF, Normal};

/// Source-faithful Phase 4D sgRNA p-value result
#[derive(Debug, Clone, PartialEq)]
pub struct SgRnaPvalueResult {
    pub guide_index: usize,
    pub control_mean: f64,
    pub treatment_mean: f64,
    pub adjusted_variance: f64,
    pub p_lower: f64,
    pub p_higher: f64,
    pub p_twosided: f64,
    pub lfc: f64,
    pub score_negative: f64,
    pub score_positive: f64,
    pub rank_negative: f64,
    pub rank_positive: f64,
}

pub fn calculate_sgrna_pvalues(
    counts: &CountMatrix,
    treatment_indices: &[usize],
    nb_model: &NbModel,
) -> Result<Vec<SgRnaPvalueResult>, MageckError> {
    let n_guides = counts.n_guides();
    if n_guides == 0 {
        return Err(MageckError::Math("CountMatrix is empty".into()));
    }

    if treatment_indices.is_empty() {
        return Err(MageckError::Math("No treatment samples selected".into()));
    }

    // 1. Validate treatment indices
    for &t in treatment_indices {
        if t >= counts.n_samples() {
            return Err(MageckError::Math("Treatment index out of bounds".into()));
        }
    }

    let mut trt_dedup = treatment_indices.to_vec();
    trt_dedup.sort_unstable();
    trt_dedup.dedup();
    if trt_dedup.len() != treatment_indices.len() {
        return Err(MageckError::Math("Duplicate treatment indices".into()));
    }

    let n_trt = treatment_indices.len() as f64;

    let raw_control_means = nb_model.mean();
    let final_dispersions = nb_model.final_dispersion();

    if raw_control_means.len() != n_guides {
        return Err(MageckError::Math("NbModel mean length mismatch".into()));
    }
    if final_dispersions.len() != n_guides {
        return Err(MageckError::Math(
            "NbModel dispersion length mismatch".into(),
        ));
    }

    for &m in raw_control_means {
        if !m.is_finite() || m < 0.0 {
            return Err(MageckError::Math("Invalid raw control mean".into()));
        }
    }
    for &d in final_dispersions {
        if !d.is_finite() || d <= 0.0 {
            return Err(MageckError::Math("Invalid final dispersion".into()));
        }
    }

    // 2. Find min positive control mean
    let mut min_pos_control = f64::INFINITY;
    for &m in raw_control_means {
        if m > 0.0 && m < min_pos_control {
            min_pos_control = m;
        }
    }
    if min_pos_control.is_infinite() {
        return Err(MageckError::Math("No positive control mean found".into()));
    }

    // Clamped control means (tabc_mean)
    let tabc_mean: Vec<f64> = raw_control_means
        .iter()
        .map(|&m| {
            if m > min_pos_control {
                m
            } else {
                min_pos_control
            }
        })
        .collect();

    // 3. Calculate raw treatment means (ttmean)
    let mut ttmean = Vec::with_capacity(n_guides);
    let mut min_pos_trt = f64::INFINITY;
    for g in 0..n_guides {
        let mut sum = 0.0;
        for &t in treatment_indices {
            let val = counts.entry(t, g);
            if !val.is_finite() || val < 0.0 {
                return Err(MageckError::Math("Invalid treatment count".into()));
            }
            sum += val;
        }
        if !sum.is_finite() {
            return Err(MageckError::Math("Treatment sum is not finite".into()));
        }
        let mean = sum / n_trt;
        if !mean.is_finite() {
            return Err(MageckError::Math("Treatment mean is not finite".into()));
        }
        ttmean.push(mean);
        if mean > 0.0 && mean < min_pos_trt {
            min_pos_trt = mean;
        }
    }
    if min_pos_trt.is_infinite() {
        return Err(MageckError::Math("No positive treatment mean found".into()));
    }

    // 4. Calculate adjusted variances and p-values
    let norm = Normal::new(0.0, 1.0)
        .map_err(|_| MageckError::Math("Failed to construct Normal distribution".into()))?;

    let mut p_lower = Vec::with_capacity(n_guides);
    let mut p_higher = Vec::with_capacity(n_guides);
    let mut tabc_adjvar = Vec::with_capacity(n_guides);
    let mut lfc = Vec::with_capacity(n_guides);

    for i in 0..n_guides {
        let var0 = tabc_mean[i] + tabc_mean[i].powi(2) * final_dispersions[i];
        if !var0.is_finite() || var0 <= 0.0 {
            return Err(MageckError::Math("Invalid adjusted variance".into()));
        }
        tabc_adjvar.push(var0);

        let mean1_adj = if ttmean[i] > min_pos_trt {
            ttmean[i]
        } else {
            min_pos_trt
        };

        let t_theta = (mean1_adj - tabc_mean[i]) / var0.sqrt();
        let t_theta_0 = (0.0 - tabc_mean[i]) / var0.sqrt();

        let t_p_lower = norm.cdf(t_theta);
        let t_p_0_lower = norm.cdf(t_theta_0);

        let p_low = (t_p_lower - t_p_0_lower) / (1.0 - t_p_0_lower);

        let t_p_higher = norm.sf(t_theta);
        let p_high = t_p_higher / (1.0 - t_p_0_lower);

        if !p_low.is_finite() || p_low < -1e-9 || p_low > 1.0 + 1e-9 {
            return Err(MageckError::Math("Lower probability out of bounds".into()));
        }
        if !p_high.is_finite() || p_high < -1e-9 || p_high > 1.0 + 1e-9 {
            return Err(MageckError::Math("Higher probability out of bounds".into()));
        }

        p_lower.push(p_low.clamp(0.0, 1.0));
        p_higher.push(p_high.clamp(0.0, 1.0));

        let lfc_val = (ttmean[i] + 1.0).log2() - (tabc_mean[i] + 1.0).log2();
        if !lfc_val.is_finite() {
            return Err(MageckError::Math("LFC is not finite".into()));
        }
        lfc.push(lfc_val);
    }

    // 5. Score Reconstruction
    let mut tt_p_lower_small_nz = f64::INFINITY;
    for &p in &p_lower {
        let p01 = p * 0.1;
        if p01 > 0.0 && p01 < tt_p_lower_small_nz {
            tt_p_lower_small_nz = p01;
        }
    }
    if tt_p_lower_small_nz.is_infinite() {
        return Err(MageckError::Math(
            "No positive lower-tail probability found".into(),
        ));
    }

    let mut score_negative = Vec::with_capacity(n_guides);
    let mut score_positive = Vec::with_capacity(n_guides);
    let mut p_twosided = Vec::with_capacity(n_guides);

    for i in 0..n_guides {
        let p_low_rpnz = p_lower[i].max(tt_p_lower_small_nz);
        let ln_p = p_low_rpnz.ln();

        let var0 = tabc_adjvar[i];
        let tt_theta_0 = (ttmean[i] - tabc_mean[i]) / var0.sqrt();

        let tt_theta_low = match get_qnorm(ln_p) {
            Ok(Some(interp)) => interp,
            Ok(None) => tt_theta_0,
            Err(e) => return Err(e),
        };

        let tt_theta = if ttmean[i] < tabc_mean[i] {
            tt_theta_low
        } else {
            tt_theta_0
        };

        if !tt_theta.is_finite() {
            return Err(MageckError::Math(
                "Reconstructed score is not finite".into(),
            ));
        }

        score_negative.push(tt_theta);
        score_positive.push(-tt_theta);

        let p2 = if p_lower[i] < p_higher[i] {
            2.0 * p_lower[i]
        } else {
            2.0 * p_higher[i]
        };
        p_twosided.push(p2.clamp(0.0, 1.0));
    }

    // 6. Rank Generation
    let rank_negative = stats::assign_ranks(&score_negative);
    let rank_positive = stats::assign_ranks(&score_positive);

    for &r in &rank_negative {
        if !r.is_finite() || r <= 0.0 || r > 1.0 {
            return Err(MageckError::Math("Rank is out of bounds".into()));
        }
    }
    for &r in &rank_positive {
        if !r.is_finite() || r <= 0.0 || r > 1.0 {
            return Err(MageckError::Math("Rank is out of bounds".into()));
        }
    }

    // 7. Collect results
    let mut results = Vec::with_capacity(n_guides);
    for i in 0..n_guides {
        results.push(SgRnaPvalueResult {
            guide_index: i,
            control_mean: tabc_mean[i],
            treatment_mean: ttmean[i],
            adjusted_variance: tabc_adjvar[i],
            p_lower: p_lower[i],
            p_higher: p_higher[i],
            p_twosided: p_twosided[i],
            lfc: lfc[i],
            score_negative: score_negative[i],
            score_positive: score_positive[i],
            rank_negative: rank_negative[i],
            rank_positive: rank_positive[i],
        });
    }

    Ok(results)
}

use std::sync::OnceLock;

#[derive(Debug)]
pub struct QNormTable {
    pub logp_vec: Vec<f64>,
    pub norm_qnorm_vec: Vec<f64>,
    pub qmin: f64,
    pub qmax: f64,
    pub qstep: f64,
}

static QNORM_TABLE: OnceLock<Result<QNormTable, String>> = OnceLock::new();

pub fn parse_qnorm_table(data: &str) -> Result<QNormTable, String> {
    let mut lines = data.lines();

    let header_line = lines
        .next()
        .ok_or_else(|| "Missing header line".to_string())?;
    let mut header_parts = header_line.split_whitespace();
    let h1 = header_parts.next();
    let h2 = header_parts.next();
    let h3 = header_parts.next();

    match (h1, h2, h3) {
        (Some("LOGP"), Some("NORM_QNORM"), None) => {}
        (Some(a), Some(b), None) => return Err(format!("Invalid header fields: {} {}", a, b)),
        _ => return Err("Header must contain exactly LOGP and NORM_QNORM".to_string()),
    }

    let mut logp_vec = Vec::with_capacity(24000);
    let mut norm_qnorm_vec = Vec::with_capacity(24000);

    for (i, line) in lines.enumerate() {
        let i = i + 1; // 1-indexed (header is 0)
        let mut parts = line.split_whitespace();
        let logp_str = match parts.next() {
            Some(s) => s,
            None => {
                if line.trim().is_empty() {
                    continue; // allow trailing empty lines, but if altering shape, initializer catches it
                }
                return Err(format!("Missing logp column at row {}", i));
            }
        };
        let qnorm_str = parts
            .next()
            .ok_or_else(|| format!("Missing qnorm column at row {}", i))?;
        if parts.next().is_some() {
            return Err(format!("Extra column found at row {}", i));
        }

        let logp: f64 = logp_str
            .parse()
            .map_err(|_| format!("Invalid logp at row {}", i))?;
        let qnorm: f64 = qnorm_str
            .parse()
            .map_err(|_| format!("Invalid qnorm at row {}", i))?;

        if !logp.is_finite() {
            return Err(format!("logp is not finite at row {}", i));
        }
        if !qnorm.is_finite() {
            return Err(format!("qnorm is not finite at row {}", i));
        }

        logp_vec.push(logp);
        norm_qnorm_vec.push(qnorm);
    }

    if logp_vec.len() < 2 {
        return Err("Table has fewer than 2 data rows".to_string());
    }

    let qmin = *logp_vec
        .first()
        .ok_or_else(|| "Failed to get qmin".to_string())?;
    let qmax = *logp_vec
        .last()
        .ok_or_else(|| "Failed to get qmax".to_string())?;

    let logp_second = logp_vec
        .get(1)
        .ok_or_else(|| "Failed to get second logp".to_string())?;
    let qstep = *logp_second - qmin;

    if !qstep.is_finite() || qstep <= 0.0 {
        return Err("qstep must be finite and positive".to_string());
    }

    let mut prev_logp = qmin;
    for (i, logp) in logp_vec.iter().enumerate().skip(1) {
        if *logp <= prev_logp {
            return Err(format!("logp not strictly ascending at row {}", i));
        }
        let step = *logp - prev_logp;
        if (step - qstep).abs() > 1e-12 {
            return Err(format!(
                "Inconsistent qstep at row {}: expected {}, got {}",
                i, qstep, step
            ));
        }
        prev_logp = *logp;
    }

    Ok(QNormTable {
        logp_vec,
        norm_qnorm_vec,
        qmin,
        qmax,
        qstep,
    })
}

pub fn get_qnorm_table() -> Result<&'static QNormTable, MageckError> {
    let res = QNORM_TABLE.get_or_init(|| {
        let data = include_str!("qnorm_table.txt");
        let table = parse_qnorm_table(data)?;

        if table.logp_vec.len() != 24000 {
            return Err(format!(
                "Embedded QNorm table must have exactly 24000 rows, found {}",
                table.logp_vec.len()
            ));
        }
        if table.logp_vec.len() != table.norm_qnorm_vec.len() {
            return Err("Embedded QNorm table vector lengths mismatch".to_string());
        }
        if (table.qmin - (-24.0)).abs() > 1e-12 {
            return Err(format!(
                "Embedded QNorm table qmin mismatch: {}",
                table.qmin
            ));
        }
        if (table.qmax - (-0.001)).abs() > 1e-12 {
            return Err(format!(
                "Embedded QNorm table qmax mismatch: {}",
                table.qmax
            ));
        }
        if (table.qstep - 0.001).abs() > 1e-12 {
            return Err(format!(
                "Embedded QNorm table qstep mismatch: {}",
                table.qstep
            ));
        }

        Ok(table)
    });
    match res {
        Ok(table) => Ok(table),
        Err(e) => Err(MageckError::Math(format!(
            "QNorm table initialization failed: {}",
            e
        ))),
    }
}

pub fn get_qnorm(lp: f64) -> Result<Option<f64>, MageckError> {
    if !lp.is_finite() {
        return Err(MageckError::Math("lp is not finite".into()));
    }

    let table = get_qnorm_table()?;

    if lp < table.qmin || lp > table.qmax {
        return Ok(None);
    }

    let ind_min_raw = ((lp - table.qmin) / table.qstep).floor();
    if !ind_min_raw.is_finite() || ind_min_raw < 0.0 {
        return Err(MageckError::Math("Invalid ind_min_raw".into()));
    }

    let max_safe_index = table.logp_vec.len() as f64 + 1.0;
    if ind_min_raw > max_safe_index {
        return Err(MageckError::Math(
            "ind_min_raw exceeds safe index range".into(),
        ));
    }

    let ind_usize = ind_min_raw as usize;
    let table_len = table.logp_vec.len();

    let ind_min_clamped = if table_len >= 2 {
        ind_usize.min(table_len - 2)
    } else {
        return Err(MageckError::Math("Table too short for clamping".into()));
    };

    let lp_val = table
        .logp_vec
        .get(ind_min_clamped)
        .ok_or_else(|| MageckError::Math("Failed to get lp_val".into()))?;
    let qnorm_val = table
        .norm_qnorm_vec
        .get(ind_min_clamped)
        .ok_or_else(|| MageckError::Math("Failed to get qnorm_val".into()))?;
    let qnorm_next_val = table
        .norm_qnorm_vec
        .get(ind_min_clamped + 1)
        .ok_or_else(|| MageckError::Math("Failed to get qnorm_next_val".into()))?;

    let frac_value = ((lp - *lp_val) / table.qstep).clamp(0.0, 1.0);
    let interp_value = *qnorm_val * (1.0 - frac_value) + *qnorm_next_val * frac_value;
    Ok(Some(interp_value))
}
