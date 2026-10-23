//! Bounded, single-gene translation of MAGeCK2's classical EM/IRLS loop.

use super::likelihood::nb_log_likelihood_terms;
use crate::error::MageckError;

const MAX_ITERATIONS: usize = 1000;
const MAX_ABS_BETA: f64 = 100.0;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum IrlsStatus {
    Converged,
    ConvergedSameModelRetry,
    NonfiniteIterate,
    SingularSystem,
    InverseResidualFailed,
    InferenceInvalid,
    IterationLimitReached,
    BetaRejected,
}
impl IrlsStatus {
    pub(crate) fn converged(self) -> bool {
        matches!(self, Self::Converged | Self::ConvergedSameModelRetry)
    }
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Converged => "CONVERGED_PRIMARY",
            Self::ConvergedSameModelRetry => "CONVERGED_SAME_MODEL_RETRY",
            Self::IterationLimitReached => "ITERATION_LIMIT",
            Self::BetaRejected => "COEFFICIENT_BOUND",
            Self::NonfiniteIterate => "NONFINITE_ITERATE",
            Self::SingularSystem => "SINGULAR_SYSTEM",
            Self::InverseResidualFailed => "INVERSE_RESIDUAL_FAILED",
            Self::InferenceInvalid => "INFERENCE_INVALID",
        }
    }
}
#[derive(Clone, Debug, Default)]
pub(crate) struct UpdateDiagnostics {
    pub(crate) accepted: usize,
    pub(crate) d: Option<f64>,
    pub(crate) r: Option<f64>,
    pub(crate) streak: usize,
    pub(crate) failure: Option<String>,
}
#[derive(Clone, Copy, Debug)]
pub(crate) enum DispersionInput<'a> {
    Scalar(f64),
    Vector(&'a [f64]),
}
#[derive(Debug)]
pub(crate) struct IrlsInput<'a> {
    pub(crate) relaxation: f64,
    pub(crate) n_samples: usize,
    pub(crate) n_guides: usize,
    pub(crate) n_conditions: usize,
    pub(crate) counts: &'a [f64],
    pub(crate) extended_design: &'a [f64],
    pub(crate) beta_init: &'a [f64],
    pub(crate) size_factors: &'a [f64],
    pub(crate) dispersion: DispersionInput<'a>,
    pub(crate) alpha_val: f64,
    pub(crate) guide_efficiencies: Option<&'a [f64]>,
    pub(crate) update_efficiency: bool,
}
#[derive(Clone, Debug)]
pub(crate) struct IrlsResult {
    pub(crate) diagnostics: UpdateDiagnostics,
    pub(crate) status: IrlsStatus,
    pub(crate) iterations: usize,
    pub(crate) beta: Vec<f64>,
    pub(crate) efficiencies: Option<Vec<f64>>,
    pub(crate) condition_se: Vec<f64>,
    pub(crate) means: Vec<f64>,
    pub(crate) residuals: Vec<f64>,
}

fn numerical(context: impl Into<String>) -> MageckError {
    MageckError::Math(context.into())
}
fn checked_product(a: usize, b: usize, label: &str) -> Result<usize, MageckError> {
    a.checked_mul(b)
        .ok_or_else(|| numerical(format!("IRLS {label} dimension overflow")))
}
fn finite(v: f64, what: &str) -> Result<(), MageckError> {
    if v.is_finite() {
        Ok(())
    } else {
        Err(numerical(format!("IRLS non-finite {what}")))
    }
}

impl IrlsInput<'_> {
    fn validate(&self) -> Result<(usize, usize), MageckError> {
        if ![1.0, 0.5, 0.25].contains(&self.relaxation) {
            return Err(numerical("IRLS invalid relaxation"));
        }
        if self.beta_init.iter().any(|x| x.abs() > MAX_ABS_BETA) {
            return Err(numerical("IRLS initialization exceeds coefficient bound"));
        }
        if self.n_samples == 0 || self.n_guides == 0 || self.n_conditions == 0 {
            return Err(numerical("IRLS dimensions must be non-zero"));
        }
        let sample_rows = self
            .n_samples
            .checked_mul(2)
            .and_then(|x| x.checked_sub(1))
            .ok_or_else(|| numerical("IRLS observation dimension overflow"))?;
        let n = checked_product(sample_rows, self.n_guides, "observation")?;
        let p = self
            .n_guides
            .checked_add(self.n_conditions)
            .ok_or_else(|| numerical("IRLS coefficient dimension overflow"))?;
        if self.counts.len() != n
            || self.size_factors.len() != n
            || self.beta_init.len() != p
            || self.extended_design.len() != checked_product(n, p, "design")?
        {
            return Err(numerical("IRLS input shape mismatch"));
        }
        for &x in self.counts {
            if !x.is_finite() || x < 0.0 || x.fract() != 0.0 {
                return Err(numerical(
                    "IRLS counts must be finite non-negative integers",
                ));
            }
        }
        for &x in self.size_factors {
            if !x.is_finite() || x <= 0.0 {
                return Err(numerical("IRLS size factors must be finite and positive"));
            }
        }
        for &x in self.beta_init {
            finite(x, "starting beta")?;
        }
        for &x in self.extended_design {
            finite(x, "design entry")?;
        }
        if !self.alpha_val.is_finite() || self.alpha_val < 0.0 {
            return Err(numerical(
                "IRLS ridge alpha must be finite and non-negative",
            ));
        }
        match &self.dispersion {
            DispersionInput::Scalar(x) => {
                if !x.is_finite() || *x <= 0.0 {
                    return Err(numerical("IRLS dispersion must be finite and positive"));
                }
            }
            DispersionInput::Vector(v) => {
                if v.len() != n {
                    return Err(numerical("IRLS dispersion shape mismatch"));
                }
                for &x in *v {
                    if !x.is_finite() || x <= 0.0 {
                        return Err(numerical("IRLS dispersion must be finite and positive"));
                    }
                }
            }
        }
        if let Some(e) = &self.guide_efficiencies {
            if e.len() != self.n_guides {
                return Err(numerical("IRLS efficiency shape mismatch"));
            }
            for &x in *e {
                if !x.is_finite() || !(0.0..=1.0).contains(&x) {
                    return Err(numerical("IRLS efficiency must be finite in [0, 1]"));
                }
            }
        }
        Ok((n, p))
    }
}

fn mat_vec(x: &[f64], rows: usize, cols: usize, v: &[f64]) -> Result<Vec<f64>, MageckError> {
    let mut out = vec![0.0; rows];
    for r in 0..rows {
        let mut s = 0.0;
        for c in 0..cols {
            s += x[r * cols + c] * v[c];
            finite(s, "matrix-vector intermediate")?;
        }
        out[r] = s;
    }
    Ok(out)
}
fn normal_system(
    x: &[f64],
    w: &[f64],
    z: &[f64],
    n: usize,
    p: usize,
    ridge: f64,
) -> Result<(Vec<f64>, Vec<f64>), MageckError> {
    let mut a = vec![0.0; p * p];
    let mut rhs = vec![0.0; p];
    for r in 0..n {
        for i in 0..p {
            // Each extended-design row has one nonzero guide coefficient.
            // Skip the other guide columns while keeping accumulation order
            // and all nonzero multiplication associations exactly unchanged.
            if x[r * p + i] == 0.0 {
                continue;
            }
            rhs[i] += x[r * p + i] * w[r] * z[r];
            finite(rhs[i], "normal RHS")?;
            for j in 0..p {
                a[i * p + j] += x[r * p + i] * w[r] * x[r * p + j];
                finite(a[i * p + j], "normal matrix")?;
            }
        }
    }
    for i in 0..p {
        a[i * p + i] += ridge;
    }
    Ok((a, rhs))
}
fn weighted_square_system(x: &[f64], w: &[f64], n: usize, p: usize) -> Vec<f64> {
    let mut q = vec![0.0; p * p];
    for r in 0..n {
        for i in 0..p {
            if x[r * p + i] == 0.0 {
                continue;
            }
            for j in 0..p {
                q[i * p + j] += x[r * p + i] * w[r] * w[r] * x[r * p + j];
            }
        }
    }
    q
}
fn mat_mul(a: &[f64], ar: usize, ac: usize, b: &[f64], bc: usize) -> Result<Vec<f64>, MageckError> {
    let mut out = vec![0.0; ar * bc];
    for i in 0..ar {
        for k in 0..ac {
            for j in 0..bc {
                out[i * bc + j] += a[i * ac + k] * b[k * bc + j];
                finite(out[i * bc + j], "matrix product")?;
            }
        }
    }
    Ok(out)
}
fn infinity_norm(a: &[f64], rows: usize, cols: usize) -> Result<f64, MageckError> {
    let mut out: f64 = 0.0;
    for r in 0..rows {
        let mut s: f64 = 0.0;
        for c in 0..cols {
            finite(a[r * cols + c], "matrix input")?;
            s += a[r * cols + c].abs();
            finite(s, "matrix norm")?;
        }
        out = out.max(s);
    }
    Ok(out)
}
fn invert(a: &[f64], dim: usize, iteration: usize) -> Result<Vec<f64>, MageckError> {
    if a.len() != dim * dim {
        return Err(numerical("IRLS inverse shape mismatch"));
    }
    let norm = infinity_norm(a, dim, dim)?;
    let threshold = 64.0 * f64::EPSILON * norm.max(1.0);
    let mut aug = vec![0.0; dim * 2 * dim];
    for r in 0..dim {
        for c in 0..dim {
            aug[r * 2 * dim + c] = a[r * dim + c];
        }
        aug[r * 2 * dim + dim + r] = 1.0;
    }
    for col in 0..dim {
        let mut best: Option<(usize, f64)> = None;
        for r in col..dim {
            let mut scale: f64 = 0.0;
            for c in col..dim {
                scale = scale.max(aug[r * 2 * dim + c].abs());
            }
            if !scale.is_finite() || scale == 0.0 {
                return Err(numerical(format!(
                    "IRLS singular normal system at iteration {iteration}, dimension {dim}"
                )));
            }
            let score = aug[r * 2 * dim + col].abs() / scale;
            finite(score, "pivot score")?;
            if best.map_or(true, |(b, bs): (usize, f64)| {
                score > bs || (score == bs && r < b)
            }) {
                best = Some((r, score));
            }
        }
        let pivot_row = best.unwrap().0;
        if pivot_row != col {
            for c in 0..2 * dim {
                aug.swap(col * 2 * dim + c, pivot_row * 2 * dim + c);
            }
        }
        let pivot = aug[col * 2 * dim + col];
        if !pivot.is_finite() || pivot.abs() <= threshold {
            return Err(numerical(format!(
                "IRLS singular normal system at iteration {iteration}, dimension {dim}"
            )));
        }
        for c in 0..2 * dim {
            aug[col * 2 * dim + c] /= pivot;
            finite(aug[col * 2 * dim + c], "inverse normalization")?;
        }
        for r in 0..dim {
            if r != col {
                let q = aug[r * 2 * dim + col];
                for c in 0..2 * dim {
                    aug[r * 2 * dim + c] -= q * aug[col * 2 * dim + c];
                    finite(aug[r * 2 * dim + c], "inverse elimination")?;
                }
            }
        }
    }
    let mut inv = vec![0.0; dim * dim];
    for r in 0..dim {
        for c in 0..dim {
            inv[r * dim + c] = aug[r * 2 * dim + dim + c];
        }
    }
    let product = mat_mul(a, dim, dim, &inv, dim)?;
    let mut residual = vec![0.0; dim * dim];
    for i in 0..dim {
        for j in 0..dim {
            residual[i * dim + j] = product[i * dim + j] - if i == j { 1.0 } else { 0.0 };
        }
    }
    let bound = 256.0 * f64::EPSILON * (norm * infinity_norm(&inv, dim, dim)?).max(1.0);
    if infinity_norm(&residual, dim, dim)? > bound {
        return Err(numerical(format!(
            "IRLS inverse residual failure at iteration {iteration}, dimension {dim}"
        )));
    }
    Ok(inv)
}

fn efficiency_multiplier(i: usize, n_samples: usize, n_guides: usize, efficiencies: &[f64]) -> f64 {
    if i < n_guides {
        1.0
    } else if i < n_samples * n_guides {
        efficiencies[i % n_guides]
    } else {
        1.0 - efficiencies[i % n_guides]
    }
}

fn treatment_posterior_minimum(
    null_terms: &[f64],
    terms: &[f64],
    n_samples: usize,
    n_guides: usize,
    guide: usize,
) -> Result<f64, MageckError> {
    let mut minimum = f64::INFINITY;
    for sample in 0..(n_samples - 1) {
        let index = (sample + 1) * n_guides + guide;
        let mut difference = null_terms[index] - terms[index];
        finite(difference, "efficiency likelihood difference")?;
        if difference > 100.0 {
            difference = 100.0;
        }
        minimum = minimum.min(difference);
    }
    Ok(minimum)
}

struct State {
    mu: Vec<f64>,
    residual: Vec<f64>,
    w: Vec<f64>,
    a: Vec<f64>,
    rhs: Vec<f64>,
    r: f64,
}
fn state(
    input: &IrlsInput<'_>,
    beta: &[f64],
    efficiencies: Option<&[f64]>,
    n: usize,
    p: usize,
) -> Result<State, MageckError> {
    let logmu = mat_vec(input.extended_design, n, p, beta)?;
    let mut mu = vec![0.0; n];
    let mut residual = vec![0.0; n];
    let mut z = vec![0.0; n];
    let mut w = vec![0.0; n];
    let mut score = vec![0.0; p];
    for i in 0..n {
        mu[i] = input.size_factors[i] * logmu[i].exp();
        if !mu[i].is_finite() || mu[i] <= 0.0 {
            return Err(numerical("IRLS invalid mean"));
        }
        residual[i] = input.counts[i] - mu[i];
        let e = efficiencies.map_or(1.0, |e| {
            efficiency_multiplier(i, input.n_samples, input.n_guides, e)
        });
        let working = residual[i] * e / mu[i];
        z[i] = logmu[i] + working;
        finite(z[i], "working response")?;
        let alpha = match input.dispersion {
            DispersionInput::Scalar(x) => x,
            DispersionInput::Vector(v) => v[i],
        };
        w[i] = 1.0 / (1.0 / mu[i] + alpha);
        if efficiencies.is_some() && i >= input.n_guides {
            w[i] *= e.max(0.01);
        }
        if !w[i].is_finite() || w[i] <= 0.0 {
            return Err(numerical("IRLS invalid weight"));
        }
    }
    let dispersion = match input.dispersion {
        DispersionInput::Scalar(ref x) => std::slice::from_ref(x),
        DispersionInput::Vector(v) => v,
    };
    for t in nb_log_likelihood_terms(input.counts, &mu, dispersion)? {
        finite(t, "likelihood")?;
    }
    let floor = w.iter().copied().fold(0.0, f64::max) / 100.0;
    for q in &mut w {
        *q = q.max(floor);
    }
    for i in 0..n {
        let e = efficiencies.map_or(1.0, |e| {
            efficiency_multiplier(i, input.n_samples, input.n_guides, e)
        });
        for j in 0..p {
            score[j] += input.extended_design[i * p + j] * w[i] * (residual[i] * e / mu[i]);
            finite(score[j], "score")?;
        }
    }
    let (a, rhs) = normal_system(input.extended_design, &w, &z, n, p, input.alpha_val)?;
    let mut r: f64 = 0.0;
    for j in 0..p {
        let denominator = 1.0
            + rhs[j].abs()
            + (0..p)
                .map(|l| a[j * p + l].abs() * beta[l].abs())
                .sum::<f64>();
        let q = (score[j] - input.alpha_val * beta[j]).abs() / denominator;
        finite(denominator, "residual denominator")?;
        finite(q, "scaled residual")?;
        r = r.max(q);
    }
    Ok(State {
        mu,
        residual,
        w,
        a,
        rhs,
        r,
    })
}
fn updated_efficiencies(
    input: &IrlsInput<'_>,
    beta: &[f64],
    n: usize,
    p: usize,
) -> Result<Vec<f64>, MageckError> {
    let mut null = beta.to_vec();
    null[input.n_guides..].fill(0.0);
    let means = |b: &[f64]| -> Result<Vec<f64>, MageckError> {
        mat_vec(input.extended_design, n, p, b)?
            .iter()
            .zip(input.size_factors)
            .map(|(l, s)| {
                let x = l.exp() * s;
                if x.is_finite() && x > 0.0 {
                    Ok(x)
                } else {
                    Err(numerical("IRLS invalid efficiency mean"))
                }
            })
            .collect()
    };
    let dispersion = match input.dispersion {
        DispersionInput::Scalar(ref x) => std::slice::from_ref(x),
        DispersionInput::Vector(v) => v,
    };
    let terms = nb_log_likelihood_terms(input.counts, &means(beta)?, dispersion)?;
    let null_terms = nb_log_likelihood_terms(input.counts, &means(&null)?, dispersion)?;
    (0..input.n_guides)
        .map(|g| {
            let q = 1.0
                / (1.0
                    + treatment_posterior_minimum(
                        &null_terms,
                        &terms,
                        input.n_samples,
                        input.n_guides,
                        g,
                    )?
                    .exp());
            finite(q, "efficiency")?;
            Ok(q)
        })
        .collect()
}
fn full_change(old: &[f64], new: &[f64]) -> f64 {
    old.iter()
        .zip(new)
        .map(|(a, b)| (b - a).abs() / (1e-8 + 1e-6 * a.abs().max(b.abs())))
        .fold(0.0, f64::max)
}
fn advance_streak(streak: usize, d: f64, r: f64, stable_eff: bool) -> usize {
    if d.is_finite() && r.is_finite() && d <= 1.0 && r <= 1e-9 && stable_eff {
        streak + 1
    } else {
        0
    }
}
fn relaxed_candidate(
    beta: &[f64],
    raw: &[f64],
    relaxation: f64,
) -> Result<Option<Vec<f64>>, MageckError> {
    if raw.iter().any(|x| !x.is_finite()) {
        return Err(numerical("IRLS non-finite raw direction"));
    }
    let candidate: Vec<_> = beta
        .iter()
        .zip(raw)
        .map(|(old, raw)| old + relaxation * (raw - old))
        .collect();
    if candidate.iter().any(|x| !x.is_finite()) {
        return Err(numerical("IRLS non-finite candidate"));
    }
    Ok(if candidate.iter().any(|x| x.abs() > MAX_ABS_BETA) {
        None
    } else {
        Some(candidate)
    })
}
fn classify(error: &MageckError) -> IrlsStatus {
    let s = error.to_string();
    if s.contains("singular") {
        IrlsStatus::SingularSystem
    } else if s.contains("inverse residual") {
        IrlsStatus::InverseResidualFailed
    } else {
        IrlsStatus::NonfiniteIterate
    }
}
pub(crate) fn em_whileloop(input: &IrlsInput<'_>) -> Result<IrlsResult, MageckError> {
    em_with_limit(input, MAX_ITERATIONS)
}
fn em_with_limit(input: &IrlsInput<'_>, limit: usize) -> Result<IrlsResult, MageckError> {
    let (n, p) = input.validate()?;
    let mut beta = input.beta_init.to_vec();
    let mut efficiencies = input.guide_efficiencies.map(ToOwned::to_owned);
    // An invalid initial state or initial solve cannot be rescued by relaxation.
    let initial = (|| {
        let current = state(input, &beta, efficiencies.as_deref(), n, p)?;
        invert(&current.a, p, 0)?;
        Ok(current)
    })();
    let mut current = match initial {
        Ok(current) => current,
        Err(error) => {
            return Ok(IrlsResult {
                status: classify(&error),
                iterations: 0,
                beta,
                efficiencies,
                condition_se: Vec::new(),
                means: Vec::new(),
                residuals: Vec::new(),
                diagnostics: UpdateDiagnostics {
                    failure: Some(format!("initial state/solve: {error}")),
                    ..Default::default()
                },
            })
        }
    };
    let mut diagnostics = UpdateDiagnostics::default();
    for iteration in 1..=limit {
        let update = (|| {
            let inv = invert(&current.a, p, iteration)?;
            let raw = mat_vec(&inv, p, p, &current.rhs)?;
            let Some(candidate) = relaxed_candidate(&beta, &raw, input.relaxation)? else {
                return Ok(None);
            };
            let next_eff = if efficiencies.is_some() && input.update_efficiency {
                Some(updated_efficiencies(input, &candidate, n, p)?)
            } else {
                efficiencies.clone()
            };
            let next = state(input, &candidate, next_eff.as_deref(), n, p)?;
            Ok(Some((candidate, next_eff, next)))
        })();
        let failure = match update {
            Ok(Some((candidate, next_eff, next))) => {
                let d = full_change(&beta, &candidate);
                let stable_eff = match (&efficiencies, &next_eff) {
                    (Some(a), Some(b)) => a.iter().zip(b).all(|(x, y)| (x - y).abs() <= 1e-8),
                    _ => true,
                };
                diagnostics.d = Some(d);
                diagnostics.r = Some(next.r);
                diagnostics.accepted += 1;
                diagnostics.streak = advance_streak(diagnostics.streak, d, next.r, stable_eff);
                beta = candidate;
                efficiencies = next_eff;
                current = next;
                if diagnostics.streak == 3 {
                    let inference = (|| {
                        let inv = invert(&current.a, p, iteration)?;
                        let q = weighted_square_system(input.extended_design, &current.w, n, p);
                        let cov = mat_mul(&mat_mul(&inv, p, p, &q, p)?, p, p, &inv, p)?;
                        let g = input.n_guides;
                        let b = input.n_conditions;
                        let mut se = vec![0.0; b * b];
                        for i in 0..b {
                            for j in 0..b {
                                se[i * b + j] = cov[(g + i) * p + g + j];
                            }
                        }
                        let se = final_scale(
                            input.extended_design,
                            n,
                            p,
                            g,
                            b,
                            &current.w,
                            &inv,
                            &current.mu,
                            &current.residual,
                            se,
                        )?;
                        for i in 0..b {
                            if !se[i * b + i].is_finite() || se[i * b + i] <= 0.0 {
                                return Err(numerical("IRLS invalid final covariance diagonal"));
                            }
                        }
                        Ok(se)
                    })();
                    let (status, condition_se) = match inference {
                        Ok(se) => (IrlsStatus::Converged, se),
                        Err(e) => {
                            diagnostics.failure = Some(format!("final inference: {e}"));
                            (IrlsStatus::InferenceInvalid, Vec::new())
                        }
                    };
                    return Ok(IrlsResult {
                        status,
                        iterations: iteration,
                        beta,
                        efficiencies,
                        condition_se,
                        means: current.mu,
                        residuals: current.residual,
                        diagnostics,
                    });
                }
                None
            }
            Ok(None) => Some((
                IrlsStatus::BetaRejected,
                "relaxed coefficient bound".to_owned(),
            )),
            Err(e) => Some((classify(&e), e.to_string())),
        };
        if let Some((status, reason)) = failure {
            diagnostics.failure = Some(format!("candidate {iteration}: {reason}"));
            return Ok(IrlsResult {
                status,
                iterations: iteration,
                beta,
                efficiencies,
                condition_se: Vec::new(),
                means: current.mu,
                residuals: current.residual,
                diagnostics,
            });
        }
    }
    diagnostics.failure = Some(format!("candidate limit {limit}"));
    Ok(IrlsResult {
        status: IrlsStatus::IterationLimitReached,
        iterations: limit,
        beta,
        efficiencies,
        condition_se: Vec::new(),
        means: current.mu,
        residuals: current.residual,
        diagnostics,
    })
}
fn final_scale(
    x: &[f64],
    n: usize,
    p: usize,
    _g: usize,
    b: usize,
    w: &[f64],
    inv: &[f64],
    mu: &[f64],
    res: &[f64],
    mut se: Vec<f64>,
) -> Result<Vec<f64>, MageckError> {
    let mut h = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for a in 0..p {
                for c in 0..p {
                    s += x[i * p + a] * inv[a * p + c] * x[j * p + c];
                }
            }
            h[i * n + j] = w[i].sqrt() * s * w[j].sqrt();
            finite(h[i * n + j], "hat matrix")?;
        }
    }
    let mut trace = 0.0;
    for i in 0..n {
        let mut hh = 0.0;
        for j in 0..n {
            hh += h[i * n + j] * h[i * n + j];
        }
        trace += 2.0 * h[i * n + i] - hh;
    }
    let dof = n as f64 - trace;
    if !dof.is_finite() || dof <= 0.0 {
        return Err(numerical(
            "IRLS degrees of freedom must be finite and positive",
        ));
    }
    let mut rss = 0.0;
    for i in 0..n {
        rss += (res[i] / mu[i]).powi(2);
    }
    let scale = rss / dof;
    finite(scale, "residual scale")?;
    for q in &mut se {
        *q *= scale;
        finite(*q, "standard error")?;
    }
    debug_assert_eq!(se.len(), b * b);
    Ok(se)
}

#[cfg(test)]
mod tests {
    use super::*;
    enum TestDispersion {
        Scalar(f64),
        Vector(Vec<f64>),
    }
    struct TestIrlsInput {
        n_samples: usize,
        n_guides: usize,
        n_conditions: usize,
        counts: Vec<f64>,
        extended_design: Vec<f64>,
        beta_init: Vec<f64>,
        size_factors: Vec<f64>,
        dispersion: TestDispersion,
        alpha_val: f64,
        guide_efficiencies: Option<Vec<f64>>,
        update_efficiency: bool,
    }
    impl TestIrlsInput {
        fn as_input(&self) -> IrlsInput<'_> {
            IrlsInput {
                relaxation: 1.0,
                n_samples: self.n_samples,
                n_guides: self.n_guides,
                n_conditions: self.n_conditions,
                counts: &self.counts,
                extended_design: &self.extended_design,
                beta_init: &self.beta_init,
                size_factors: &self.size_factors,
                dispersion: match &self.dispersion {
                    TestDispersion::Scalar(value) => DispersionInput::Scalar(*value),
                    TestDispersion::Vector(values) => DispersionInput::Vector(values),
                },
                alpha_val: self.alpha_val,
                guide_efficiencies: self.guide_efficiencies.as_deref(),
                update_efficiency: self.update_efficiency,
            }
        }
    }
    fn fixed() -> TestIrlsInput {
        TestIrlsInput {
            n_samples: 2,
            n_guides: 2,
            n_conditions: 1,
            counts: vec![101., 81., 61., 41., 61., 41.],
            extended_design: vec![
                1., 0., 0., 0., 1., 0., 1., 0., 1., 0., 1., 1., 1., 0., 0., 0., 1., 0.,
            ],
            beta_init: vec![101f64.ln(), 81f64.ln(), -0.5866949211069701],
            size_factors: vec![1.; 6],
            dispersion: TestDispersion::Scalar(0.1),
            alpha_val: 0.01,
            guide_efficiencies: None,
            update_efficiency: false,
        }
    }
    fn close(a: f64, b: f64) {
        assert!((a - b).abs() <= 2e-10 + 2e-10 * b.abs(), "{a} {b}");
    }
    fn assert_close_slice(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (&a, &b) in actual.iter().zip(expected) {
            close(a, b);
        }
    }
    fn meanvar() -> TestIrlsInput {
        let mut extended_design = Vec::new();
        for condition in [[0., 0.], [1., 0.], [0., 1.], [0., 0.], [0., 0.]] {
            for guide in 0..3 {
                for col in 0..3 {
                    extended_design.push(if guide == col { 1. } else { 0. });
                }
                extended_design.extend(condition);
            }
        }
        TestIrlsInput {
            n_samples: 3,
            n_guides: 3,
            n_conditions: 2,
            counts: vec![
                121., 91., 61., 81., 51., 41., 51., 36., 31., 81., 51., 41., 51., 36., 31.,
            ],
            extended_design,
            beta_init: vec![
                4.795790545596741,
                4.51085950651685,
                4.110873864173311,
                -0.5490454130689122,
                -0.7102675226560533,
            ],
            size_factors: vec![
                1., 1., 1., 1.1, 1.1, 1.1, 0.9, 0.9, 0.9, 1.1, 1.1, 1.1, 0.9, 0.9, 0.9,
            ],
            dispersion: TestDispersion::Vector(vec![0.01; 15]),
            alpha_val: 0.01,
            guide_efficiencies: None,
            update_efficiency: false,
        }
    }

    #[test]
    fn sparse_systems_are_bit_identical_to_dense() {
        let (n, p) = (4, 6);
        let x = [
            1., 0., 0., 0., -0.5, 1., 0., 1., 0., 0., 1., -0.25, 0., 0., 1., 0., 0., 0.75, 0., 0.,
            0., 1., -1., 0.,
        ];
        let w = [0.125, 0.7, 1.25, 2.3];
        let z = [-3.5, 0.25, 2.7, -1.2];
        let ridge = 0.01;
        let mut dense_a = vec![0.0; p * p];
        let mut dense_rhs = vec![0.0; p];
        let mut dense_q = vec![0.0; p * p];
        // Original dense arithmetic, including every zero multiplication.
        for r in 0..n {
            for i in 0..p {
                dense_rhs[i] += x[r * p + i] * w[r] * z[r];
                for j in 0..p {
                    dense_a[i * p + j] += x[r * p + i] * w[r] * x[r * p + j];
                    dense_q[i * p + j] += x[r * p + i] * w[r] * w[r] * x[r * p + j];
                }
            }
        }
        for i in 0..p {
            dense_a[i * p + i] += ridge;
        }
        let (a, rhs) = normal_system(&x, &w, &z, n, p, ridge).unwrap();
        let q = weighted_square_system(&x, &w, n, p);
        for (sparse, dense) in [(&a, &dense_a), (&rhs, &dense_rhs), (&q, &dense_q)] {
            assert_eq!(
                sparse.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                dense.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
            );
        }
        let sparse_inverse = invert(&a, p, 2).unwrap();
        let dense_inverse = invert(&dense_a, p, 2).unwrap();
        let sparse_beta = mat_vec(&sparse_inverse, p, p, &rhs).unwrap();
        let dense_beta = mat_vec(&dense_inverse, p, p, &dense_rhs).unwrap();
        assert_eq!(
            sparse_beta.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            dense_beta.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
    }
    #[test]
    fn fixed_extended_design_has_legacy_duplicate_baseline_rows() {
        assert_eq!(
            fixed().extended_design,
            vec![1., 0., 0., 0., 1., 0., 1., 0., 1., 0., 1., 1., 1., 0., 0., 0., 1., 0.,]
        );
    }

    #[test]
    fn efficiency_partition_ends_treatment_at_samples_times_guides() {
        let efficiencies = [0.8, 0.3];
        let multipliers: Vec<_> = (0..10)
            .map(|i| efficiency_multiplier(i, 3, 2, &efficiencies))
            .collect();
        assert_close_slice(
            &multipliers,
            &[1., 1., 0.8, 0.3, 0.8, 0.3, 0.2, 0.7, 0.2, 0.7],
        );
    }

    #[test]
    fn posterior_minimum_uses_only_original_treatment_rows() {
        let null_terms = vec![0.; 10];
        let mut terms = vec![0.; 10];
        terms[2] = -2.;
        terms[4] = -3.;
        terms[6] = -99.;
        assert_eq!(
            treatment_posterior_minimum(&null_terms, &terms, 3, 2, 0).unwrap(),
            2.
        );
        terms[6] = -1_000.;
        assert_eq!(
            treatment_posterior_minimum(&null_terms, &terms, 3, 2, 0).unwrap(),
            2.
        );
    }
    #[test]
    fn fixed_alpha_direct_oracle() {
        let fixture = fixed();
        let r = em_whileloop(&fixture.as_input()).unwrap();
        assert_eq!(r.status, IrlsStatus::Converged);
        assert_eq!(r.iterations, 7);
        assert_close_slice(
            &r.beta,
            &[4.409979456412462, 4.089897634487935, -0.3357395722824612],
        );
        assert_close_slice(
            &r.means,
            &[
                82.26777341163684,
                59.73377671252765,
                58.805764400989986,
                42.69825539779795,
                82.26777341163684,
                59.73377671252765,
            ],
        );
        assert_close_slice(
            &r.residuals,
            &[
                18.73222658836316,
                21.26622328747235,
                2.194235599010014,
                -1.6982553977979506,
                -21.26777341163684,
                -18.73377671252765,
            ],
        );
        assert_close_slice(&r.condition_se, &[0.0865552907745817]);
        assert_eq!(r.diagnostics.streak, 3);
        assert!(r.diagnostics.r.unwrap() <= 1e-9);
    }
    #[test]
    fn efficiency_direct_oracle() {
        let mut x = fixed();
        x.counts = vec![101., 81., 31., 76., 31., 76.];
        x.beta_init[2] = -0.6162619439317928;
        x.guide_efficiencies = Some(vec![1., 0.4]);
        x.update_efficiency = true;
        let r = em_whileloop(&x.as_input()).unwrap();
        assert_eq!(r.status, IrlsStatus::Converged);
        assert_eq!(r.iterations, 10);
        assert_close_slice(
            &r.beta,
            &[4.608616250418932, 4.360662134625523, -1.173016410220073],
        );
        assert_close_slice(
            &r.means,
            &[
                100.34520089372944,
                78.30896834010183,
                31.050032038008013,
                24.231313049027545,
                100.34520089372944,
                78.30896834010183,
            ],
        );
        assert_close_slice(
            &r.residuals,
            &[
                0.6547991062705592,
                2.6910316598981723,
                -0.05003203800801259,
                51.76868695097245,
                -69.34520089372944,
                -2.3089683401018277,
            ],
        );
        assert_close_slice(&r.condition_se, &[3.1554962758588356]);
        assert_eq!(r.diagnostics.streak, 3);
        assert!(r.diagnostics.r.unwrap() <= 1e-9);
        assert_close_slice(
            r.efficiencies.as_deref().unwrap(),
            &[0.984905188862146, 0.0004772592168971316],
        );
    }
    #[test]
    fn meanvar_alpha_direct_oracle() {
        let fixture = meanvar();
        let r = em_whileloop(&fixture.as_input()).unwrap();
        assert_eq!(r.status, IrlsStatus::Converged);
        assert_eq!(r.iterations, 7);
        assert_close_slice(
            &r.beta,
            &[
                4.435143962310108,
                4.062897588884515,
                3.8087753670906204,
                -0.17924772617587384,
                -0.35480592369873265,
            ],
        );
        assert_close_slice(
            &r.means,
            &[
                84.36426928919911,
                58.142540410344694,
                45.095179995831856,
                77.57199058311384,
                53.46140771667902,
                41.46450751548615,
                53.24892001383578,
                36.69832631511967,
                28.46311184627309,
                92.80069621811903,
                63.95679445137917,
                49.604697995415044,
                75.9278423602792,
                52.32828636931023,
                40.585661996248675,
            ],
        );
        assert_close_slice(
            &r.residuals,
            &[
                36.635730710800885,
                32.857459589655306,
                15.904820004168144,
                3.4280094168861552,
                -2.4614077166790196,
                -0.4645075154861473,
                -2.2489200138357788,
                -0.6983263151196724,
                2.5368881537269097,
                -11.800696218119029,
                -12.956794451379167,
                -8.604697995415044,
                -24.9278423602792,
                -16.328286369310227,
                -9.585661996248675,
            ],
        );
        assert_close_slice(
            &r.condition_se,
            &[
                0.0454111651245559,
                0.011357547802493402,
                0.011357547802493402,
                0.045670173251520334,
            ],
        );
        assert_eq!(r.diagnostics.streak, 3);
        assert!(r.diagnostics.r.unwrap() <= 1e-9);
    }
    #[test]
    fn validation_and_inverse_contract() {
        let mut x = fixed();
        x.counts[0] = 1.5;
        assert!(em_whileloop(&x.as_input()).is_err());
        assert!(invert(&[1., 1., 1., 1.], 2, 1).is_err());
        let q = invert(&[2., 0., 0., 3.], 2, 1).unwrap();
        assert_eq!(q, [0.5, 0., 0., 1. / 3.]);
    }
    #[test]
    fn deterministic() {
        let fixture_a = fixed();
        let fixture_b = fixed();
        let a = em_whileloop(&fixture_a.as_input()).unwrap();
        let b = em_whileloop(&fixture_b.as_input()).unwrap();
        assert_eq!(
            a.beta.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            b.beta.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
    }
    #[test]
    fn full_vector_residual_streak_and_finite_bound_contract() {
        // A zero condition change cannot conceal a moving guide intercept.
        assert!(full_change(&[5., 0.], &[5.1, 0.]) > 1.);
        assert_eq!(advance_streak(2, 0.5, 1e-9, true), 3);
        for (d, r, e) in [
            (1.01, 0., true),
            (0., 1.01e-9, true),
            (0., 0., false),
            (f64::NAN, 0., true),
            (0., f64::INFINITY, true),
        ] {
            assert_eq!(advance_streak(2, d, r, e), 0);
        }
        assert_eq!(
            relaxed_candidate(&[0.], &[100.], 1.).unwrap(),
            Some(vec![100.])
        );
        assert!(relaxed_candidate(&[0.], &[100.000001], 1.)
            .unwrap()
            .is_none());
        assert_eq!(
            relaxed_candidate(&[0.], &[160.], 0.5).unwrap(),
            Some(vec![80.])
        );
        assert!(relaxed_candidate(&[0.], &[f64::NAN], 0.25).is_err());
        let mut fixture = fixed();
        fixture.beta_init[0] = 101.;
        assert!(em_whileloop(&fixture.as_input()).is_err());
        fixture = fixed();
        fixture.beta_init[0] = 99.;
        fixture.size_factors[0] = 1e300;
        let failed = em_whileloop(&fixture.as_input()).unwrap();
        assert_eq!(failed.status, IrlsStatus::NonfiniteIterate);
        assert_eq!(failed.iterations, 0);
        assert!(failed
            .diagnostics
            .failure
            .unwrap()
            .contains("initial state/solve"));
    }
    #[test]
    fn update_limit_accepts_a_completed_streak_on_last_request() {
        let mut fixture = fixed();
        let reference = em_whileloop(&fixture.as_input()).unwrap();
        fixture.beta_init = reference.beta;
        let failed = em_with_limit(&fixture.as_input(), 2).unwrap();
        assert_eq!(failed.status, IrlsStatus::IterationLimitReached);
        assert_eq!(failed.iterations, 2);
        assert_eq!(failed.diagnostics.streak, 2);
        assert!(failed.condition_se.is_empty());
        let success = em_with_limit(&fixture.as_input(), 3).unwrap();
        assert_eq!(success.status, IrlsStatus::Converged);
        assert_eq!(success.iterations, 3);
        // Explicit 1,000-request counter boundary using the same pass transition.
        let mut streak = 0;
        for request in 1..=1000 {
            streak = advance_streak(streak, if request < 998 { 2. } else { 0. }, 0., true);
            if request < 1000 {
                assert!(streak < 3);
            }
        }
        assert_eq!(streak, 3);
    }
}
