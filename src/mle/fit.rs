//! Private, single-gene preparation and retry control for the classical MLE fit.

use super::irls::{em_whileloop, DispersionInput, IrlsInput, IrlsStatus};
use crate::error::MageckError;

const RELAXATIONS: [f64; 3] = [1.0, 0.5, 0.25];

#[derive(Clone, Debug)]
pub(super) struct FitInput {
    pub(super) counts: Vec<f64>,
    pub(super) n_samples: usize,
    pub(super) n_guides: usize,
    pub(super) design: Vec<f64>,
    pub(super) n_conditions: usize,
    pub(super) size_factors: Vec<f64>,
    pub(super) dispersion: FitDispersion,
    pub(super) alpha_val: f64,
    pub(super) estimate_efficiency: bool,
    pub(super) restart_beta: Option<Vec<f64>>,
    pub(super) restart_efficiencies: Option<Vec<f64>>,
    pub(super) update_efficiency: bool,
}

#[derive(Clone, Debug)]
pub(super) enum FitDispersion {
    Scalar(f64),
    Vector(Vec<f64>),
}

impl FitDispersion {
    fn as_irls_input(&self) -> DispersionInput<'_> {
        match self {
            Self::Scalar(value) => DispersionInput::Scalar(*value),
            Self::Vector(values) => DispersionInput::Vector(values),
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct FitResult {
    pub(super) history: Vec<FitAttempt>,
    pub(super) status: IrlsStatus,
    pub(super) attempts: usize,
    iterations: usize,
    pub(super) beta: Vec<f64>,
    pub(super) efficiencies: Option<Vec<f64>>,
    pub(super) condition_se: Vec<f64>,
    pub(super) means: Vec<f64>,
    pub(super) residuals: Vec<f64>,
    counts: Vec<f64>,
    size_factors: Vec<f64>,
    extended_design: Vec<f64>,
    extended_rows: usize,
    extended_cols: usize,
}

#[derive(Clone, Debug)]
pub(super) struct FitAttempt {
    pub(super) relaxation: f64,
    pub(super) alpha: Vec<f64>,
    pub(super) lambda: f64,
    pub(super) initial_beta: Vec<f64>,
    pub(super) initial_efficiencies: Option<Vec<f64>>,
    pub(super) result: super::irls::IrlsResult,
}
#[cfg(test)]
impl FitResult {
    pub(super) fn gene_result_fixture(
        status: IrlsStatus,
        attempts: usize,
        beta: Vec<f64>,
        efficiencies: Option<Vec<f64>>,
        condition_se: Vec<f64>,
        means: Vec<f64>,
        residuals: Vec<f64>,
    ) -> Self {
        Self {
            history: Vec::new(),
            status,
            attempts,
            iterations: 0,
            beta,
            efficiencies,
            condition_se,
            means,
            residuals,
            counts: Vec::new(),
            size_factors: Vec::new(),
            extended_design: Vec::new(),
            extended_rows: 0,
            extended_cols: 0,
        }
    }
}

fn numerical(message: impl Into<String>) -> MageckError {
    MageckError::Math(message.into())
}

fn product(a: usize, b: usize, label: &str) -> Result<usize, MageckError> {
    a.checked_mul(b)
        .ok_or_else(|| numerical(format!("fit {label} dimension overflow")))
}

fn finite(v: f64, label: &str) -> Result<(), MageckError> {
    if v.is_finite() {
        Ok(())
    } else {
        Err(numerical(format!("fit non-finite {label}")))
    }
}

// This is deliberately task-local: the ratified IRLS surface does not expose its
// inversion helper.  Its pivoting and residual policy matches the IRLS helper.
fn invert(a: &[f64], dim: usize) -> Result<Vec<f64>, MageckError> {
    if a.len() != product(dim, dim, "normal")? {
        return Err(numerical("fit inverse shape mismatch"));
    }
    let norm = a
        .chunks(dim)
        .map(|r| r.iter().map(|x| x.abs()).sum::<f64>())
        .fold(0.0, f64::max);
    if !norm.is_finite() {
        return Err(numerical("fit non-finite normal matrix"));
    }
    let threshold = 64.0 * f64::EPSILON * norm.max(1.0);
    let width = dim
        .checked_mul(2)
        .ok_or_else(|| numerical("fit inverse dimension overflow"))?;
    let mut aug = vec![0.0; product(dim, width, "inverse")?];
    for row in 0..dim {
        for col in 0..dim {
            aug[row * width + col] = a[row * dim + col];
        }
        aug[row * width + dim + row] = 1.0;
    }
    for col in 0..dim {
        let mut best: Option<(usize, f64)> = None;
        for row in col..dim {
            let scale = (col..dim)
                .map(|c| aug[row * width + c].abs())
                .fold(0.0, f64::max);
            if !scale.is_finite() || scale == 0.0 {
                return Err(numerical("fit singular condition normal system"));
            }
            let score = aug[row * width + col].abs() / scale;
            if best.map_or(true, |(old, value)| {
                score > value || (score == value && row < old)
            }) {
                best = Some((row, score));
            }
        }
        let row = best.expect("nonempty pivot range").0;
        if row != col {
            for c in 0..width {
                aug.swap(col * width + c, row * width + c);
            }
        }
        let pivot = aug[col * width + col];
        if !pivot.is_finite() || pivot.abs() <= threshold {
            return Err(numerical("fit singular condition normal system"));
        }
        for c in 0..width {
            aug[col * width + c] /= pivot;
            finite(aug[col * width + c], "inverse")?;
        }
        for row in 0..dim {
            if row == col {
                continue;
            }
            let q = aug[row * width + col];
            for c in 0..width {
                aug[row * width + c] -= q * aug[col * width + c];
                finite(aug[row * width + c], "inverse")?;
            }
        }
    }
    let mut inverse = vec![0.0; product(dim, dim, "inverse")?];
    for row in 0..dim {
        for col in 0..dim {
            inverse[row * dim + col] = aug[row * width + dim + col];
        }
    }
    Ok(inverse)
}

fn mat_vec(a: &[f64], rows: usize, cols: usize, v: &[f64]) -> Result<Vec<f64>, MageckError> {
    if a.len() != product(rows, cols, "matrix")? || v.len() != cols {
        return Err(numerical("fit matrix shape mismatch"));
    }
    (0..rows)
        .map(|r| {
            let x = (0..cols).map(|c| a[r * cols + c] * v[c]).sum::<f64>();
            finite(x, "matrix product")?;
            Ok(x)
        })
        .collect()
}

pub(super) fn setup_and_fit(input: &FitInput) -> Result<FitResult, MageckError> {
    setup_and_fit_with_runner(input, em_whileloop)
}

fn setup_and_fit_with_runner<F>(input: &FitInput, mut runner: F) -> Result<FitResult, MageckError>
where
    F: FnMut(&IrlsInput<'_>) -> Result<super::irls::IrlsResult, MageckError>,
{
    let s = input.n_samples;
    let g = input.n_guides;
    let b = input.n_conditions;
    if s < 2 || g == 0 || b == 0 {
        return Err(numerical(
            "fit dimensions must include at least two samples and one guide/condition",
        ));
    }
    if input.counts.len() != product(s, g, "count")?
        || input.design.len() != product(s, b + 1, "design")?
        || input.size_factors.len() != s
    {
        return Err(numerical("fit input shape mismatch"));
    }
    if !input.alpha_val.is_finite() || input.alpha_val < 0.0 {
        return Err(numerical("fit ridge alpha must be finite and non-negative"));
    }
    for &x in &input.counts {
        if !x.is_finite() || x <= 0.0 || x.fract() != 0.0 {
            return Err(numerical("fit counts must be positive finite integers"));
        }
    }
    for &x in &input.size_factors {
        if !x.is_finite() || x <= 0.0 {
            return Err(numerical("fit size factors must be finite and positive"));
        }
    }
    for &x in &input.design {
        finite(x, "design")?;
    }
    if (0..s).any(|r| input.design[r * (b + 1)] != 1.0)
        || input.design[0..b + 1]
            != [1.0]
                .iter()
                .chain(std::iter::repeat(&0.0))
                .take(b + 1)
                .copied()
                .collect::<Vec<_>>()
    {
        return Err(numerical("fit design requires literal row zero baseline"));
    }
    let baseline_rows: Vec<usize> = (0..s)
        .filter(|&r| {
            input.design[r * (b + 1) + 1..(r + 1) * (b + 1)]
                .iter()
                .all(|&x| x == 0.0)
        })
        .collect();
    if baseline_rows.is_empty() {
        return Err(numerical("fit design has no baseline-only row"));
    }
    let n = product(
        s.checked_mul(2)
            .and_then(|x| x.checked_sub(1))
            .ok_or_else(|| numerical("fit observation dimension overflow"))?,
        g,
        "observation",
    )?;
    match &input.dispersion {
        FitDispersion::Scalar(x) if x.is_finite() && *x > 0.0 => {}
        FitDispersion::Vector(xs)
            if xs.len() == n && xs.iter().all(|x| x.is_finite() && *x > 0.0) => {}
        _ => {
            return Err(numerical(
                "fit dispersion must be positive and observation-aligned",
            ))
        }
    }
    let p = g
        .checked_add(b)
        .ok_or_else(|| numerical("fit coefficient dimension overflow"))?;
    if let Some(beta) = &input.restart_beta {
        if beta.len() != p || beta.iter().any(|x| !x.is_finite()) {
            return Err(numerical("fit invalid restart beta"));
        }
    }
    if let Some(e) = &input.restart_efficiencies {
        if e.len() != g || e.iter().any(|x| !x.is_finite() || !(0.0..=1.0).contains(x)) {
            return Err(numerical("fit invalid restart efficiency"));
        }
    }

    let mut counts = input.counts.clone();
    counts.extend_from_slice(&input.counts[g..]);
    let mut sizes = Vec::with_capacity(n);
    for &x in &input.size_factors {
        sizes.extend(std::iter::repeat(x).take(g));
    }
    for &x in &input.size_factors[1..] {
        sizes.extend(std::iter::repeat(x).take(g));
    }
    let reduced: Vec<f64> = (1..s)
        .flat_map(|r| {
            input.design[r * (b + 1) + 1..(r + 1) * (b + 1)]
                .iter()
                .copied()
        })
        .collect();
    let mut extended = Vec::with_capacity(product(n, p, "extended design")?);
    for sample in 0..(2 * s - 1) {
        for guide in 0..g {
            for coefficient in 0..p {
                extended.push(if coefficient < g {
                    (coefficient == guide) as u8 as f64
                } else if sample < s {
                    if sample == 0 {
                        0.0
                    } else {
                        reduced[(sample - 1) * b + coefficient - g]
                    }
                } else {
                    0.0
                });
            }
        }
    }
    let guide_beta: Vec<f64> = (0..g)
        .map(|guide| {
            baseline_rows
                .iter()
                .map(|&row| input.counts[row * g + guide].ln() - input.size_factors[row].ln())
                .sum::<f64>()
                / baseline_rows.len() as f64
        })
        .collect();
    let treatment_means: Vec<f64> = (1..s)
        .map(|row| {
            (0..g)
                .map(|guide| {
                    input.counts[row * g + guide].ln()
                        - input.size_factors[row].ln()
                        - guide_beta[guide]
                })
                .sum::<f64>()
                / g as f64
        })
        .collect();
    let mut normal = vec![0.0; product(b, b, "condition normal")?];
    let mut rhs = vec![0.0; b];
    for row in 0..s - 1 {
        for i in 0..b {
            rhs[i] += reduced[row * b + i] * treatment_means[row];
            for j in 0..b {
                normal[i * b + j] += reduced[row * b + i] * reduced[row * b + j];
            }
        }
    }
    for i in 0..b {
        normal[i * b + i] += input.alpha_val;
    }
    let inverse = invert(&normal, b)?;
    let condition_beta = mat_vec(&inverse, b, b, &rhs)?;
    let mut initial_beta = guide_beta;
    initial_beta.extend(condition_beta);
    if let Some(beta) = &input.restart_beta {
        initial_beta = beta.clone();
    }
    if initial_beta.iter().any(|x| !x.is_finite()) {
        return Err(numerical("fit non-finite resulting initialization"));
    }
    let efficiencies = if input.estimate_efficiency {
        Some(
            input
                .restart_efficiencies
                .clone()
                .unwrap_or_else(|| vec![1.0; g]),
        )
    } else {
        None
    };
    let mut history = Vec::new();
    for (index, relaxation) in RELAXATIONS.into_iter().enumerate() {
        let mut r = runner(&IrlsInput {
            relaxation,
            n_samples: s,
            n_guides: g,
            n_conditions: b,
            counts: &counts,
            extended_design: &extended,
            beta_init: &initial_beta,
            size_factors: &sizes,
            dispersion: input.dispersion.as_irls_input(),
            alpha_val: input.alpha_val,
            guide_efficiencies: efficiencies.as_deref(),
            update_efficiency: input.update_efficiency,
        })?;
        history.push(FitAttempt {
            relaxation,
            alpha: match &input.dispersion {
                FitDispersion::Scalar(x) => vec![*x],
                FitDispersion::Vector(v) => v.clone(),
            },
            lambda: input.alpha_val,
            initial_beta: initial_beta.clone(),
            initial_efficiencies: efficiencies.clone(),
            result: r.clone(),
        });
        let retryable = matches!(
            r.status,
            IrlsStatus::BetaRejected | IrlsStatus::IterationLimitReached
        ) || (r.diagnostics.accepted > 0
            && matches!(
                r.status,
                IrlsStatus::NonfiniteIterate
                    | IrlsStatus::SingularSystem
                    | IrlsStatus::InverseResidualFailed
            ));
        if r.status.converged() || !retryable || index == RELAXATIONS.len() - 1 {
            if r.status.converged() && index > 0 {
                r.status = IrlsStatus::ConvergedSameModelRetry;
            }
            return Ok(FitResult {
                history,
                status: r.status,
                attempts: index + 1,
                iterations: r.iterations,
                beta: r.beta,
                efficiencies: r.efficiencies,
                condition_se: r.condition_se,
                means: r.means,
                residuals: r.residuals,
                counts,
                size_factors: sizes,
                extended_design: extended,
                extended_rows: n,
                extended_cols: p,
            });
        }
    }
    unreachable!("bounded attempts always return")
}

#[cfg(test)]
mod tests {
    use super::super::irls::IrlsResult;
    use super::*;
    fn input() -> FitInput {
        FitInput {
            counts: vec![100., 80., 60., 40., 90., 70., 45., 35., 110., 95., 70., 50.],
            n_samples: 3,
            n_guides: 4,
            design: vec![1., 0., 1., 1., 1., 1.],
            n_conditions: 1,
            size_factors: vec![1., 1., 1.],
            dispersion: FitDispersion::Scalar(0.01),
            alpha_val: 0.01,
            estimate_efficiency: false,
            restart_beta: None,
            restart_efficiencies: None,
            update_efficiency: true,
        }
    }
    #[test]
    fn standard_setup_and_fit() {
        let r = setup_and_fit(&input()).unwrap();
        assert_eq!(r.counts.len(), 20);
        assert_eq!(r.size_factors.len(), 20);
        assert_eq!((r.extended_rows, r.extended_cols), (20, 5));
        assert_eq!(r.efficiencies, None);
        assert!(r.beta.iter().all(|x| x.is_finite()));
        let expected = [
            4.603_564_991_104_803,
            4.405_098_595_757_647,
            4.058_804_911_650_852,
            3.735_990_305_704_632,
            0.003_552_630_807_881_318_5,
        ];
        for (actual, expected) in r.beta.iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-10, "{actual} != {expected}");
        }
    }
    #[test]
    fn strict_validation() {
        let mut x = input();
        x.counts[0] = 0.;
        assert!(setup_and_fit(&x).is_err());
        x = input();
        x.design[1] = 1.;
        assert!(setup_and_fit(&x).is_err());
        x = input();
        x.restart_efficiencies = Some(vec![1.1; 4]);
        assert!(setup_and_fit(&x).is_err());
    }

    fn assert_beta(input: FitInput, expected: &[f64]) {
        let actual = setup_and_fit(&input).unwrap().beta;
        for (&actual, &expected) in actual.iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-10, "{actual} != {expected}");
        }
    }

    #[test]
    fn direct_nonunit_and_multiple_baseline_oracles() {
        let mut nonunit = input();
        nonunit.counts = vec![
            120., 80., 40., 20., 250., 130., 70., 30., 70., 50., 28., 18.,
        ];
        nonunit.size_factors = vec![1., 2., 0.5];
        nonunit.alpha_val = 0.025;
        nonunit.dispersion = FitDispersion::Scalar(0.025);
        assert_beta(
            nonunit,
            &[
                4.852129428842234,
                4.359058298974036,
                3.7219857744037252,
                3.02987743158994,
                0.01596730432252102,
            ],
        );
        let multiple = FitInput {
            counts: vec![
                100., 64., 36., 25., 121., 81., 49., 36., 80., 48., 30., 20., 50., 35., 20., 12.,
            ],
            n_samples: 4,
            n_guides: 4,
            design: vec![1., 0., 0., 1., 0., 0., 1., 1., 0., 1., 0., 1.],
            n_conditions: 2,
            size_factors: vec![1., 1.1, 0.8, 1.3],
            dispersion: FitDispersion::Scalar(0.02),
            alpha_val: 0.02,
            estimate_efficiency: false,
            restart_beta: None,
            restart_efficiencies: None,
            update_efficiency: true,
        };
        assert_beta(
            multiple,
            &[
                4.51344379378081,
                4.080574212603518,
                3.562237939485547,
                3.184129421943147,
                0.054617783229513585,
                -0.8480907455660498,
            ],
        );
    }

    #[test]
    fn direct_restart_and_vector_alpha_oracles() {
        let restart = FitInput {
            counts: vec![100., 80., 60., 40., 95., 70., 55., 35., 90., 65., 50., 30.],
            restart_beta: Some(vec![4., 3.8, 3.5, 3.2, -0.12]),
            restart_efficiencies: Some(vec![0.2, 0.4, 0.6, 0.8]),
            estimate_efficiency: true,
            ..input()
        };
        assert_beta(
            restart,
            &[
                4.598552815099011,
                4.335660823574315,
                4.063404889584742,
                3.6243797575238084,
                -0.10693763026984876,
            ],
        );
        let mut vector = input();
        vector.dispersion = FitDispersion::Vector(vec![0.330718636511549; 20]);
        assert_beta(
            vector,
            &[
                4.598583574879413,
                4.400227615794217,
                4.054195629099106,
                3.7315300566191776,
                0.008665131173940615,
            ],
        );
    }

    #[test]
    fn retry_sequence_resets_beta_and_efficiency() {
        let mut fixture = input();
        fixture.estimate_efficiency = true;
        fixture.dispersion = FitDispersion::Vector(vec![0.07; 20]);
        let mut calls = 0;
        let mut immutable_pointers = None;
        let mut dispersion_pointer = None;
        let result = setup_and_fit_with_runner(&fixture, |call| {
            calls += 1;
            assert_eq!(call.alpha_val, 0.01);
            assert_eq!(call.relaxation, [1., 0.5, 0.25][calls as usize - 1]);
            assert_eq!(call.beta_init[0], 100_f64.ln());
            let pointers = (
                call.counts.as_ptr(),
                call.extended_design.as_ptr(),
                call.beta_init.as_ptr(),
                call.size_factors.as_ptr(),
            );
            if let Some(expected) = immutable_pointers {
                assert_eq!(pointers, expected);
            } else {
                immutable_pointers = Some(pointers);
            }
            match call.dispersion {
                DispersionInput::Vector(values) => {
                    assert!(values.iter().all(|&value| value == 0.07));
                    if let Some(expected) = dispersion_pointer {
                        assert_eq!(values.as_ptr(), expected);
                    } else {
                        dispersion_pointer = Some(values.as_ptr());
                    }
                }
                DispersionInput::Scalar(_) => panic!("vector dispersion fixture changed shape"),
            }
            let expected_efficiency = 1.0;
            assert_eq!(
                call.guide_efficiencies,
                Some([expected_efficiency; 4].as_slice())
            );
            let efficiencies = Some(vec![calls as f64 / 10.0; 4]);
            Ok(IrlsResult {
                diagnostics: Default::default(),
                status: if calls == 3 {
                    IrlsStatus::Converged
                } else {
                    IrlsStatus::BetaRejected
                },
                iterations: calls as usize,
                beta: call.beta_init.to_vec(),
                efficiencies,
                condition_se: vec![0.0],
                means: vec![1.0; 20],
                residuals: vec![0.0; 20],
            })
        })
        .unwrap();
        assert_eq!(calls, 3);
        assert_eq!(result.attempts, 3);
        assert_eq!(result.status, IrlsStatus::ConvergedSameModelRetry);
        assert_eq!(result.efficiencies, Some(vec![0.3; 4]));
    }

    #[test]
    fn retry_exhaustion_is_ineligible_after_three_attempts() {
        let mut calls = 0;
        let result = setup_and_fit_with_runner(&input(), |call| {
            calls += 1;
            Ok(IrlsResult {
                diagnostics: Default::default(),
                status: IrlsStatus::BetaRejected,
                iterations: calls,
                beta: vec![calls as f64; call.beta_init.len()],
                efficiencies: None,
                condition_se: vec![0.0],
                means: vec![1.0; 20],
                residuals: vec![0.0; 20],
            })
        })
        .unwrap();
        assert_eq!(calls, 3);
        assert_eq!(result.attempts, 3);
        assert!(!result.status.converged());
        assert_eq!(result.beta[0], 3.0);
    }
    #[test]
    fn retry_failure_classes_and_initial_error_do_not_change_model() {
        for (status, accepted, expected_calls) in [
            (IrlsStatus::InferenceInvalid, 5, 1),
            (IrlsStatus::NonfiniteIterate, 0, 1),
            (IrlsStatus::SingularSystem, 0, 1),
            (IrlsStatus::NonfiniteIterate, 1, 3),
            (IrlsStatus::SingularSystem, 1, 3),
            (IrlsStatus::InverseResidualFailed, 1, 3),
            (IrlsStatus::IterationLimitReached, 1000, 3),
        ] {
            let mut calls = 0;
            let result = setup_and_fit_with_runner(&input(), |call| {
                calls += 1;
                assert_eq!(call.alpha_val, 0.01);
                assert!(matches!(call.dispersion, DispersionInput::Scalar(0.01)));
                Ok(IrlsResult {
                    status,
                    iterations: accepted,
                    diagnostics: super::super::irls::UpdateDiagnostics {
                        accepted,
                        ..Default::default()
                    },
                    beta: call.beta_init.to_vec(),
                    efficiencies: None,
                    condition_se: Vec::new(),
                    means: Vec::new(),
                    residuals: Vec::new(),
                })
            })
            .unwrap();
            assert_eq!(calls, expected_calls);
            assert!(!result.status.converged());
            assert_eq!(result.history.len(), expected_calls);
        }
        let mut calls = 0;
        let result = setup_and_fit_with_runner(&input(), |_| {
            calls += 1;
            Err(numerical("initial normal solve failure"))
        });
        assert!(result.is_err());
        assert_eq!(calls, 1);
    }
}
