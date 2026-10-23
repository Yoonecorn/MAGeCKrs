//! Private fitted-result projection and two-sided Wald statistics.

use super::irls::IrlsStatus;
use crate::error::MageckError;

// Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier.
//
// Permission to use, copy, modify, and distribute this software is freely
// granted, provided that this notice is preserved.
//
// Deterministic Rust translation of the Cephes Math Library Release 2.0
// `ndtr.c` rational-polynomial normal-tail implementation.
const CEPHES_P: [f64; 9] = [
    2.46196981473530512524e-10,
    5.64189564831068821977e-1,
    7.46321056442269912687e0,
    4.86371970985681366614e1,
    1.96520832956077098242e2,
    5.26445194995477358631e2,
    9.34528527171957607540e2,
    1.02755188689515710272e3,
    5.57535335369399327526e2,
];
const CEPHES_Q: [f64; 8] = [
    1.32281951154744992508e1,
    8.67072140885989742329e1,
    3.54937778887819891062e2,
    9.75708501743205489753e2,
    1.82390916687909736289e3,
    2.24633760818710981792e3,
    1.65666309194161350182e3,
    5.57535340817727675546e2,
];
const CEPHES_R: [f64; 6] = [
    5.64189583547755073984e-1,
    1.27536670759978104416e0,
    5.01905042251180477414e0,
    6.16021097993053585195e0,
    7.40974269950448939160e0,
    2.97886665372100240670e0,
];
const CEPHES_S: [f64; 6] = [
    2.26052863220117276590e0,
    9.39603524938001434673e0,
    1.20489539808096656605e1,
    1.70814450747565897222e1,
    9.60896809063285878198e0,
    3.36907645100081516050e0,
];
const CEPHES_T: [f64; 5] = [
    9.60497373987051638749e0,
    9.00260197203842689217e1,
    2.23200534594684319226e3,
    7.00332514112805075473e3,
    5.55923013010394962768e4,
];
const CEPHES_U: [f64; 5] = [
    3.35617141647503099647e1,
    5.21357949780152679795e2,
    4.59432382970980127987e3,
    2.26290000613890934246e4,
    4.92673942608635921086e4,
];
const CEPHES_MAXLOG: f64 = -709.782712893383996732;

fn cephes_polevl(x: f64, coefficients: &[f64]) -> f64 {
    coefficients[1..]
        .iter()
        .fold(coefficients[0], |value, coefficient| {
            value * x + coefficient
        })
}

fn cephes_p1evl(x: f64, coefficients: &[f64]) -> f64 {
    coefficients[1..]
        .iter()
        .fold(x + coefficients[0], |value, coefficient| {
            value * x + coefficient
        })
}

fn cephes_erf(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        return -cephes_erf(-x);
    }
    if x > 1.0 {
        return 1.0 - cephes_erfc(x);
    }
    let squared = x * x;
    x * cephes_polevl(squared, &CEPHES_T) / cephes_p1evl(squared, &CEPHES_U)
}

fn cephes_erfc(a: f64) -> f64 {
    if a.is_nan() {
        return f64::NAN;
    }
    let absolute = a.abs();
    if absolute < 1.0 {
        return 1.0 - cephes_erf(a);
    }
    let exponent = -a * a;
    if exponent < CEPHES_MAXLOG {
        return if a < 0.0 { 2.0 } else { 0.0 };
    }
    let (numerator, denominator) = if absolute < 8.0 {
        (
            cephes_polevl(absolute, &CEPHES_P),
            cephes_p1evl(absolute, &CEPHES_Q),
        )
    } else {
        (
            cephes_polevl(absolute, &CEPHES_R),
            cephes_p1evl(absolute, &CEPHES_S),
        )
    };
    let mut value = exponent.exp() * numerator / denominator;
    if a < 0.0 {
        value = 2.0 - value;
    }
    if value != 0.0 {
        value
    } else if a < 0.0 {
        2.0
    } else {
        0.0
    }
}

fn cephes_ndtr(a: f64) -> f64 {
    if a.is_nan() {
        return f64::NAN;
    }
    let x = a * std::f64::consts::FRAC_1_SQRT_2;
    let absolute = x.abs();
    if absolute < std::f64::consts::FRAC_1_SQRT_2 {
        0.5 + 0.5 * cephes_erf(x)
    } else {
        let mut value = 0.5 * cephes_erfc(absolute);
        if x > 0.0 {
            value = 1.0 - value;
        }
        value
    }
}

fn wald_p_value(z_score: f64) -> f64 {
    2.0 * cephes_ndtr(-z_score.abs())
}

#[derive(Clone, Debug)]
pub(crate) struct WaldInput<'a> {
    pub(crate) full_beta: &'a [f64],
    pub(crate) n_guides: usize,
    pub(crate) n_conditions: usize,
    pub(crate) condition_covariance: &'a [f64],
    pub(crate) status: Option<IrlsStatus>,
    pub(crate) attempts: usize,
    pub(crate) means: &'a [f64],
    pub(crate) residuals: &'a [f64],
    pub(crate) efficiencies: Option<&'a [f64]>,
    pub(crate) condition_order: &'a [String],
}

#[derive(Clone, Debug)]
pub(crate) struct WaldResult {
    pub(crate) condition_beta: Vec<f64>,
    pub(crate) condition_covariance: Vec<f64>,
    pub(crate) standard_error: Vec<f64>,
    pub(crate) z_score: Vec<f64>,
    pub(crate) two_sided_p_value: Vec<f64>,
    pub(crate) status: IrlsStatus,
    pub(crate) attempts: usize,
    pub(crate) means: Vec<f64>,
    pub(crate) residuals: Vec<f64>,
    pub(crate) efficiencies: Option<Vec<f64>>,
    pub(crate) condition_order: Vec<String>,
}

fn numerical(message: impl Into<String>) -> MageckError {
    MageckError::Math(message.into())
}

fn product(a: usize, b: usize, label: &str) -> Result<usize, MageckError> {
    a.checked_mul(b)
        .ok_or_else(|| numerical(format!("wald {label} dimension overflow")))
}

fn finite_all(values: &[f64], label: &str) -> Result<(), MageckError> {
    if values.iter().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err(numerical(format!("wald non-finite {label}")))
    }
}

pub(crate) fn project_wald(input: WaldInput<'_>) -> Result<WaldResult, MageckError> {
    if input.status.is_some_and(|s| !s.converged()) {
        return Err(MageckError::Math(
            "Wald requires converged eligible fit".into(),
        ));
    }

    let g = input.n_guides;
    let b = input.n_conditions;
    if g == 0 || b == 0 {
        return Err(numerical(
            "wald guide and condition counts must be non-zero",
        ));
    }
    let beta_len = g
        .checked_add(b)
        .ok_or_else(|| numerical("wald coefficient dimension overflow"))?;
    if input.full_beta.len() != beta_len {
        return Err(numerical("wald beta shape mismatch"));
    }
    if input.condition_covariance.len() != product(b, b, "covariance")? {
        return Err(numerical("wald covariance shape mismatch"));
    }
    if input.condition_order.len() != b {
        return Err(numerical("wald condition-order shape mismatch"));
    }
    let status = input
        .status
        .ok_or_else(|| numerical("wald missing or invalid fitted status"))?;
    if input.attempts == 0 {
        return Err(numerical("wald attempts must be non-zero"));
    }
    finite_all(input.full_beta, "beta")?;
    finite_all(input.condition_covariance, "covariance")?;
    finite_all(input.means, "means")?;
    if input.means.iter().any(|x| *x <= 0.0) {
        return Err(numerical("wald means must be positive"));
    }
    finite_all(input.residuals, "residuals")?;
    if let Some(efficiencies) = input.efficiencies {
        if efficiencies.len() != g {
            return Err(numerical("wald efficiency shape mismatch"));
        }
        finite_all(efficiencies, "efficiencies")?;
    }

    let condition_beta = input.full_beta[g..].to_vec();
    let mut standard_error = Vec::with_capacity(b);
    let mut z_score = Vec::with_capacity(b);
    let mut two_sided_p_value = Vec::with_capacity(b);
    for (index, &beta) in condition_beta.iter().enumerate() {
        let variance = input.condition_covariance[index * b + index];
        if variance <= 0.0 {
            return Err(numerical("wald covariance diagonal must be positive"));
        }
        let se = variance.sqrt();
        if !se.is_finite() || se <= 0.0 {
            return Err(numerical("wald non-finite standard error"));
        }
        let z = beta / se;
        if !z.is_finite() {
            return Err(numerical("wald non-finite z score"));
        }
        let p_value = wald_p_value(z);
        if !p_value.is_finite() || !(0.0..=1.0).contains(&p_value) {
            return Err(numerical("wald p-value outside unit interval"));
        }
        standard_error.push(se);
        z_score.push(z);
        two_sided_p_value.push(p_value);
    }

    Ok(WaldResult {
        condition_beta,
        condition_covariance: input.condition_covariance.to_vec(),
        standard_error,
        z_score,
        two_sided_p_value,
        status,
        attempts: input.attempts,
        means: input.means.to_vec(),
        residuals: input.residuals.to_vec(),
        efficiencies: input.efficiencies.map(ToOwned::to_owned),
        condition_order: input.condition_order.to_vec(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(names: &[&str]) -> Vec<String> {
        names.iter().map(|name| (*name).to_owned()).collect()
    }

    fn assert_close(actual: &[f64], expected: &[f64], tolerance: f64) {
        assert_eq!(actual.len(), expected.len());
        for (&actual, &expected) in actual.iter().zip(expected) {
            assert!(
                (actual - expected).abs() <= tolerance,
                "{actual:.17e} != {expected:.17e}"
            );
        }
    }

    #[test]
    fn positive_one_condition() {
        let order = names(&["treatment_A"]);
        let result = project_wald(WaldInput {
            full_beta: &[3.0, 3.1, 0.75],
            n_guides: 2,
            n_conditions: 1,
            condition_covariance: &[0.0625],
            status: Some(IrlsStatus::Converged),
            attempts: 1,
            means: &[20.0, 22.0, 18.0, 21.0],
            residuals: &[0.1, -0.2, 0.3, -0.1],
            efficiencies: Some(&[1.0, 0.8]),
            condition_order: &order,
        })
        .unwrap();
        assert_close(&result.condition_beta, &[0.75], 0.0);
        assert_close(&result.condition_covariance, &[0.0625], 0.0);
        assert_close(&result.standard_error, &[0.25], 1e-15);
        assert_close(&result.z_score, &[3.0], 1e-15);
        assert_close(&result.two_sided_p_value, &[0.002699796063260186], 1e-14);
        assert_eq!(result.status, IrlsStatus::Converged);
        assert_eq!(result.attempts, 1);
        assert_eq!(result.condition_order, order);
        assert_eq!(result.means, vec![20.0, 22.0, 18.0, 21.0]);
        assert_eq!(result.residuals, vec![0.1, -0.2, 0.3, -0.1]);
        assert_eq!(result.efficiencies, Some(vec![1.0, 0.8]));
    }

    #[test]
    fn ordered_three_conditions() {
        let order = names(&["negative_first", "small_middle", "positive_last"]);
        let result = project_wald(WaldInput {
            full_beta: &[1.5, 1.6, -1.2, 0.05, 2.0],
            n_guides: 2,
            n_conditions: 3,
            condition_covariance: &[0.04, 0.01, -0.005, 0.01, 0.25, 0.03, -0.005, 0.03, 0.01],
            status: Some(IrlsStatus::Converged),
            attempts: 1,
            means: &[5.0, 5.4, 4.8, 5.3, 6.1, 6.4, 3.2, 3.5],
            residuals: &[0.2, -0.1, 0.3, -0.2, 0.1, -0.1, 0.4, -0.3],
            efficiencies: Some(&[0.9, 0.6]),
            condition_order: &order,
        })
        .unwrap();
        assert_close(&result.condition_beta, &[-1.2, 0.05, 2.0], 0.0);
        assert_close(
            &result.condition_covariance,
            &[0.04, 0.01, -0.005, 0.01, 0.25, 0.03, -0.005, 0.03, 0.01],
            0.0,
        );
        assert_close(&result.standard_error, &[0.2, 0.5, 0.1], 1e-15);
        assert_close(&result.z_score, &[-6.0, 0.1, 20.0], 1e-14);
        assert_close(
            &result.two_sided_p_value,
            &[
                1.973175290075403e-09,
                0.920344325445942,
                5.50724823721231e-89,
            ],
            1e-13,
        );
        assert_eq!(result.condition_order, order);
    }

    #[test]
    fn retry_exhaustion_retained() {
        let order = names(&["retry_negative", "retry_positive"]);
        let result = project_wald(WaldInput {
            full_beta: &[2.1, 1.9, 2.2, -0.4, 0.9],
            n_guides: 3,
            n_conditions: 2,
            condition_covariance: &[0.16, -0.02, -0.02, 0.09],
            status: Some(IrlsStatus::IterationLimitReached),
            attempts: 11,
            means: &[8.0, 7.5, 8.2, 6.0, 5.8, 6.1],
            residuals: &[0.3, -0.3, 0.2, 0.4, -0.2, -0.1],
            efficiencies: Some(&[0.7, 0.8, 1.0]),
            condition_order: &order,
        });
        assert!(result.is_err());
    }

    #[test]
    fn strict_validation_rejects_all_documented_boundaries() {
        let order = names(&["condition"]);
        let valid = || WaldInput {
            full_beta: &[1.0, 0.5],
            n_guides: 1,
            n_conditions: 1,
            condition_covariance: &[0.25],
            status: Some(IrlsStatus::Converged),
            attempts: 1,
            means: &[1.0],
            residuals: &[0.0],
            efficiencies: Some(&[1.0]),
            condition_order: &order,
        };
        assert!(project_wald(WaldInput {
            n_conditions: 0,
            condition_order: &[],
            ..valid()
        })
        .is_err());
        assert!(project_wald(WaldInput {
            full_beta: &[1.0],
            ..valid()
        })
        .is_err());
        assert!(project_wald(WaldInput {
            condition_covariance: &[0.25, 0.0],
            ..valid()
        })
        .is_err());
        assert!(project_wald(WaldInput {
            full_beta: &[f64::NAN, 0.5],
            ..valid()
        })
        .is_err());
        assert!(project_wald(WaldInput {
            condition_covariance: &[f64::INFINITY],
            ..valid()
        })
        .is_err());
        assert!(project_wald(WaldInput {
            condition_covariance: &[-0.01],
            ..valid()
        })
        .is_err());
        assert!(project_wald(WaldInput {
            condition_covariance: &[0.0],
            ..valid()
        })
        .is_err());
        assert!(project_wald(WaldInput {
            full_beta: &[1.0e308, 1.0e308],
            condition_covariance: &[1.0e-308],
            ..valid()
        })
        .is_err());
        assert!(project_wald(WaldInput {
            status: None,
            ..valid()
        })
        .is_err());
        assert!(project_wald(WaldInput {
            condition_order: &[],
            ..valid()
        })
        .is_err());
        assert!(project_wald(WaldInput {
            efficiencies: Some(&[]),
            ..valid()
        })
        .is_err());
    }
}
