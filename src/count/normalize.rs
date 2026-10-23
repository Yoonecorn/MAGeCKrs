use crate::cli::NormMethod;
use crate::data::count_matrix::CountMatrix;
use crate::error::MageckError;

/// Normalizes the count matrix according to the specified method.
/// Returns a new `CountMatrix` with the normalized counts.
pub fn normalize(
    matrix: &CountMatrix,
    method: NormMethod,
    control_guides: Option<&[usize]>,
) -> Result<CountMatrix, MageckError> {
    if matrix.is_empty() || matrix.n_samples() == 0 || matrix.n_guides() == 0 {
        return Ok(CountMatrix::from_vec(
            matrix.as_slice().to_vec(),
            matrix.n_samples(),
            matrix.n_guides(),
        )
        .unwrap());
    }

    // Validate inputs
    for &val in matrix.as_slice() {
        if !val.is_finite() || val < 0.0 {
            return Err(MageckError::Config(
                "Matrix contains negative, NaN, or infinite counts".to_string(),
            ));
        }
    }

    match method {
        NormMethod::None => Ok(CountMatrix::from_vec(
            matrix.as_slice().to_vec(),
            matrix.n_samples(),
            matrix.n_guides(),
        )
        .unwrap()),
        NormMethod::Total => normalize_total(matrix),
        NormMethod::Median => normalize_median_ratio(matrix, None),
        NormMethod::Control => {
            if let Some(ctrls) = control_guides {
                // Validate control guides
                for &g in ctrls {
                    if g >= matrix.n_guides() {
                        return Err(MageckError::Config(format!(
                            "Control guide index {} out of bounds",
                            g
                        )));
                    }
                }
                normalize_median_ratio(matrix, Some(ctrls))
            } else {
                Err(MageckError::Config(
                    "Control guides must be provided for control normalization".to_string(),
                ))
            }
        }
    }
}

fn normalize_total(matrix: &CountMatrix) -> Result<CountMatrix, MageckError> {
    let n_samples = matrix.n_samples();
    let mut sample_totals = Vec::with_capacity(n_samples);
    let mut total_sum = 0.0;

    for s in 0..n_samples {
        let sum: f64 = matrix.sample_slice(s).iter().sum();
        if !sum.is_finite() {
            return Err(MageckError::Config(format!(
                "Sample {} total count is not finite",
                s
            )));
        }
        if sum == 0.0 {
            return Err(MageckError::Config(format!(
                "Sample {} has zero total count; cannot compute total-count size factor",
                s
            )));
        }
        sample_totals.push(sum);
        total_sum += sum;
    }

    if !total_sum.is_finite() || total_sum == 0.0 {
        return Err(MageckError::Config(
            "Total count across all samples is not positive finite".to_string(),
        ));
    }

    let mean_total = total_sum / n_samples as f64;
    let mut sf = Vec::with_capacity(n_samples);
    for sum in sample_totals {
        sf.push(sum / mean_total);
    }

    let mut out = CountMatrix::new(n_samples, matrix.n_guides());
    for (s, &size_factor) in sf.iter().enumerate() {
        let in_slice = matrix.sample_slice(s);
        let out_slice = out.sample_slice_mut(s);
        for (in_val, out_val) in in_slice.iter().zip(out_slice.iter_mut()) {
            *out_val = in_val / size_factor;
        }
    }

    Ok(out)
}

fn normalize_median_ratio(
    matrix: &CountMatrix,
    subset_guides: Option<&[usize]>,
) -> Result<CountMatrix, MageckError> {
    let n_samples = matrix.n_samples();
    let n_guides = matrix.n_guides();

    let guide_iter: Vec<usize> = match subset_guides {
        Some(guides) => guides.to_vec(),
        None => (0..n_guides).collect(),
    };

    let mut valid_guides = Vec::new();
    let mut geometric_means = Vec::new();

    for &g in &guide_iter {
        let mut has_zero = false;
        let mut sum_log = 0.0;
        for s in 0..n_samples {
            let val = matrix.entry(s, g);
            if val <= 0.0 {
                has_zero = true;
                break;
            }
            sum_log += val.ln();
        }

        if !has_zero {
            let gm = (sum_log / n_samples as f64).exp();
            if gm > 0.0 && gm.is_finite() {
                valid_guides.push(g);
                geometric_means.push(gm);
            }
        }
    }

    if valid_guides.is_empty() {
        return Err(MageckError::Config(
            "Every guide has at least one zero (or subnormal) count; cannot compute median-ratio size factors"
                .to_string(),
        ));
    }

    let mut size_factors = Vec::with_capacity(n_samples);
    for s in 0..n_samples {
        let mut ratios = Vec::with_capacity(valid_guides.len());
        for (i, &g) in valid_guides.iter().enumerate() {
            ratios.push(matrix.entry(s, g) / geometric_means[i]);
        }
        let sf = median(&mut ratios);
        if sf <= 0.0 || !sf.is_finite() {
            return Err(MageckError::Config(format!(
                "Size factor for sample {} is not a positive finite number (got {})",
                s, sf
            )));
        }
        size_factors.push(sf);
    }

    let mut out = CountMatrix::new(n_samples, n_guides);
    for (s, &sf) in size_factors.iter().enumerate() {
        let in_slice = matrix.sample_slice(s);
        let out_slice = out.sample_slice_mut(s);
        for (in_val, out_val) in in_slice.iter().zip(out_slice.iter_mut()) {
            *out_val = in_val / sf;
        }
    }

    Ok(out)
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_unstable_by(|a, b| a.total_cmp(b));
    let len = values.len();
    if len == 0 {
        return f64::NAN;
    }
    if len % 2 == 1 {
        values[len / 2]
    } else {
        let a = values[len / 2 - 1];
        let b = values[len / 2];
        a + (b - a) / 2.0
    }
}
