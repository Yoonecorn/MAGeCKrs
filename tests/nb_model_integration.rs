use mageckrs::data::count_matrix::CountMatrix;
use mageckrs::mageck_test::{fit_meanvar, POISSON_DISPERSION, POISSON_SIZE};

fn make_matrix(counts: &[&[f64]]) -> CountMatrix {
    let n_samples = counts.len();
    let n_guides = counts[0].len();
    let mut flat = Vec::with_capacity(n_samples * n_guides);
    for row in counts {
        flat.extend_from_slice(row);
    }
    CountMatrix::from_vec(flat, n_samples, n_guides).unwrap()
}

#[test]
fn test_empty_control_list_rejected() {
    let matrix = make_matrix(&[&[1.0, 2.0], &[3.0, 4.0]]);
    let res = fit_meanvar(&matrix, &[]);
    assert!(res.is_err(), "Empty control list should be rejected");
}

#[test]
fn test_one_control_rejected() {
    let matrix = make_matrix(&[&[1.0, 2.0], &[3.0, 4.0]]);
    let res = fit_meanvar(&matrix, &[0]);
    assert!(res.is_err(), "Single control should be rejected");
}

#[test]
fn test_duplicate_control_rejected() {
    let matrix = make_matrix(&[&[1.0, 2.0], &[3.0, 4.0], &[5.0, 6.0]]);
    let res = fit_meanvar(&matrix, &[0, 0]);
    assert!(res.is_err(), "Duplicate controls should be rejected");
}

#[test]
fn test_out_of_range_control_rejected() {
    let matrix = make_matrix(&[&[1.0, 2.0], &[3.0, 4.0]]);
    let res = fit_meanvar(&matrix, &[0, 2]);
    assert!(res.is_err(), "Out of range control should be rejected");
}

#[test]
fn test_two_controls_mean_and_variance() {
    let matrix = make_matrix(&[&[2.0], &[4.0]]);
    let model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    assert_eq!(model.mean()[0], 3.0);
    assert_eq!(model.variance()[0], 2.0);
}

#[test]
fn test_three_controls_mean_and_variance() {
    let matrix = make_matrix(&[&[2.0], &[4.0], &[6.0]]);
    let model = fit_meanvar(&matrix, &[0, 1, 2]).unwrap();
    assert_eq!(model.mean()[0], 4.0);
    assert_eq!(model.variance()[0], 4.0);
}

#[test]
fn test_population_variance_not_used() {
    let matrix = make_matrix(&[&[2.0], &[4.0], &[6.0]]);
    let model = fit_meanvar(&matrix, &[0, 1, 2]).unwrap();
    assert_ne!(model.variance()[0], 8.0 / 3.0);
    assert_eq!(model.variance()[0], 4.0);
}

#[test]
fn test_raw_dispersion_eligibility() {
    let matrix = make_matrix(&[
        &[0.0, 5.0, 2.0, 2.0],
        &[0.0, 5.0, 4.0, 10.0],
        &[0.0, 5.0, 6.0, 18.0],
    ]);
    let model = fit_meanvar(&matrix, &[0, 1, 2]).unwrap();
    assert_eq!(model.raw_dispersion()[0], None);
    assert_eq!(model.raw_dispersion()[1], None);
    assert_eq!(model.raw_dispersion()[2], None);
    let raw3 = model.raw_dispersion()[3].unwrap();
    assert!((raw3 - 0.54).abs() < 1e-9);
}

#[test]
fn test_zero_mean_fallback() {
    let matrix = make_matrix(&[&[0.0], &[0.0]]);
    let model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    assert_eq!(model.mean()[0], 0.0);
    assert_eq!(model.raw_dispersion()[0], None);
    assert_eq!(model.final_dispersion()[0], POISSON_DISPERSION);
    assert_eq!(model.size()[0], POISSON_SIZE);
}

#[test]
fn test_guide_order_preservation() {
    let matrix = make_matrix(&[&[10.0, 0.0, 5.0], &[20.0, 0.0, 5.0]]);
    let model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    assert_eq!(model.mean(), &[15.0, 0.0, 5.0]);
}

#[test]
fn test_only_control_columns_influence_fit() {
    let matrix = make_matrix(&[&[10.0], &[1000.0], &[20.0]]);
    let model = fit_meanvar(&matrix, &[0, 2]).unwrap();
    assert_eq!(model.mean()[0], 15.0);
}

#[test]
fn test_changing_non_control_values_does_not_affect_model() {
    let m1 = make_matrix(&[&[10.0], &[1000.0], &[20.0]]);
    let model1 = fit_meanvar(&m1, &[0, 2]).unwrap();
    let m2 = make_matrix(&[&[10.0], &[9999.0], &[20.0]]);
    let model2 = fit_meanvar(&m2, &[0, 2]).unwrap();
    assert_eq!(model1.mean(), model2.mean());
    assert_eq!(model1.variance(), model2.variance());
}

#[test]
fn test_zero_centered_x_variance_uses_median() {
    let matrix = make_matrix(&[&[2.0, 2.0, 2.0], &[10.0, 10.0, 10.0], &[18.0, 18.0, 18.0]]);
    let model = fit_meanvar(&matrix, &[0, 1, 2]).unwrap();
    assert_eq!(model.eligible_guides_count(), 3);
    assert!(
        model.regression_coefficients().is_none(),
        "Regression should fail due to zero variance in x"
    );
    assert!(model.median_dispersion() > 0.0);
    assert_eq!(model.fitted_dispersion()[0], model.median_dispersion());
}

#[test]
fn test_ols_slope_intercept() {
    // G0: mean 10, v 20 (d = sqrt(20) = 4.47213595499958) => raw 0.1
    // G1: mean 20, v 100 (d = 10) => raw 0.2
    // G2: mean 30, v 210 (d = sqrt(210) = 14.491376746189438) => raw 0.2
    // G3: mean 40, v 520 (d = sqrt(520) = 22.80350850198276) => raw 0.3
    let matrix = make_matrix(&[
        &[10.0, 20.0, 30.0, 40.0],
        &[
            10.0 + 4.47213595499958,
            20.0 + 10.0,
            30.0 + 14.491376746189438,
            40.0 + 22.80350850198276,
        ],
        &[
            10.0 - 4.47213595499958,
            20.0 - 10.0,
            30.0 - 14.491376746189438,
            40.0 - 22.80350850198276,
        ],
    ]);
    let model = fit_meanvar(&matrix, &[0, 1, 2]).unwrap();
    assert_eq!(model.eligible_guides_count(), 4);
    let coeffs = model.regression_coefficients().expect("OLS should fit");

    let xs = [10.0_f64.ln(), 20.0_f64.ln(), 30.0_f64.ln(), 40.0_f64.ln()];
    let ys = [0.1_f64.ln(), 0.2_f64.ln(), 0.2_f64.ln(), 0.3_f64.ln()];
    let mean_x = xs.iter().sum::<f64>() / 4.0;
    let mean_y = ys.iter().sum::<f64>() / 4.0;

    let mut ss_x = 0.0;
    let mut cov_xy = 0.0;
    for i in 0..4 {
        ss_x += (xs[i] - mean_x).powi(2);
        cov_xy += (xs[i] - mean_x) * (ys[i] - mean_y);
    }
    let expected_slope = cov_xy / ss_x;
    let expected_intercept = mean_y - expected_slope * mean_x;

    assert!((coeffs.slope - expected_slope).abs() < 1e-5);
    assert!((coeffs.intercept - expected_intercept).abs() < 1e-5);
}

#[test]
fn test_median_fallback_0_guides() {
    let matrix = make_matrix(&[&[1.0, 2.0], &[1.0, 2.0]]);
    let model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    assert_eq!(model.eligible_guides_count(), 0);
    assert_eq!(model.median_dispersion(), POISSON_DISPERSION);
}

#[test]
fn test_median_fallback_1_guide() {
    // 1 eligible guide: mean=10, var=60 (d=sqrt(60)=7.745966692414834) => raw=0.5
    let matrix = make_matrix(&[
        &[10.0, 1.0],
        &[10.0 + 7.745966692414834, 1.0],
        &[10.0 - 7.745966692414834, 1.0],
    ]);
    let model = fit_meanvar(&matrix, &[0, 1, 2]).unwrap();
    assert_eq!(model.eligible_guides_count(), 1);
    assert!((model.median_dispersion() - 0.5).abs() < 1e-5);
}

#[test]
fn test_median_fallback_2_guides() {
    // G0: mean 10, var 60 (d=7.745966692414834) => raw 0.5
    // G1: mean 20, var 140 (d=11.832159566199232) => raw 0.3
    let matrix = make_matrix(&[
        &[10.0, 20.0],
        &[10.0 + 7.745966692414834, 20.0 + 11.832159566199232],
        &[10.0 - 7.745966692414834, 20.0 - 11.832159566199232],
    ]);
    let model = fit_meanvar(&matrix, &[0, 1, 2]).unwrap();
    assert_eq!(model.eligible_guides_count(), 2);
    // Median of [0.3, 0.5] is 0.4
    assert!((model.median_dispersion() - 0.4).abs() < 1e-5);
}

#[test]
fn test_conservative_maximum() {
    let matrix = make_matrix(&[
        &[10.0, 20.0, 30.0, 40.0],
        &[
            10.0 + 4.47213595499958,
            20.0 + 10.0,
            30.0 + 20.0,
            40.0 + 22.80350850198276,
        ],
        &[
            10.0 - 4.47213595499958,
            20.0 - 10.0,
            30.0 - 20.0,
            40.0 - 22.80350850198276,
        ],
    ]);
    let model = fit_meanvar(&matrix, &[0, 1, 2]).unwrap();
    assert!(model.regression_coefficients().is_some());
    for g in 0..4 {
        let fitted = model.fitted_dispersion()[g];
        let final_d = model.final_dispersion()[g];
        if let Some(raw) = model.raw_dispersion()[g] {
            assert_eq!(final_d, raw.max(fitted));
        } else {
            assert_eq!(final_d, fitted);
        }
    }
}

#[test]
fn test_positive_mean_without_raw_dispersion_uses_fitted() {
    let matrix = make_matrix(&[
        &[10.0, 20.0, 30.0, 5.0],
        &[
            10.0 + 4.47213595499958,
            20.0 + 10.0,
            30.0 + 14.491376746189438,
            5.0,
        ],
        &[
            10.0 - 4.47213595499958,
            20.0 - 10.0,
            30.0 - 14.491376746189438,
            5.0,
        ],
    ]);
    let model = fit_meanvar(&matrix, &[0, 1, 2]).unwrap();
    assert!(model.regression_coefficients().is_some());
    assert_eq!(model.raw_dispersion()[3], None);
    assert_eq!(model.final_dispersion()[3], model.fitted_dispersion()[3]);
}

#[test]
fn test_all_final_dispersions_and_sizes_positive_and_finite() {
    let matrix = make_matrix(&[&[10.0, 0.0, 1.0], &[12.0, 0.0, 1.0], &[8.0, 0.0, 1.0]]);
    let model = fit_meanvar(&matrix, &[0, 1, 2]).unwrap();
    for g in 0..3 {
        assert!(model.final_dispersion()[g].is_finite());
        assert!(model.final_dispersion()[g] > 0.0);
        assert!(model.size()[g].is_finite());
        assert!(model.size()[g] > 0.0);
        assert!((model.size()[g] * model.final_dispersion()[g] - 1.0).abs() < 1e-9);
    }
}

#[test]
fn test_result_vector_lengths() {
    let matrix = make_matrix(&[&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]]);
    let model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    assert_eq!(model.mean().len(), 3);
    assert_eq!(model.variance().len(), 3);
    assert_eq!(model.raw_dispersion().len(), 3);
    assert_eq!(model.fitted_dispersion().len(), 3);
    assert_eq!(model.final_dispersion().len(), 3);
    assert_eq!(model.size().len(), 3);
}

#[test]
fn test_invalid_counts() {
    let matrix_neg = make_matrix(&[&[-1.0], &[-2.0]]);
    assert!(fit_meanvar(&matrix_neg, &[0, 1]).is_err());
    let matrix_nan = make_matrix(&[&[f64::NAN], &[2.0]]);
    assert!(fit_meanvar(&matrix_nan, &[0, 1]).is_err());
    let matrix_inf = make_matrix(&[&[f64::INFINITY], &[2.0]]);
    assert!(fit_meanvar(&matrix_inf, &[0, 1]).is_err());
    let matrix_neg_inf = make_matrix(&[&[f64::NEG_INFINITY], &[2.0]]);
    assert!(fit_meanvar(&matrix_neg_inf, &[0, 1]).is_err());
}

#[test]
fn test_zero_guide_matrix() {
    let matrix = CountMatrix::from_vec(vec![], 2, 0).unwrap();
    let res = fit_meanvar(&matrix, &[0, 1]);
    assert!(res.is_err(), "Zero-guide matrix should be rejected");
}

#[test]
fn test_source_matrix_unchanged() {
    let original = make_matrix(&[&[1.0, 2.0], &[3.0, 4.0]]);
    let to_use = make_matrix(&[&[1.0, 2.0], &[3.0, 4.0]]);
    let _ = fit_meanvar(&to_use, &[0, 1]).unwrap();
    assert_eq!(original.as_slice(), to_use.as_slice());
}

#[test]
fn test_invalid_predicted_dispersion_fallback() {
    let matrix = make_matrix(&[
        &[10.0, 20.0, 30.0, 1e-300],
        &[
            10.0 + 4.47213595499958,
            20.0 + 10.0,
            30.0 + 14.491376746189438,
            1e-300,
        ],
        &[
            10.0 - 4.47213595499958,
            20.0 - 10.0,
            30.0 - 14.491376746189438,
            1e-300,
        ],
    ]);
    let model = fit_meanvar(&matrix, &[0, 1, 2]).unwrap();
    assert!(model.regression_coefficients().is_some());
    let fitted3 = model.fitted_dispersion()[3];
    assert!(fitted3 > 0.0, "Fitted dispersion must be positive");
    assert!(fitted3.is_finite(), "Fitted dispersion must be finite");
}
