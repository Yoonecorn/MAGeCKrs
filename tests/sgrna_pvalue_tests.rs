use mageckrs::data::CountMatrix;
use mageckrs::error::MageckError;
use mageckrs::mageck_test::nb_model::fit_meanvar;
use mageckrs::mageck_test::sgrna_pvalue::calculate_sgrna_pvalues;
use mageckrs::mageck_test::sgrna_pvalue::{get_qnorm, get_qnorm_table, parse_qnorm_table};
use statrs::distribution::{ContinuousCDF, Normal};

fn create_matrix(samples: usize, guides: usize, data: &[f64]) -> CountMatrix {
    assert_eq!(data.len(), guides * samples);
    let mut matrix = CountMatrix::new(samples, guides);
    for g in 0..guides {
        for s in 0..samples {
            matrix.set(s, g, data[g * samples + s]);
        }
    }
    matrix
}

#[test]
fn test_source_faithful_differential() {
    // 5 guides, 4 samples (0, 1: control, 2, 3: treatment)
    let data = vec![
        100.0, 100.0, 50.0, 50.0, // Guide 0: control=100, trt=50 (Depleted)
        100.0, 100.0, 200.0, 200.0, // Guide 1: control=100, trt=200 (Enriched)
        50.0, 50.0, 50.0, 50.0, // Guide 2: control=50, trt=50 (Equal)
        10.0, 10.0, 0.0, 0.0, // Guide 3: control=10, trt=0 (Zero trt)
        0.0, 0.0, 10.0, 10.0, // Guide 4: control=0, trt=10 (Zero control)
    ];
    let matrix = create_matrix(4, 5, &data);
    let nb_model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    let res = calculate_sgrna_pvalues(&matrix, &[2, 3], &nb_model).unwrap();

    assert_eq!(res.len(), 5);

    // Guide Order Test (Repair 13)
    let guide_indices: Vec<usize> = res.iter().map(|r| r.guide_index).collect();
    assert_eq!(guide_indices, vec![0, 1, 2, 3, 4]);

    // Guide 0: Depleted
    // Repair 12: Strict p-value fixture based on crisprFunction.py
    // mean=100, disp=0.0001 -> var=101.0 -> t_theta = -50 / sqrt(101) = -4.9751859
    assert!((res[0].control_mean - 100.0).abs() < 1e-9);
    assert!((res[0].treatment_mean - 50.0).abs() < 1e-9);
    let expected_p_lower = 3.259_251_753_695_64e-7;
    let p_lower_abs_err = (res[0].p_lower - expected_p_lower).abs();
    let p_lower_rel_err = p_lower_abs_err / expected_p_lower;
    assert!(p_lower_abs_err < 1e-12 || p_lower_rel_err < 1e-5);

    // Check score negative
    assert!((res[0].score_negative - (-4.9751859)).abs() < 1e-4);

    // Guide 1: Enriched
    assert!((res[1].p_higher).abs() < 1e-9);
    assert!((res[1].score_positive - (-9.9503719)).abs() < 1e-5);

    // Guide 4: Zero control is clamped to min positive (10.0)
    assert!((res[4].control_mean - 10.0).abs() < 1e-9); // Clamped from 0 to 10
    assert!((res[4].treatment_mean - 10.0).abs() < 1e-9);
    assert!((res[4].score_negative - 0.0).abs() < 1e-9);

    // Guide 3: Zero treatment is clamped to min positive inside getNormalPValue, but ttmean remains 0.0
    assert!((res[3].treatment_mean - 0.0).abs() < 1e-9);
    assert!(res[3].score_negative < 0.0);

    // Check rank directions
    assert!((res[0].rank_negative - 0.2).abs() < 1e-9);
    assert!((res[3].rank_negative - 0.4).abs() < 1e-9);
    assert!((res[2].rank_negative - 0.7).abs() < 1e-9);
    assert!((res[4].rank_negative - 0.7).abs() < 1e-9);
    assert!((res[1].rank_negative - 1.0).abs() < 1e-9);

    assert!((res[1].rank_positive - 0.2).abs() < 1e-9);
    assert!((res[2].rank_positive - 0.5).abs() < 1e-9);
    assert!((res[4].rank_positive - 0.5).abs() < 1e-9);
    assert!((res[3].rank_positive - 0.8).abs() < 1e-9);
    assert!((res[0].rank_positive - 1.0).abs() < 1e-9);
}

#[test]
fn test_invalid_variance() {
    let mut matrix = create_matrix(4, 1, &[100.0, 100.0, 100.0, 100.0]);
    matrix.set(0, 0, f64::NAN); // Induce an error
    let err = fit_meanvar(&matrix, &[0, 1]).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));
}

#[test]
fn test_no_positive_treatment_mean() {
    let data = vec![100.0, 100.0, 0.0, 0.0];
    let matrix = create_matrix(4, 1, &data);
    let nb_model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    let err = calculate_sgrna_pvalues(&matrix, &[2, 3], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));
}

#[test]
fn test_no_positive_control_mean() {
    let data = vec![0.0, 0.0, 100.0, 100.0];
    let matrix = create_matrix(4, 1, &data);
    let nb_model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    let err = calculate_sgrna_pvalues(&matrix, &[2, 3], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));
}

#[test]
fn test_lfc_faithfulness() {
    let data = vec![100.0, 100.0, 50.0, 50.0];
    let matrix = create_matrix(4, 1, &data);
    let nb_model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    let res = calculate_sgrna_pvalues(&matrix, &[2, 3], &nb_model).unwrap();

    let expected_lfc = (51.0_f64).log2() - (101.0_f64).log2();
    assert!((res[0].lfc - expected_lfc).abs() < 1e-9);
}

#[test]
fn test_fractional_counts() {
    let data = vec![1.5, 2.5, 0.5, 0.5];
    let matrix = create_matrix(4, 1, &data);
    let nb_model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    let res = calculate_sgrna_pvalues(&matrix, &[2, 3], &nb_model).unwrap();

    assert_eq!(res[0].control_mean, 2.0);
    assert_eq!(res[0].treatment_mean, 0.5);

    let expected_lfc = (1.5_f64).log2() - (3.0_f64).log2();
    assert!((res[0].lfc - expected_lfc).abs() < 1e-9);
    assert!(res[0].p_lower < 0.5);
}

#[test]
fn test_treatment_indices_validation() {
    let matrix = create_matrix(4, 5, &[100.0; 20]);
    let nb_model = fit_meanvar(&matrix, &[0, 1]).unwrap();

    let err = calculate_sgrna_pvalues(&matrix, &[], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));

    let err = calculate_sgrna_pvalues(&matrix, &[4], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));

    let err = calculate_sgrna_pvalues(&matrix, &[5], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));

    let err = calculate_sgrna_pvalues(&matrix, &[2, 2], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));

    let res = calculate_sgrna_pvalues(&matrix, &[2], &nb_model).unwrap();
    assert_eq!(res.len(), 5);
}

#[test]
fn test_mismatched_nb_model() {
    let matrix1 = create_matrix(4, 5, &[100.0; 20]);
    let nb_model = fit_meanvar(&matrix1, &[0, 1]).unwrap();

    let matrix2 = create_matrix(4, 4, &[100.0; 16]);
    let err = calculate_sgrna_pvalues(&matrix2, &[2, 3], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));
}

#[test]
fn test_invalid_treatment_values() {
    let mut matrix = create_matrix(3, 1, &[10.0, 10.0, 10.0]);

    matrix.set(2, 0, -10.0);
    let nb_model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    let err = calculate_sgrna_pvalues(&matrix, &[2], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));

    matrix.set(2, 0, f64::NAN);
    let err = calculate_sgrna_pvalues(&matrix, &[2], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));

    matrix.set(2, 0, f64::INFINITY);
    let err = calculate_sgrna_pvalues(&matrix, &[2], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));

    matrix.set(2, 0, f64::NEG_INFINITY);
    let err = calculate_sgrna_pvalues(&matrix, &[2], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));
}

#[test]
fn test_rank_boundaries_and_determinism() {
    let data = vec![
        100.0, 100.0, 50.0, 100.0, 100.0, 50.0, 100.0, 100.0, 50.0, 100.0, 100.0, 50.0,
    ];
    let matrix = create_matrix(3, 4, &data);
    let nb_model = fit_meanvar(&matrix, &[0, 1]).unwrap();

    let res1 = calculate_sgrna_pvalues(&matrix, &[2], &nb_model).unwrap();
    let res2 = calculate_sgrna_pvalues(&matrix, &[2], &nb_model).unwrap();

    assert_eq!(res1, res2);

    for r in res1 {
        assert!(r.rank_negative > 0.0 && r.rank_negative <= 1.0);
        assert!(r.rank_positive > 0.0 && r.rank_positive <= 1.0);
    }
}

#[test]
fn test_empty_matrix() {
    let matrix = create_matrix(4, 0, &[]);
    let dummy_matrix = create_matrix(4, 5, &[100.0; 20]);
    let nb_model = fit_meanvar(&dummy_matrix, &[0, 1]).unwrap();
    let err = calculate_sgrna_pvalues(&matrix, &[2, 3], &nb_model).unwrap_err();
    assert!(matches!(err, MageckError::Math(_)));
}

#[test]
fn test_unequal_replicates_and_zeros() {
    let data = vec![40.0, 60.0, 40.0, 60.0, 0.0, 20.0, 0.0, 20.0];
    let matrix = create_matrix(4, 2, &data);
    let nb_model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    let res = calculate_sgrna_pvalues(&matrix, &[2, 3], &nb_model).unwrap();

    assert_eq!(res[0].control_mean, 50.0);
    assert_eq!(res[0].treatment_mean, 50.0);
    assert_eq!(res[1].control_mean, 10.0);
    assert_eq!(res[1].treatment_mean, 10.0);
}

#[test]
fn test_probability_relationships_and_lfc_independence() {
    let data = vec![100.0, 100.0, 50.0, 50.0];
    let matrix = create_matrix(4, 1, &data);
    let nb_model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    let res = calculate_sgrna_pvalues(&matrix, &[2, 3], &nb_model).unwrap();

    let p_lower = res[0].p_lower;
    let p_higher = res[0].p_higher;

    assert!((p_lower + p_higher - 1.0).abs() < 1e-9);

    let expected_p2 = 2.0 * p_lower.min(p_higher);
    assert!((res[0].p_twosided - expected_p2).abs() < 1e-9);

    assert!(res[0].lfc < 0.0);
    assert!(res[0].p_lower > 0.0);
}

#[test]
fn test_qnorm_table_properties() {
    let table = get_qnorm_table().unwrap();
    assert_eq!(table.logp_vec.len(), 24000);
    assert_eq!(table.norm_qnorm_vec.len(), 24000);
    assert_eq!(table.qmin, -24.0);
    assert_eq!(table.qmax, -0.001);
    assert!((table.qstep - 0.001).abs() < 1e-9);
}

#[test]
fn test_qnorm_exact_points() {
    let v_first = get_qnorm(-24.0).unwrap().unwrap();
    assert!((v_first - -6.5092982359543).abs() < 1e-12);

    // From pinned reference commit, qnorm_table.txt for LOGP -12.000
    let v_mid = get_qnorm(-12.0).unwrap().unwrap();
    let expected_mid = -4.37240686011874;
    assert!((v_mid - expected_mid).abs() < 1e-12);

    let v_last = get_qnorm(-0.001).unwrap().unwrap();
    assert!((v_last - 3.09038078691668).abs() < 1e-12);
}

#[test]
fn test_qnorm_interpolation() {
    let table = get_qnorm_table().unwrap();

    let v1 = table.norm_qnorm_vec[11999]; // -12.001
    let v2 = table.norm_qnorm_vec[12000]; // -12.000
    let expected_half = v1 * 0.5 + v2 * 0.5;
    let v_half = get_qnorm(-12.0005).unwrap().unwrap();
    assert!((v_half - expected_half).abs() < 1e-12);

    let expected_non_half = v1 * 0.7 + v2 * 0.3;
    let v_non_half = get_qnorm(-12.0007).unwrap().unwrap();
    assert!((v_non_half - expected_non_half).abs() < 1e-12);
}

#[test]
fn test_qnorm_boundaries_and_fallback() {
    // 10. Below qmin
    assert_eq!(get_qnorm(-24.0001).unwrap(), None);

    // 11. Inclusive qmin
    assert!(get_qnorm(-24.0).unwrap().is_some());

    // 12. Inclusive qmax
    assert!(get_qnorm(-0.001).unwrap().is_some());

    // 13. Above qmax
    assert_eq!(get_qnorm(-0.0009).unwrap(), None);

    // 14. Below qmin depletion fallback
    let data = vec![1e8, 1e8, 0.0, 0.0, 1000.0, 1000.0, 250.0, 250.0];
    let matrix = create_matrix(4, 2, &data);
    let nb_model = fit_meanvar(&matrix, &[0, 1]).unwrap();
    let res = calculate_sgrna_pvalues(&matrix, &[2, 3], &nb_model).unwrap();

    // Reconstruct the global minimum replacement rule:
    // 1. Collect every positive p_lower
    // 2. Calculate global_min_replacement = min(positive p_lower * 0.1)
    let global_min_replacement = res
        .iter()
        .filter_map(|r| {
            if r.p_lower > 0.0 {
                Some(r.p_lower * 0.1)
            } else {
                None
            }
        })
        .fold(f64::INFINITY, f64::min);

    // 3. For the selected guide (guide 0), calculate lookup_probability
    let selected_p_lower = res[0].p_lower;
    let lookup_probability = selected_p_lower.max(global_min_replacement);

    // 4. Calculate ln_p = ln(lookup_probability)
    let ln_p = lookup_probability.ln();

    // 5. Assert ln_p < qmin
    assert!(ln_p < -24.0);

    // 6. Assert get_qnorm(ln_p) == None
    assert_eq!(get_qnorm(ln_p).unwrap(), None);

    // 7. Calculate theta_0
    let var0 = res[0].adjusted_variance;
    let theta_0 = (res[0].treatment_mean - res[0].control_mean) / var0.sqrt();

    // 8. Assert negative score equals theta_0
    assert!((res[0].score_negative - theta_0).abs() < 1e-9);

    // 15. Enrichment classification
    // Enrichment uses theta_0 directly because treatment_mean >= control_mean, not because of QNorm fallback!
    let data2 = vec![100.0, 100.0, 1e8, 1e8];
    let matrix2 = create_matrix(4, 1, &data2);
    let nb_model2 = fit_meanvar(&matrix2, &[0, 1]).unwrap();
    let res2 = calculate_sgrna_pvalues(&matrix2, &[2, 3], &nb_model2).unwrap();
    let var0_2 = res2[0].adjusted_variance;
    let theta_0_2 = (res2[0].treatment_mean - res2[0].control_mean) / var0_2.sqrt();
    assert!((res2[0].score_negative - theta_0_2).abs() < 1e-9);
}

#[test]
fn test_proof_of_interpolation_over_inverse_cdf() {
    let lp = -0.0015;
    let norm = Normal::new(0.0, 1.0).unwrap();

    let table_val = get_qnorm(lp).unwrap().unwrap();
    let cdf_val = norm.inverse_cdf(lp.exp());

    // We expect the interpolated table_val to match interpolation of the exact table points.
    // -0.0015 is halfway between -0.002 and -0.001.
    // Let's assert it deviates by more than a meaningful threshold (0.01) from the continuous cdf.
    assert!((table_val - cdf_val).abs() > 0.01);

    // Demonstrate it strictly uses the literal table endpoints:
    let expected_exact_interp = (2.87847711980085 + 3.09038078691668) / 2.0;
    assert!((table_val - expected_exact_interp).abs() < 1e-12);
    // Continuous inverse CDF should fail the strict tolerance
    assert!((cdf_val - expected_exact_interp).abs() > 1e-12);
}

#[test]
fn test_parse_qnorm_table_errors() {
    let err1 = parse_qnorm_table("logp qnorm\n-24.0").unwrap_err();
    assert!(err1.contains("Invalid header fields"));

    let err2 = parse_qnorm_table("LOGP NORM_QNORM\n-24.0").unwrap_err();
    assert!(err2.contains("Missing qnorm column"));

    let err_missing = parse_qnorm_table("LOGP NORM_QNORM\n -6.5").unwrap_err();
    assert!(err_missing.contains("Missing qnorm column"));

    let err_extra = parse_qnorm_table("LOGP NORM_QNORM\n-24.0 -6.5 EXTRA").unwrap_err();
    assert!(err_extra.contains("Extra column found"));

    let err3 = parse_qnorm_table("LOGP NORM_QNORM\n-24.0 abc").unwrap_err();
    assert!(err3.contains("Invalid qnorm"));

    let err4 = parse_qnorm_table("LOGP NORM_QNORM\n-24.0 -6.5").unwrap_err();
    assert!(err4.contains("fewer than 2 data rows"));

    let err5 =
        parse_qnorm_table("LOGP NORM_QNORM\n-24.0 -6.5\n-23.0 -6.4\n-23.0 -6.4").unwrap_err();
    assert!(err5.contains("not strictly ascending"));

    let err6 =
        parse_qnorm_table("LOGP NORM_QNORM\n-24.0 -6.5\n-23.0 -6.4\n-21.0 -6.2").unwrap_err();
    assert!(err6.contains("Inconsistent qstep"));

    // Spacing deviation that would pass 1e-4 but fails 1e-12:
    let err_dev =
        parse_qnorm_table("LOGP NORM_QNORM\n-24.000 -6.5\n-23.999 -6.4\n-23.9979999 -6.3")
            .unwrap_err();
    assert!(err_dev.contains("Inconsistent qstep"));
}

#[test]
fn test_get_qnorm_infinity_nan() {
    let err_nan = get_qnorm(f64::NAN).unwrap_err();
    assert!(matches!(err_nan, MageckError::Math(_)));

    let err_pos_inf = get_qnorm(f64::INFINITY).unwrap_err();
    assert!(matches!(err_pos_inf, MageckError::Math(_)));

    let err_neg_inf = get_qnorm(f64::NEG_INFINITY).unwrap_err();
    assert!(matches!(err_neg_inf, MageckError::Math(_)));
}
