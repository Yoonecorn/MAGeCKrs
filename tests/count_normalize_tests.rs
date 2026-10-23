use mageckrs::cli::NormMethod;
use mageckrs::count::normalize::normalize;
use mageckrs::data::count_matrix::CountMatrix;
use mageckrs::error::MageckError;

fn is_normalization_error(res: &Result<CountMatrix, MageckError>) -> bool {
    matches!(res, Err(MageckError::Config(_)))
}

fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

fn assert_matrix_eq(a: &CountMatrix, b: &CountMatrix, tol: f64) {
    assert_eq!(a.n_samples(), b.n_samples());
    assert_eq!(a.n_guides(), b.n_guides());
    for s in 0..a.n_samples() {
        for g in 0..a.n_guides() {
            assert!(
                approx_eq(a.entry(s, g), b.entry(s, g), tol),
                "Mismatch at ({}, {}): {} vs {}",
                s,
                g,
                a.entry(s, g),
                b.entry(s, g)
            );
        }
    }
}

#[test]
fn test_empty_matrix() {
    let mat = CountMatrix::new(0, 0);
    let res = normalize(&mat, NormMethod::Median, None).unwrap();
    assert_eq!(res.n_samples(), 0);
    assert_eq!(res.n_guides(), 0);
}

#[test]
fn test_zero_samples() {
    let mat = CountMatrix::new(0, 10);
    let res = normalize(&mat, NormMethod::Median, None).unwrap();
    assert_eq!(res.n_samples(), 0);
    assert_eq!(res.n_guides(), 10);
}

#[test]
fn test_zero_guides() {
    let mat = CountMatrix::new(5, 0);
    let res = normalize(&mat, NormMethod::Median, None).unwrap();
    assert_eq!(res.n_samples(), 5);
    assert_eq!(res.n_guides(), 0);
}

#[test]
fn test_one_sample_one_guide() {
    let mut mat = CountMatrix::new(1, 1);
    mat.set(0, 0, 10.0);
    let res = normalize(&mat, NormMethod::Median, None).unwrap();
    // GM = 10, Ratio = 1, SF = 1, Output = 10
    assert!(approx_eq(res.entry(0, 0), 10.0, 1e-9));
}

#[test]
fn test_all_zero_matrix() {
    let mat = CountMatrix::new(2, 3); // All zeros
    let res = normalize(&mat, NormMethod::Median, None);
    assert!(is_normalization_error(&res)); // Every guide excluded
}

#[test]
fn test_negative_count_rejection() {
    let mut mat = CountMatrix::new(1, 1);
    mat.set(0, 0, -1.0);
    let res = normalize(&mat, NormMethod::Median, None);
    assert!(is_normalization_error(&res));
}

#[test]
fn test_nan_rejection() {
    let mut mat = CountMatrix::new(1, 1);
    mat.set(0, 0, f64::NAN);
    let res = normalize(&mat, NormMethod::Median, None);
    assert!(is_normalization_error(&res));
}

#[test]
fn test_inf_rejection() {
    let mut mat = CountMatrix::new(1, 1);
    mat.set(0, 0, f64::INFINITY);
    let res = normalize(&mat, NormMethod::Median, None);
    assert!(is_normalization_error(&res));

    mat.set(0, 0, f64::NEG_INFINITY);
    let res = normalize(&mat, NormMethod::Median, None);
    assert!(is_normalization_error(&res));
}

#[test]
fn test_normalization_none() {
    let mut mat = CountMatrix::new(2, 2);
    mat.set(0, 0, 1.0);
    mat.set(0, 1, 2.0);
    mat.set(1, 0, 3.0);
    mat.set(1, 1, 4.0);
    let res = normalize(&mat, NormMethod::None, None).unwrap();
    assert_matrix_eq(&mat, &res, 1e-9);
}

#[test]
fn test_total_normalization() {
    // 2 samples, 2 guides
    // s0: 10, 30 -> total 40
    // s1: 40, 40 -> total 80
    // mean total = 60
    // SF0 = 40/60 = 2/3
    // SF1 = 80/60 = 4/3
    let mut mat = CountMatrix::new(2, 2);
    mat.set(0, 0, 10.0);
    mat.set(0, 1, 30.0);
    mat.set(1, 0, 40.0);
    mat.set(1, 1, 40.0);

    let res = normalize(&mat, NormMethod::Total, None).unwrap();
    // expected:
    // s0: 10 / (2/3) = 15, 30 / (2/3) = 45 -> total 60
    // s1: 40 / (4/3) = 30, 40 / (4/3) = 30 -> total 60
    assert!(approx_eq(res.entry(0, 0), 15.0, 1e-9));
    assert!(approx_eq(res.entry(0, 1), 45.0, 1e-9));
    assert!(approx_eq(res.entry(1, 0), 30.0, 1e-9));
    assert!(approx_eq(res.entry(1, 1), 30.0, 1e-9));
}

#[test]
fn test_median_ratio_simple() {
    // 2 samples, 3 guides.
    // G0: 2, 8 -> GM = 4. Ratios: s0=0.5, s1=2.0
    // G1: 4, 16 -> GM = 8. Ratios: s0=0.5, s1=2.0
    // G2: 8, 32 -> GM = 16. Ratios: s0=0.5, s1=2.0
    // SF0 = 0.5, SF1 = 2.0
    let mut mat = CountMatrix::new(2, 3);
    mat.set(0, 0, 2.0);
    mat.set(0, 1, 4.0);
    mat.set(0, 2, 8.0);
    mat.set(1, 0, 8.0);
    mat.set(1, 1, 16.0);
    mat.set(1, 2, 32.0);

    let res = normalize(&mat, NormMethod::Median, None).unwrap();

    // sf0 = 0.5. s0 values: 4, 8, 16
    assert!(approx_eq(res.entry(0, 0), 4.0, 1e-9));
    assert!(approx_eq(res.entry(0, 1), 8.0, 1e-9));
    assert!(approx_eq(res.entry(0, 2), 16.0, 1e-9));

    // sf1 = 2.0. s1 values: 4, 8, 16
    assert!(approx_eq(res.entry(1, 0), 4.0, 1e-9));
    assert!(approx_eq(res.entry(1, 1), 8.0, 1e-9));
    assert!(approx_eq(res.entry(1, 2), 16.0, 1e-9));
}

#[test]
fn test_median_odd_and_even_ratios() {
    // 1 sample, 4 guides (even)
    // G0=1, G1=2, G2=3, G3=4 -> GM=self -> Ratios=1,1,1,1 -> SF=1.0

    // 2 samples, 4 guides
    // G0: 1, 2 (GM=sqrt(2)) -> ratios s0: 1/sqrt(2), s1: 2/sqrt(2)
    // let's do a case where ratios are distinct for even median
    let mut mat = CountMatrix::new(2, 4);
    // G0: s0=10, s1=10 -> GM=10. s0=1, s1=1
    // G1: s0=20, s1=20 -> GM=20. s0=1, s1=1
    // G2: s0=30, s1=30 -> GM=30. s0=1, s1=1
    // G3: s0=40, s1=40 -> GM=40. s0=1, s1=1
    // Ratios are all 1. Let's make them differ:

    mat.set(0, 0, 1.0);
    mat.set(1, 0, 10.0); // Ratios s0: 1/sqrt(10), s1: 10/sqrt(10)
    mat.set(0, 1, 2.0);
    mat.set(1, 1, 10.0); // Ratios s0: 2/sqrt(20), s1: 10/sqrt(20)
    mat.set(0, 2, 4.0);
    mat.set(1, 2, 10.0); // Ratios s0: 4/sqrt(40), s1: 10/sqrt(40)
    mat.set(0, 3, 8.0);
    mat.set(1, 3, 10.0); // Ratios s0: 8/sqrt(80), s1: 10/sqrt(80)

    // for s0, ratios are 1/3.16, 2/4.47, 4/6.32, 8/8.94
    // 0.316, 0.447, 0.632, 0.894
    // median should be (0.447 + 0.632)/2
    let res = normalize(&mat, NormMethod::Median, None).unwrap();
    let sf0_expected = ((2.0 / 20.0f64.sqrt()) + (4.0 / 40.0f64.sqrt())) / 2.0;

    assert!(approx_eq(res.entry(0, 0), 1.0 / sf0_expected, 1e-9));
}

#[test]
fn test_zero_exclusion_rule() {
    let mut mat = CountMatrix::new(2, 3);
    mat.set(0, 0, 10.0);
    mat.set(1, 0, 20.0); // valid
    mat.set(0, 1, 0.0);
    mat.set(1, 1, 10.0); // s0 has zero -> exclude G1
    mat.set(0, 2, 10.0);
    mat.set(1, 2, 0.0); // s1 has zero -> exclude G2

    // Only G0 is valid. GM0 = sqrt(200) = 14.1421356237
    // Ratios for s0: 10 / 14.1421356237 = 0.70710678
    // Ratios for s1: 20 / 14.1421356237 = 1.41421356

    let res = normalize(&mat, NormMethod::Median, None).unwrap();
    let gm0 = 200.0f64.sqrt();
    let sf0 = 10.0 / gm0;
    let sf1 = 20.0 / gm0;

    assert!(approx_eq(res.entry(0, 0), 10.0 / sf0, 1e-9));
    // The excluded guides still get normalized by the sample's size factor
    assert!(approx_eq(res.entry(0, 1), 0.0 / sf0, 1e-9));
    assert!(approx_eq(res.entry(1, 1), 10.0 / sf1, 1e-9));
}

#[test]
fn test_control_normalization() {
    let mut mat = CountMatrix::new(2, 4);
    // G0, G1: control guides. G2, G3: test guides.
    mat.set(0, 0, 2.0);
    mat.set(1, 0, 8.0); // gm=4. r0=0.5, r1=2.0
    mat.set(0, 1, 4.0);
    mat.set(1, 1, 16.0); // gm=8. r0=0.5, r1=2.0

    mat.set(0, 2, 100.0);
    mat.set(1, 2, 100.0); // Should be ignored in SF calc
    mat.set(0, 3, 0.0);
    mat.set(1, 3, 10.0); // Should be ignored

    let controls = vec![0, 1];
    let res = normalize(&mat, NormMethod::Control, Some(&controls)).unwrap();

    // SF0 = 0.5, SF1 = 2.0
    assert!(approx_eq(res.entry(0, 2), 100.0 / 0.5, 1e-9));
    assert!(approx_eq(res.entry(1, 2), 100.0 / 2.0, 1e-9));
}

#[test]
fn test_all_guides_excluded() {
    let mut mat = CountMatrix::new(2, 2);
    mat.set(0, 0, 0.0);
    mat.set(1, 0, 10.0);
    mat.set(0, 1, 10.0);
    mat.set(1, 1, 0.0);
    let res = normalize(&mat, NormMethod::Median, None);
    assert!(is_normalization_error(&res));
}

#[test]
fn test_large_counts_no_overflow() {
    let mut mat = CountMatrix::new(2, 1);
    mat.set(0, 0, 1e150);
    mat.set(1, 0, 1e150);
    let res = normalize(&mat, NormMethod::Median, None).unwrap();
    // geometric mean is 1e150. ratio is 1. sf is 1. output 1e150
    assert!(approx_eq(res.entry(0, 0), 1e150, 1e141));
}

#[test]
fn test_original_matrix_unchanged() {
    let mut mat = CountMatrix::new(2, 2);
    mat.set(0, 0, 10.0);
    mat.set(0, 1, 20.0);
    mat.set(1, 0, 30.0);
    mat.set(1, 1, 40.0);
    let _res = normalize(&mat, NormMethod::Total, None).unwrap();

    assert_eq!(mat.entry(0, 0), 10.0);
    assert_eq!(mat.entry(1, 1), 40.0);
}

#[test]
fn test_sample_and_guide_order_preservation() {
    let mut mat = CountMatrix::new(3, 3);
    for s in 0..3 {
        for g in 0..3 {
            mat.set(s, g, (s * 10 + g + 1) as f64);
        }
    }
    let res = normalize(&mat, NormMethod::None, None).unwrap();
    for s in 0..3 {
        for g in 0..3 {
            assert_eq!(res.entry(s, g), (s * 10 + g + 1) as f64);
        }
    }
}

#[test]
fn test_no_premature_rounding() {
    let mut mat = CountMatrix::new(1, 1);
    mat.set(0, 0, 10.0);
    // SF will be 1.0. If it were rounded it would be 10.0.
    // Let's create a case where SF is a fraction.
    let mut mat2 = CountMatrix::new(2, 1);
    mat2.set(0, 0, 1.0);
    mat2.set(1, 0, 3.0);
    let res = normalize(&mat2, NormMethod::Median, None).unwrap();
    // gm = sqrt(3).
    // s0 sf = 1/sqrt(3).
    // out s0 = 1 / (1/sqrt(3)) = sqrt(3) = 1.7320508075...
    assert!(res.entry(0, 0) != 1.7320); // Not rounded to 4 decimals
    assert!(approx_eq(res.entry(0, 0), 3.0f64.sqrt(), 1e-12));
}

#[test]
fn test_transpose_fail_case() {
    // If guide and sample axes are transposed, the size factors will differ.
    // 2 samples, 3 guides
    let mut mat = CountMatrix::new(2, 3);
    mat.set(0, 0, 2.0);
    mat.set(0, 1, 4.0);
    mat.set(0, 2, 8.0); // sum 14
    mat.set(1, 0, 8.0);
    mat.set(1, 1, 16.0);
    mat.set(1, 2, 32.0); // sum 56

    // total normalization sf0 = 14 / ((14+56)/2) = 14/35 = 0.4
    // sf1 = 56/35 = 1.6
    // out s0,g0 = 2/0.4 = 5
    let res = normalize(&mat, NormMethod::Total, None).unwrap();
    assert!(approx_eq(res.entry(0, 0), 5.0, 1e-9));
    assert!(approx_eq(res.entry(1, 0), 5.0, 1e-9));
}

#[test]
fn test_divide_not_multiply() {
    let mut mat = CountMatrix::new(2, 1);
    mat.set(0, 0, 10.0);
    mat.set(1, 0, 40.0);
    // Total sum = 50. Mean total = 25.
    // sf0 = 10 / 25 = 0.4
    // If we mistakenly multiply: 10 * 0.4 = 4.
    // If we divide: 10 / 0.4 = 25.
    let res = normalize(&mat, NormMethod::Total, None).unwrap();
    assert!(approx_eq(res.entry(0, 0), 25.0, 1e-9));
}

#[test]
fn test_zero_total_sample_with_positive() {
    let mut mat = CountMatrix::new(2, 1);
    mat.set(0, 0, 0.0);
    mat.set(1, 0, 10.0);
    let res = normalize(&mat, NormMethod::Total, None);
    assert!(is_normalization_error(&res));
}

#[test]
fn test_all_zero_totals() {
    let mut mat = CountMatrix::new(2, 1);
    mat.set(0, 0, 0.0);
    mat.set(1, 0, 0.0);
    let res = normalize(&mat, NormMethod::Total, None);
    assert!(is_normalization_error(&res));
}

#[test]
fn test_sum_overflow_total() {
    let mut mat = CountMatrix::new(1, 2);
    mat.set(0, 0, f64::MAX);
    mat.set(0, 1, f64::MAX);
    let res = normalize(&mat, NormMethod::Total, None);
    assert!(is_normalization_error(&res));
}

#[test]
fn test_empty_control_selection() {
    let mut mat = CountMatrix::new(1, 1);
    mat.set(0, 0, 10.0);
    let controls: Vec<usize> = vec![];
    let res = normalize(&mat, NormMethod::Control, Some(&controls));
    assert!(is_normalization_error(&res));
}

#[test]
fn test_duplicate_control_indices() {
    let mut mat = CountMatrix::new(2, 3);
    mat.set(0, 0, 2.0);
    mat.set(1, 0, 8.0);
    // Use guide 0 twice
    let controls = vec![0, 0];
    let res = normalize(&mat, NormMethod::Control, Some(&controls)).unwrap();
    // GM0 = 4. ratios: s0=0.5, s1=2.0
    // median of [0.5, 0.5] is 0.5.
    // If duplicates are allowed, it weights the median, which is correct behavior for the current impl.
    assert!(approx_eq(res.entry(0, 0), 4.0, 1e-9));
}

#[test]
fn test_out_of_range_control_index() {
    let mut mat = CountMatrix::new(1, 1);
    mat.set(0, 0, 10.0);
    let controls = vec![1]; // index out of bounds
    let res = normalize(&mat, NormMethod::Control, Some(&controls));
    assert!(is_normalization_error(&res));
}

#[test]
fn test_all_control_guides_excluded() {
    let mut mat = CountMatrix::new(2, 2);
    mat.set(0, 0, 0.0); // excluded
    mat.set(1, 0, 10.0);
    let controls = vec![0];
    let res = normalize(&mat, NormMethod::Control, Some(&controls));
    assert!(is_normalization_error(&res));
}

#[test]
fn test_one_eligible_median_guide() {
    let mut mat = CountMatrix::new(2, 2);
    mat.set(0, 0, 2.0);
    mat.set(1, 0, 8.0);
    mat.set(0, 1, 0.0); // excluded
    mat.set(1, 1, 10.0);
    let res = normalize(&mat, NormMethod::Median, None).unwrap();
    // GM0 = 4. SF0 = 0.5, SF1 = 2.0
    assert!(approx_eq(res.entry(0, 0), 4.0, 1e-9));
    assert!(approx_eq(res.entry(1, 1), 5.0, 1e-9));
}

#[test]
fn test_nonfinite_computed_size_factor() {
    let mut mat = CountMatrix::new(2, 1);
    mat.set(0, 0, f64::MIN_POSITIVE);
    mat.set(1, 0, 1.0);
    // subnormal or zero geometric mean handled?
    // gm = sqrt(1e-308) which is ~ 1e-154.
    // wait, what if we use something that underflows exp()?
    // sum_log / 2.0 = ln(min_pos) / 2.0. ln(min_pos) is -744. /2 = -372.
    // exp(-372) underflows to 0.0!
    // gm will be 0.0, which means it will be excluded!
    // WAIT, f64::exp(-372) is ~1e-161, which is fully representable and normal!
    let res = normalize(&mat, NormMethod::Median, None).unwrap();
    assert!(res.entry(0, 0) > 0.0);
    assert!(res.entry(1, 0) > 0.0);
}

#[test]
fn test_none_mode_invalid_value_policy() {
    let mut mat = CountMatrix::new(1, 1);
    mat.set(0, 0, f64::NAN);
    let res = normalize(&mat, NormMethod::None, None);
    assert!(is_normalization_error(&res)); // Should be rejected before None normalization
}

#[test]
fn test_negative_zero_behavior() {
    let mut mat = CountMatrix::new(1, 1);
    mat.set(0, 0, -0.0); // valid (equal to 0.0)
    let res = normalize(&mat, NormMethod::Median, None);
    // since it's 0, it's excluded. so all guides excluded
    assert!(is_normalization_error(&res));
}

#[test]
fn test_subnormal_positive_values() {
    let mut mat = CountMatrix::new(1, 1);
    mat.set(0, 0, f64::from_bits(1)); // smallest subnormal
    let res = normalize(&mat, NormMethod::Median, None).unwrap();
    // Rust's f64::exp(-744.4) accurately produces the subnormal 5e-324.
    // GM is correct. Ratio is 1.0. SF is 1.0.
    // It works perfectly!
    assert_eq!(res.entry(0, 0), f64::from_bits(1));
}

#[test]
fn test_original_matrix_unchanged_after_errors() {
    let mut mat = CountMatrix::new(1, 1);
    mat.set(0, 0, 10.0);
    let controls = vec![1]; // error
    let _res = normalize(&mat, NormMethod::Control, Some(&controls));
    assert_eq!(mat.entry(0, 0), 10.0);
}

#[test]
fn test_infinite_size_factor_ratio() {
    let mut mat2 = CountMatrix::new(2, 2);
    mat2.set(0, 0, f64::MAX);
    mat2.set(1, 0, f64::from_bits(1));
    mat2.set(0, 1, f64::MAX);
    mat2.set(1, 1, f64::from_bits(1));
    let res2 = normalize(&mat2, NormMethod::Median, None);
    assert!(is_normalization_error(&res2));
}
