//! Deterministic numerical WLS reproducer; no biological calibration claim.
use mageckrs::data::CountMatrix;
use mageckrs::mageck_test::NbModel;

#[test]
fn pooled_wls_close_means_matches_independently_centered_arithmetic() {
    let mut control = Vec::new();
    let mut treatment = Vec::new();
    for i in 0..8 {
        let mean = 1_000_000.0 + f64::from(i);
        let offset = 100_000.0 + 0.3 * f64::from(i);
        control.push(mean - offset);
        treatment.push(mean + offset);
    }
    control.extend(treatment);
    let counts = CountMatrix::from_vec(control, 2, 8).unwrap();
    let model = NbModel::fit_pooled_dispersion(&counts, 0, &[0, 1]).unwrap();
    let coefficients = model.regression_coefficients().unwrap();
    // Independent Python centered WLS on identical f64 moments/log inputs.
    // Same-input coefficient tolerance retained from the P4F24B contract.
    assert!(
        (coefficients.slope - 4.0002080061066785).abs() <= 1e-9,
        "slope {}",
        coefficients.slope
    );
    assert!(
        (coefficients.intercept - -59.176988949056216).abs() <= 1e-9,
        "intercept {}",
        coefficients.intercept
    );
}

#[test]
fn frozen_wls_workflow_and_first_failure_are_thread_deterministic() {
    use std::{fs, path::Path, process::Command};
    let fixture = |name: &str| {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/campaign_c2")
            .join(name)
    };
    for control_null in [false, true] {
        let mut expected = None;
        for workers in [1, 2, 4, 8] {
            let dir = tempfile::tempdir().unwrap();
            let prefix = dir.path().join("result");
            let mut command = Command::new(env!("CARGO_BIN_EXE_mageckrs"));
            command.args(["test", "-c", "C", "-t", "T", "--norm-method", "none"]);
            command.arg("-k").arg(fixture("close-means.count.txt"));
            command.arg("-n").arg(&prefix);
            command.args(["--threads", &workers.to_string()]);
            if control_null {
                command.arg("--control-sgrna").arg(fixture("controls.txt"));
            }
            let output = command.output().unwrap();
            assert!(output.status.success(), "{output:?}");
            let files = ["gene_summary.txt", "sgrna_summary.txt"]
                .map(|suffix| fs::read(dir.path().join(format!("result.{suffix}"))).unwrap());
            if let Some(expected) = &expected {
                assert_eq!(&files, expected);
            } else {
                expected = Some(files);
            }
        }
    }
    let mut expected_failure = None;
    for workers in [1, 2, 4, 8] {
        let dir = tempfile::tempdir().unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
            .args(["test", "-c", "C", "-t", "T", "--norm-method", "none"])
            .arg("-k")
            .arg(fixture("failure.count.txt"))
            .arg("-n")
            .arg(dir.path().join("result"))
            .args(["--threads", &workers.to_string()])
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
        let failure = (output.status.code(), output.stdout, output.stderr);
        if let Some(expected) = &expected_failure {
            assert_eq!(&failure, expected);
        } else {
            expected_failure = Some(failure);
        }
    }
}
