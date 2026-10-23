//! Exact process-boundary scientific and failure parity for explicit CPU workers.
use std::{collections::BTreeMap, fs, path::Path, process::Command};

fn invoke(args: &[String]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(args)
        .env_remove("MAGECKRS_PROFILE")
        .output()
        .unwrap()
}
fn fixture(path: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(path)
        .display()
        .to_string()
}
fn files(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    fs::read_dir(dir)
        .unwrap()
        .map(|p| {
            let p = p.unwrap().path();
            (
                p.file_name().unwrap().to_str().unwrap().to_owned(),
                fs::read(p).unwrap(),
            )
        })
        .collect()
}
fn compare(args: Vec<String>) {
    let d = tempfile::tempdir().unwrap();
    let mut base = args;
    base.extend(["-n".into(), d.path().join("result").display().to_string()]);
    let serial = invoke(&base);
    assert!(serial.status.success(), "{serial:?}");
    let expected = files(d.path());
    assert!(!expected.is_empty());
    for threads in [1, 2, 4, 8, 16] {
        for _ in 0..3 {
            for p in fs::read_dir(d.path()).unwrap() {
                fs::remove_file(p.unwrap().path()).unwrap();
            }
            let mut threaded = base.clone();
            threaded.extend(["--threads".into(), threads.to_string()]);
            let actual = invoke(&threaded);
            assert_eq!(actual.status.code(), serial.status.code());
            assert_eq!(actual.stdout, serial.stdout);
            assert_eq!(actual.stderr, serial.stderr);
            assert_eq!(files(d.path()), expected, "threads={threads}");
        }
    }
}
#[test]
fn mle_scalar_training_and_cnv_artifacts_match_repeatedly() {
    for extra in [
        vec![],
        vec!["--genes-varmodeling", "1"],
        vec![
            "--genes-varmodeling",
            "1",
            "--cnv-norm",
            "tests/fixtures/m11a/cnv.tsv",
        ],
    ] {
        let mut args = vec![
            "mle".into(),
            "-k".into(),
            fixture("tests/fixtures/m11a/counts.tsv"),
            "-d".into(),
            fixture("tests/fixtures/m11a/design.tsv"),
            "--norm-method".into(),
            "none".into(),
        ];
        for arg in extra {
            args.push(if arg.starts_with("tests/") {
                fixture(arg)
            } else {
                arg.into()
            });
        }
        compare(args);
    }
}
#[test]
fn rra_ordinary_control_and_seeds_match_across_batch_boundaries() {
    for seed in ["0", "18446744073709551615"] {
        for control in [false, true] {
            let mut args = vec![
                "test".into(),
                "-k".into(),
                fixture("tests/fixtures/m11a/counts.tsv"),
                "-t".into(),
                "TreatA,TreatB".into(),
                "-c".into(),
                "Day0".into(),
                "--norm-method".into(),
                "none".into(),
                "--permutation".into(),
                "129".into(),
                "--seed".into(),
                seed.into(),
            ];
            if control {
                args.extend([
                    "--control-sgrna".into(),
                    fixture("tests/fixtures/m11a/controls.txt"),
                ]);
            }
            compare(args);
        }
    }
}
#[test]
fn combined_forwards_analysis_threads_and_preserves_count_artifacts() {
    compare(vec![
        "run".into(),
        "-l".into(),
        fixture("examples/campaign-a/library.tsv"),
        "-f".into(),
        format!(
            "{},{}",
            fixture("examples/campaign-a/sample-0.fastq"),
            fixture("examples/campaign-a/sample-1.fastq")
        ),
        fixture("examples/campaign-a/sample-2.fastq"),
        fixture("examples/campaign-a/sample-3.fastq"),
        "--sample-label".into(),
        "C,T1,T2".into(),
        "-t".into(),
        "T1,T2".into(),
        "-c".into(),
        "C".into(),
        "--sgrna-len".into(),
        "4".into(),
        "--norm-method".into(),
        "none".into(),
        "--permutation".into(),
        "129".into(),
    ]);
}
#[test]
fn zero_threads_reject_before_missing_inputs_for_all_supported_commands() {
    for args in [
        vec!["mle", "-k", "missing", "-d", "missing"],
        vec!["test", "-k", "missing", "-t", "T", "-c", "C"],
        vec![
            "run", "-l", "missing", "-f", "missing", "-t", "T", "-c", "C",
        ],
    ] {
        let mut args: Vec<String> = args.into_iter().map(str::to_owned).collect();
        args.extend(["-n".into(), "unused".into(), "--threads".into(), "0".into()]);
        let result = invoke(&args);
        assert_eq!(result.status.code(), Some(1));
        assert_eq!(
            String::from_utf8(result.stderr).unwrap(),
            "Configuration error: --threads must be greater than zero\n"
        );
        assert!(result.stdout.is_empty());
    }
}
#[test]
fn analytical_and_publication_failures_match_serial_and_publish_nothing_new() {
    let d = tempfile::tempdir().unwrap();
    let count = d.path().join("counts.tsv");
    fs::write(
        &count,
        "sgRNA\tGene\tDay0\tTreatA\tTreatB\na\tA\t10\t11.5\t12\nb\tB\t20\t21.5\t22\n",
    )
    .unwrap();
    let prefix = d.path().join("result").display().to_string();
    let cases = [
        vec![
            "mle".into(),
            "-k".into(),
            count.display().to_string(),
            "-d".into(),
            fixture("tests/fixtures/m11a/design.tsv"),
        ],
        vec![
            "test".into(),
            "-k".into(),
            fixture("tests/fixtures/m11a/counts.tsv"),
            "-t".into(),
            "missing".into(),
            "-c".into(),
            "Day0".into(),
        ],
    ];
    for mut args in cases {
        args.extend(["-n".into(), prefix.clone()]);
        let expected = invoke(&args);
        assert!(!expected.status.success());
        for threads in [1, 2, 4, 8] {
            for _ in 0..3 {
                let mut args = args.clone();
                args.extend(["--threads".into(), threads.to_string()]);
                let actual = invoke(&args);
                assert_eq!(actual.status.code(), expected.status.code());
                assert_eq!(actual.stdout, expected.stdout);
                assert_eq!(actual.stderr, expected.stderr);
                assert_eq!(fs::read_dir(d.path()).unwrap().count(), 1);
            }
        }
    }

    // Both independently fitted genes fail after the worker pool is constructed.
    // The first source error and diagnostic-only publication must match.
    fs::write(&count, "sgRNA\tGene\tDay0\tTreatA\tTreatB\na1\tFIRST\t1e308\t1e308\t1e308\na2\tFIRST\t1e308\t1e308\t1e308\nb1\tSECOND\t1e300\t1e300\t1e300\nb2\tSECOND\t1e300\t1e300\t1e300\n").unwrap();
    let args = vec![
        "mle".into(),
        "-k".into(),
        count.display().to_string(),
        "-d".into(),
        fixture("tests/fixtures/m11a/design.tsv"),
        "--norm-method".into(),
        "none".into(),
        "-n".into(),
        prefix.clone(),
    ];
    let expected = invoke(&args);
    assert!(!expected.status.success());
    assert!(String::from_utf8_lossy(&expected.stderr).contains("gene 'FIRST'"));
    let diagnostic_path = d.path().join("result.mle_diagnostics.json");
    let expected_diagnostic = fs::read(&diagnostic_path).unwrap();
    fs::remove_file(&diagnostic_path).unwrap();
    for threads in [1, 2, 4, 8] {
        for _ in 0..3 {
            let mut args = args.clone();
            args.extend(["--threads".into(), threads.to_string()]);
            let actual = invoke(&args);
            assert_eq!(actual.stderr, expected.stderr);
            assert_eq!(actual.status.code(), expected.status.code());
            assert_eq!(fs::read_dir(d.path()).unwrap().count(), 2);
            assert_eq!(fs::read(&diagnostic_path).unwrap(), expected_diagnostic);
            fs::remove_file(&diagnostic_path).unwrap();
        }
    }
    let sentinel = d.path().join("result.gene_summary.txt");
    fs::write(&sentinel, "owner data").unwrap();
    let args = vec![
        "mle".into(),
        "-k".into(),
        fixture("tests/fixtures/m11a/counts.tsv"),
        "-d".into(),
        fixture("tests/fixtures/m11a/design.tsv"),
        "-n".into(),
        prefix,
    ];
    let expected = invoke(&args);
    assert!(!expected.status.success());
    for threads in [1, 2, 4, 8] {
        let mut args = args.clone();
        args.extend(["--threads".into(), threads.to_string()]);
        let actual = invoke(&args);
        assert_eq!(actual.stderr, expected.stderr);
        assert_eq!(actual.status.code(), expected.status.code());
        assert_eq!(fs::read_to_string(&sentinel).unwrap(), "owner data");
        assert_eq!(fs::read_dir(d.path()).unwrap().count(), 2);
    }
}
