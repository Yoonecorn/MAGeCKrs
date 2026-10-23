//! Process-boundary scientific identity for opt-in dynamic workers.
use std::{collections::BTreeMap, fs, path::Path, process::Command};
fn fixture(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(name)
        .display()
        .to_string()
}
fn invoke(args: &[String]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(args)
        .env_remove("MAGECKRS_PROFILE")
        .output()
        .unwrap()
}
fn scientific(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| {
            let name = p.file_name().unwrap().to_str().unwrap();
            !(name.ends_with(".run_info.txt") || name.ends_with(".count_run_info.txt"))
        })
        .map(|p| {
            (
                p.file_name().unwrap().to_str().unwrap().into(),
                fs::read(p).unwrap(),
            )
        })
        .collect()
}
fn compare(mut args: Vec<String>, stages: &[&str]) {
    let directory = tempfile::tempdir().unwrap();
    args.extend([
        "-n".into(),
        directory.path().join("result").display().to_string(),
    ]);
    let expected = invoke(&args);
    assert!(expected.status.success(), "{expected:?}");
    let files = scientific(directory.path());
    assert!(!files.is_empty());
    for mode in [
        "1",
        "2",
        "4",
        "8",
        "18446744073709551615",
        "auto",
        "auto:1",
        "auto:2",
        "auto:8",
    ] {
        for _ in 0..2 {
            for e in fs::read_dir(directory.path()).unwrap() {
                fs::remove_file(e.unwrap().path()).unwrap();
            }
            let mut selected = args.clone();
            selected.extend(["--threads".into(), mode.into()]);
            let actual = invoke(&selected);
            assert_eq!(actual.status.code(), expected.status.code(), "{actual:?}");
            assert_eq!(actual.stdout, expected.stdout);
            assert_eq!(actual.stderr, expected.stderr);
            assert_eq!(scientific(directory.path()), files, "mode={mode}");
            if mode.starts_with("auto") {
                let metadata = fs::read_dir(directory.path())
                    .unwrap()
                    .map(|e| e.unwrap().path())
                    .filter(|p| {
                        p.file_name()
                            .unwrap()
                            .to_str()
                            .unwrap()
                            .contains("run_info")
                    })
                    .map(|p| fs::read_to_string(p).unwrap())
                    .collect::<Vec<_>>()
                    .concat();
                for stage in stages {
                    assert!(metadata.contains(&format!("threads.{stage}.policy\tauto-v1")));
                    let prefix = format!("threads.{stage}.workers\t");
                    let workers: usize = metadata
                        .lines()
                        .find_map(|line| line.strip_prefix(&prefix))
                        .unwrap()
                        .parse()
                        .unwrap();
                    assert!((1..=8).contains(&workers));
                    if mode == "auto:1" || *stage == "count.serial" {
                        assert_eq!(workers, 1);
                    }
                    if mode == "auto:2" {
                        assert!(workers <= 2);
                    }
                }
            }
        }
    }
}
#[test]
fn rra_ordinary_control_and_seed_extremes_are_identical() {
    for seed in ["0", "18446744073709551615"] {
        for controls in [false, true] {
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
            if controls {
                args.extend([
                    "--control-sgrna".into(),
                    fixture("tests/fixtures/m11a/controls.txt"),
                ]);
            }
            compare(args, &["rra.empirical"]);
        }
    }
}
#[test]
fn mle_scalar_training_and_cnv_are_identical() {
    for training in [false, true] {
        for cnv in [false, true] {
            let mut args = vec![
                "mle".into(),
                "-k".into(),
                fixture("tests/fixtures/m11a/counts.tsv"),
                "-d".into(),
                fixture("tests/fixtures/m11a/design.tsv"),
                "--norm-method".into(),
                "none".into(),
            ];
            if training {
                args.extend(["--genes-varmodeling".into(), "3".into()]);
            }
            if cnv {
                args.extend(["--cnv-norm".into(), fixture("tests/fixtures/m11a/cnv.tsv")]);
            }
            compare(
                args,
                if training {
                    &["mle.training", "mle.fit"]
                } else {
                    &["mle.fit"]
                },
            );
        }
    }
}
#[test]
fn count_and_combined_keep_counting_serial_and_identical() {
    for command in ["count", "run"] {
        let mut args = vec![
            command.into(),
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
            "--sgrna-len".into(),
            "4".into(),
            "--norm-method".into(),
            "none".into(),
        ];
        if command == "run" {
            args.extend([
                "-t".into(),
                "T1,T2".into(),
                "-c".into(),
                "C".into(),
                "--permutation".into(),
                "129".into(),
                "--seed".into(),
                "0".into(),
            ]);
        }
        compare(
            args,
            if command == "run" {
                &["count.serial", "rra.empirical"]
            } else {
                &["count.serial"]
            },
        );
    }
}
#[test]
fn failures_and_existing_outputs_match_and_publish_nothing() {
    let directory = tempfile::tempdir().unwrap();
    let counts = directory.path().join("counts.tsv");
    fs::write(&counts, "sgRNA\tGene\tDay0\tTreatA\tTreatB\na1\tFIRST\t1e308\t1e308\t1e308\na2\tFIRST\t1e308\t1e308\t1e308\nb1\tSECOND\t1e300\t1e300\t1e300\nb2\tSECOND\t1e300\t1e300\t1e300\n").unwrap();
    let prefix = directory.path().join("result").display().to_string();
    let base = vec![
        "mle".into(),
        "-k".into(),
        counts.display().to_string(),
        "-d".into(),
        fixture("tests/fixtures/m11a/design.tsv"),
        "--norm-method".into(),
        "none".into(),
        "-n".into(),
        prefix.clone(),
    ];
    let expected = invoke(&base);
    assert!(!expected.status.success());
    let diagnostic = fs::read(format!("{prefix}.mle_diagnostics.json")).unwrap();
    for mode in ["1", "2", "4", "8", "auto", "auto:2", "auto:8"] {
        for _ in 0..2 {
            fs::remove_file(format!("{prefix}.mle_diagnostics.json")).unwrap();
            let mut args = base.clone();
            args.extend(["--threads".into(), mode.into()]);
            let actual = invoke(&args);
            assert_eq!(actual.status.code(), expected.status.code());
            assert_eq!(actual.stderr, expected.stderr);
            assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 2);
            assert_eq!(
                fs::read(format!("{prefix}.mle_diagnostics.json")).unwrap(),
                diagnostic
            );
        }
    }
    fs::remove_file(format!("{prefix}.mle_diagnostics.json")).unwrap();
    fs::write(format!("{prefix}.gene_summary.txt"), "sentinel").unwrap();
    let mut args = base;
    args[2] = fixture("tests/fixtures/m11a/counts.tsv");
    args.extend(["--threads".into(), "auto".into()]);
    assert!(!invoke(&args).status.success());
    assert_eq!(
        fs::read_to_string(format!("{prefix}.gene_summary.txt")).unwrap(),
        "sentinel"
    );
    assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 2);
}
#[test]
fn mle_auto_creates_nested_parent_after_analysis() {
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory
        .path()
        .join("new/nested/result")
        .display()
        .to_string();
    let args = vec![
        "mle".into(),
        "-k".into(),
        fixture("tests/fixtures/m11a/counts.tsv"),
        "-d".into(),
        fixture("tests/fixtures/m11a/design.tsv"),
        "--norm-method".into(),
        "none".into(),
        "-n".into(),
        prefix.clone(),
        "--threads".into(),
        "auto".into(),
    ];
    let result = invoke(&args);
    assert!(result.status.success(), "{result:?}");
    assert!(Path::new(&format!("{prefix}.gene_summary.txt")).exists());
    assert!(Path::new(&format!("{prefix}.run_info.txt")).exists());
}
#[test]
fn parser_and_zero_limits_preserve_early_errors() {
    use mageckrs::cli::ThreadMode;
    assert_eq!(
        "auto".parse::<ThreadMode>().unwrap(),
        ThreadMode::Auto { maximum: None }
    );
    assert_eq!(
        "auto:3".parse::<ThreadMode>().unwrap(),
        ThreadMode::Auto { maximum: Some(3) }
    );
    assert_eq!(
        usize::MAX.to_string().parse::<ThreadMode>().unwrap(),
        ThreadMode::Fixed(usize::MAX)
    );
    for bad in ["auto:0", "auto:", "auto:no", "-1", "no"] {
        assert!(bad.parse::<ThreadMode>().is_err());
    }
    for command in ["mle", "test", "count", "run"] {
        let mut args: Vec<String> = match command {
            "mle" => vec![command, "-k", "missing", "-d", "missing"],
            "test" => vec![command, "-k", "missing", "-c", "C", "-t", "T"],
            "count" => vec![command, "-l", "missing", "-f", "missing"],
            _ => vec![
                command, "-l", "missing", "-f", "missing", "-c", "C", "-t", "T",
            ],
        }
        .into_iter()
        .map(str::to_owned)
        .collect();
        args.extend(["-n".into(), "unused".into(), "--threads".into(), "0".into()]);
        let actual = invoke(&args);
        assert_eq!(actual.status.code(), Some(1));
        assert_eq!(
            actual.stderr,
            b"Configuration error: --threads must be greater than zero\n"
        );
    }
}
