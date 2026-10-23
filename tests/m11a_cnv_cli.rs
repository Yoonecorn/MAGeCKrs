//! M11A process-boundary validation, preservation and reproducibility checks.
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
};
fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/m11a")
        .join(name)
}
fn invoke(root: &Path, name: &str, extra: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["mle", "-k"])
        .arg(fixture("counts.tsv"))
        .arg("-d")
        .arg(fixture("design.tsv"))
        .arg("-n")
        .arg(root.join(name))
        .args(["--max-sgrnas-per-gene", "2"])
        .args(extra)
        .output()
        .unwrap()
}
fn corrected_fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/mle_numerical_contract/m11a")
        .join(name)
}
fn contents(root: &Path, name: &str, suffix: &str) -> Vec<u8> {
    fs::read(root.join(format!("{name}.{suffix}"))).unwrap()
}
fn success(output: &Output) {
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty());
}
fn failed(output: &Output, message: &str) {
    assert_eq!(output.status.code(), Some(1), "{output:?}");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains(message), "expected {message:?}: {stderr}");
    assert!(!stderr.contains("panicked"));
    assert!(output.stdout.is_empty());
}
#[test]
fn default_and_all_norm_modes_match_exact_baseline_bytes() {
    let root = tempfile::tempdir().unwrap();
    let controls = fixture("controls.txt");
    for (name, extra) in [
        ("median", vec![]),
        ("none", vec!["--norm-method", "none"]),
        ("total", vec!["--norm-method", "total"]),
        (
            "control",
            vec![
                "--norm-method",
                "control",
                "--control-sgrna",
                controls.to_str().unwrap(),
            ],
        ),
        ("meanvar", vec!["--genes-varmodeling", "2"]),
        ("all-skipped", vec!["--max-sgrnas-per-gene", "1"]),
    ] {
        // Clap rejects repeated nonrepeatable options; all-skipped supplies its own command.
        let output = if name == "all-skipped" {
            Command::new(env!("CARGO_BIN_EXE_mageckrs"))
                .args(["mle", "-k"])
                .arg(fixture("counts.tsv"))
                .arg("-d")
                .arg(fixture("design.tsv"))
                .arg("-n")
                .arg(root.path().join(name))
                .args(&extra)
                .output()
                .unwrap()
        } else {
            invoke(root.path(), name, &extra)
        };
        if name == "all-skipped" {
            failed(&output, "EMPTY_ELIGIBLE_UNIVERSE");
            assert!(!root
                .path()
                .join(format!("{name}.gene_summary.txt"))
                .exists());
            assert!(root
                .path()
                .join(format!("{name}.mle_diagnostics.json"))
                .exists());
            continue;
        }
        success(&output);
        for suffix in ["gene_summary.txt", "sgrna_summary.txt"] {
            assert_eq!(
                contents(root.path(), name, suffix),
                fs::read(corrected_fixture(&format!("corrected-{name}.{suffix}"))).unwrap(),
                "{name}.{suffix}"
            );
        }
        assert_eq!(
            output.stderr,
            fs::read(fixture(&format!("baseline-{name}.stderr"))).unwrap()
        );
        assert!(!root
            .path()
            .join(format!("{name}.cnv_run_info.txt"))
            .exists());
    }
}
#[test]
fn cnv_workflow_preserves_inference_and_repeats_byte_exactly() {
    let root = tempfile::tempdir().unwrap();
    let cnv = fixture("cnv.tsv");
    success(&invoke(root.path(), "plain", &[]));
    let ordinary = String::from_utf8(contents(root.path(), "plain", "gene_summary.txt")).unwrap();
    let mut first = None;
    for i in 0..10 {
        let name = format!("cnv{i}");
        let output = invoke(root.path(), &name, &["--cnv-norm", cnv.to_str().unwrap()]);
        success(&output);
        assert_eq!(
            output.stderr,
            fs::read(fixture("baseline-median.stderr")).unwrap()
        );
        let bytes = ["gene_summary.txt", "sgrna_summary.txt", "cnv_run_info.txt"]
            .map(|s| contents(root.path(), &name, s));
        if let Some(ref previous) = first {
            assert_eq!(&bytes, previous);
        } else {
            first = Some(bytes.clone());
        }
        assert_eq!(
            bytes[1],
            contents(root.path(), "plain", "sgrna_summary.txt")
        );
        let corrected = String::from_utf8(bytes[0].clone()).unwrap();
        assert_eq!(corrected.lines().next(), ordinary.lines().next());
        let mut changes = 0;
        for (left, right) in ordinary.lines().skip(1).zip(corrected.lines().skip(1)) {
            let left: Vec<_> = left.split('\t').collect();
            let right: Vec<_> = right.split('\t').collect();
            assert_eq!(left.len(), right.len());
            for column in 0..left.len() {
                if [2, 8].contains(&column) && left[0] != "G10" {
                    changes += usize::from(left[column] != right[column]);
                } else {
                    assert_eq!(left[column], right[column], "{} column {column}", left[0]);
                }
            }
        }
        assert!(changes >= 10);
        let info = String::from_utf8(bytes[2].clone()).unwrap();
        assert!(info.contains("not recomputed"));
        assert!(info.contains("A\tA\tmatched-genes\t10\n"));
        assert!(info.contains("B\tB\tunmatched-gene\tG10\n"));
        assert!(info.contains("\t\tskipped-gene\tSKIP\n"));
        assert!(!info.contains("matched-gene\tSKIP\n"));
        assert!(info.contains("\tDay0\tsource-sample-index\t3\n"));
        assert!(info.contains("\tExcluded\texcluded-sample\tExcluded\n"));
        assert!(info.contains("raw_count + 1.0"));
        for key in [
            "slope",
            "breakpoint",
            "sse",
            "upper-bound",
            "mu-ref",
            "minimizer",
            "objective-value",
            "iterations",
            "function-evaluations",
            "final-bracket-lower",
            "final-bracket-upper",
            "converged",
            "stop-reason",
        ] {
            assert!(info.contains(&format!("A\tA\t{key}\t")), "{key}");
        }
    }
}
#[test]
fn malformed_matrix_and_population_errors_publish_nothing() {
    for (text, message) in [
        ("", "missing header"),
        ("Gene\tA\nG0\t0\n", "SYMBOL"),
        ("SYMBOL\tA\tA\n", "duplicate profile"),
        ("SYMBOL\t\n", "profile header"),
        ("SYMBOL\tA\nG0\n", "expected 2 fields"),
        ("SYMBOL\tA\n\t0\n", "gene key"),
        ("SYMBOL\tA\nG0\t0\nG0\t1\n", "duplicate gene"),
        ("SYMBOL\tA\nG0\tNaN\n", "non-finite"),
        ("SYMBOL\tA\nG0\tbad\n", "non-numeric"),
        ("SYMBOL\ta\nG0\t0\n", "exact profile"),
        ("SYMBOL\tA\tB\nOTHER\t0\t0\n", "population is empty"),
        ("SYMBOL\tA\tB\nG0\t0\t0\n", "at least two"),
        ("SYMBOL\tA\tB\nG0\t0\t0\nG1\t0\t0\n", "Rank-deficient"),
        ("SYMBOL\tA\tB\nG0\t-1\t0\nG1\t1\t0\n", "low-side"),
        (
            "SYMBOL\tA\tB\nG0\t-2\t0\nG1\t-1\t0\n",
            "Invalid CNV breakpoint bounds",
        ),
        ("SYMBOL\tA\tB\nG0\t1024\t0\n", "overflow"),
    ] {
        let root = tempfile::tempdir().unwrap();
        let cnv = root.path().join("input.tsv");
        fs::write(&cnv, text).unwrap();
        failed(
            &invoke(
                root.path(),
                "nested/out",
                &["--cnv-norm", cnv.to_str().unwrap()],
            ),
            message,
        );
        assert_eq!(fs::read_dir(root.path()).unwrap().count(), 1, "{text:?}");
    }
    let root = tempfile::tempdir().unwrap();
    let cnv = root.path().join("invalid-utf8.tsv");
    fs::write(&cnv, b"SYMBOL\tA\nG0\t\xff\n").unwrap();
    failed(
        &invoke(
            root.path(),
            "nested/out",
            &["--cnv-norm", cnv.to_str().unwrap()],
        ),
        "I/O error:",
    );
    assert!(!root.path().join("nested").exists());
}
#[test]
fn validation_precedence_retains_design_counts_alignment_and_preparation() {
    let root = tempfile::tempdir().unwrap();
    let cnv = root.path().join("bad.cnv");
    fs::write(&cnv, "bad").unwrap();
    let count = root.path().join("count");
    let design = root.path().join("design");
    fs::write(&count, "bad").unwrap();
    fs::write(&design, "bad").unwrap();
    let run = |cnv_path: &Path, extra: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_mageckrs"))
            .args(["mle", "-k"])
            .arg(&count)
            .arg("-d")
            .arg(&design)
            .arg("-n")
            .arg(root.path().join("nested/out"))
            .arg("--cnv-norm")
            .arg(cnv_path)
            .args(extra)
            .output()
            .unwrap()
    };
    failed(
        &run(&root.path().join("missing"), &[]),
        "CNV normalization profile (--cnv-norm) file does not exist",
    );
    failed(&run(&cnv, &[]), "Design");
    fs::copy(fixture("design.tsv"), &design).unwrap();
    failed(&run(&cnv, &[]), "Count table");
    fs::write(
        &count,
        "sgRNA\tGene\tDay0\tTreatA\tTreatB\ng\tG\t1\t2.5\t3\n",
    )
    .unwrap();
    failed(&run(&cnv, &[]), "raw counts must be integers");
    fs::write(
        &count,
        "sgRNA\tGene\tDay0\tTreatA\tmissing\ng\tG\t1\t2\t3\n",
    )
    .unwrap();
    failed(&run(&cnv, &[]), "cannot be found");
    fs::copy(fixture("counts.tsv"), &count).unwrap();
    let controls = root.path().join("controls");
    fs::write(&controls, "unknown\n").unwrap();
    failed(
        &run(
            &cnv,
            &[
                "--norm-method",
                "control",
                "--control-sgrna",
                controls.to_str().unwrap(),
            ],
        ),
        "control",
    );
    failed(&run(&cnv, &[]), "CNV matrix:");
    assert!(!root.path().join("nested").exists());
}
#[test]
fn cnv_schema_precedes_an_independently_failing_mle_fit() {
    let root = tempfile::tempdir().unwrap();
    let counts = root.path().join("counts");
    fs::write(
        &counts,
        "sgRNA\tGene\tDay0\tTreatA\tTreatB\ng\tG0\t1e308\t1e308\t1e308\n",
    )
    .unwrap();
    let bad = root.path().join("cnv");
    fs::write(&bad, "bad").unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_mageckrs"));
    command
        .args(["mle", "-k"])
        .arg(&counts)
        .arg("-d")
        .arg(fixture("design.tsv"))
        .arg("-n")
        .arg(root.path().join("out"))
        .args(["--norm-method", "none"]);
    failed(
        &command.output().unwrap(),
        "IRLS initialization exceeds coefficient bound",
    );
    fs::remove_file(root.path().join("out.mle_diagnostics.json")).unwrap();
    failed(
        &command.arg("--cnv-norm").arg(bad).output().unwrap(),
        "CNV matrix:",
    );
    fs::write(root.path().join("cnv"), "SYMBOL\tA\tB\nG0\t0\t0\n").unwrap();
    failed(&command.output().unwrap(), "at least two");
    assert_eq!(fs::read_dir(root.path()).unwrap().count(), 2);
}
#[test]
fn every_existing_target_is_protected_and_default_ignores_sidecar() {
    for suffix in ["gene_summary.txt", "sgrna_summary.txt", "cnv_run_info.txt"] {
        let root = tempfile::tempdir().unwrap();
        let target = root.path().join(format!("out.{suffix}"));
        fs::write(&target, "ORIGINAL\n").unwrap();
        failed(
            &invoke(
                root.path(),
                "out",
                &["--cnv-norm", fixture("cnv.tsv").to_str().unwrap()],
            ),
            "Output file already exists:",
        );
        assert_eq!(fs::read(&target).unwrap(), b"ORIGINAL\n");
        assert_eq!(fs::read_dir(root.path()).unwrap().count(), 1);
    }
    let root = tempfile::tempdir().unwrap();
    let sidecar = root.path().join("out.cnv_run_info.txt");
    fs::write(&sidecar, "ORIGINAL\n").unwrap();
    success(&invoke(root.path(), "out", &[]));
    assert_eq!(fs::read(sidecar).unwrap(), b"ORIGINAL\n");
}
#[test]
fn cnv_row_order_does_not_change_cnv_outputs() {
    let root = tempfile::tempdir().unwrap();
    let original = fixture("cnv.tsv");
    success(&invoke(
        root.path(),
        "first",
        &["--cnv-norm", original.to_str().unwrap()],
    ));
    let data = fs::read_to_string(original).unwrap();
    let mut lines = data.lines();
    let header = lines.next().unwrap();
    let mut rows: Vec<_> = lines.collect();
    rows.reverse();
    let reordered = root.path().join("reordered");
    fs::write(&reordered, format!("{header}\n{}\n", rows.join("\n"))).unwrap();
    success(&invoke(
        root.path(),
        "second",
        &["--cnv-norm", reordered.to_str().unwrap()],
    ));
    for suffix in ["gene_summary.txt", "sgrna_summary.txt", "cnv_run_info.txt"] {
        assert_eq!(
            contents(root.path(), "first", suffix),
            contents(root.path(), "second", suffix)
        );
    }
}

#[test]
fn cnv_normalization_and_meanvar_opt_in_preserve_sgrna_bytes() {
    let root = tempfile::tempdir().unwrap();
    for (i, extra) in [
        vec!["--norm-method", "none"],
        vec!["--norm-method", "total"],
        vec![
            "--norm-method",
            "control",
            "--control-sgrna",
            fixture("controls.txt").to_str().unwrap(),
        ],
        vec!["--genes-varmodeling", "2"],
    ]
    .iter()
    .enumerate()
    {
        let plain = format!("plain{i}");
        success(&invoke(root.path(), &plain, extra));
        let mut corrected = extra.clone();
        let cnv = fixture("cnv.tsv");
        corrected.extend(["--cnv-norm", cnv.to_str().unwrap()]);
        let cnv_name = format!("corrected{i}");
        let result = invoke(root.path(), &cnv_name, &corrected);
        if extra.contains(&"control") {
            failed(&result, "only 0 low-side observations");
            assert!(!root
                .path()
                .join(format!("{cnv_name}.gene_summary.txt"))
                .exists());
            continue;
        }
        success(&result);
        assert_eq!(
            contents(root.path(), &plain, "sgrna_summary.txt"),
            contents(root.path(), &cnv_name, "sgrna_summary.txt")
        );
    }
}
#[test]
fn cnv_population_is_selected_after_large_gene_filtering() {
    let root = tempfile::tempdir().unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["mle", "-k"])
        .arg(fixture("counts.tsv"))
        .arg("-d")
        .arg(fixture("design.tsv"))
        .arg("-n")
        .arg(root.path().join("out"))
        .args(["--max-sgrnas-per-gene", "1", "--cnv-norm"])
        .arg(fixture("cnv.tsv"))
        .output()
        .unwrap();
    failed(&output, "population is empty");
    assert_eq!(fs::read_dir(root.path()).unwrap().count(), 0);
}

#[test]
fn design_sample_order_controls_fitting_when_count_columns_are_reordered() {
    let root = tempfile::tempdir().unwrap();
    let output = invoke(
        root.path(),
        "first",
        &["--cnv-norm", fixture("cnv.tsv").to_str().unwrap()],
    );
    success(&output);
    let text = fs::read_to_string(fixture("counts.tsv")).unwrap();
    let order = [0, 1, 5, 3, 2, 4];
    let mut reordered = String::new();
    for line in text.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        let fields: Vec<_> = order.iter().map(|&i| fields[i]).collect();
        reordered.push_str(&fields.join("\t"));
        reordered.push('\n');
    }
    let counts = root.path().join("reordered-counts");
    fs::write(&counts, reordered).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["mle", "-k"])
        .arg(counts)
        .arg("-d")
        .arg(fixture("design.tsv"))
        .arg("-n")
        .arg(root.path().join("second"))
        .args(["--max-sgrnas-per-gene", "2", "--cnv-norm"])
        .arg(fixture("cnv.tsv"))
        .output()
        .unwrap();
    success(&output);
    for suffix in ["gene_summary.txt", "sgrna_summary.txt"] {
        assert_eq!(
            contents(root.path(), "first", suffix),
            contents(root.path(), "second", suffix)
        );
    }
    let filter = |bytes: Vec<u8>| {
        String::from_utf8(bytes)
            .unwrap()
            .lines()
            .filter(|l| !l.contains("source-sample-index"))
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    assert_eq!(
        filter(contents(root.path(), "first", "cnv_run_info.txt")),
        filter(contents(root.path(), "second", "cnv_run_info.txt"))
    );
    let info = String::from_utf8(contents(root.path(), "second", "cnv_run_info.txt")).unwrap();
    assert!(info.contains("\tDay0\tsource-sample-index\t0\n"));
}

#[test]
fn three_conditions_require_and_correct_three_distinct_profiles() {
    let root = tempfile::tempdir().unwrap();
    let counts = root.path().join("counts");
    let text = fs::read_to_string(fixture("counts.tsv")).unwrap();
    let mut data = String::new();
    for (i, line) in text.lines().enumerate() {
        data.push_str(line);
        data.push('\t');
        data.push_str(if i == 0 {
            "TreatC"
        } else {
            line.split('\t').nth(2).unwrap()
        });
        data.push('\n');
    }
    fs::write(&counts, data).unwrap();
    let design = root.path().join("design");
    fs::write(&design,"Samples\tbaseline\tA\tB\tC\nDay0\t1\t0\t0\t0\nTreatA\t1\t1\t0\t0\nTreatB\t1\t0\t1\t0\nTreatC\t1\t0\t0\t1\n").unwrap();
    let cnv = root.path().join("cnv");
    let text = fs::read_to_string(fixture("cnv.tsv")).unwrap();
    let mut data = String::new();
    for (i, line) in text.lines().enumerate() {
        data.push_str(line);
        data.push('\t');
        data.push_str(if i == 0 {
            "C"
        } else {
            line.split('\t').nth(1).unwrap()
        });
        data.push('\n');
    }
    fs::write(&cnv, data).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["mle", "-k"])
        .arg(counts)
        .arg("-d")
        .arg(design)
        .arg("-n")
        .arg(root.path().join("out"))
        .args(["--max-sgrnas-per-gene", "2", "--cnv-norm"])
        .arg(cnv)
        .output()
        .unwrap();
    success(&output);
    let info = String::from_utf8(contents(root.path(), "out", "cnv_run_info.txt")).unwrap();
    for condition in ["A", "B", "C"] {
        assert!(info.contains(&format!("{condition}\t{condition}\tmatched-genes\t10\n")));
    }
    let header = String::from_utf8(contents(root.path(), "out", "gene_summary.txt")).unwrap();
    assert!(header.lines().next().unwrap().contains("C|beta"));
}
