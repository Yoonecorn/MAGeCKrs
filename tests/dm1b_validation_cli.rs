//! DM1B public-boundary regression tests: exact parser/alignment diagnostics,
//! component-only floating rank diagnostics, and unchanged publication behavior.
use clap::Parser;
use mageckrs::{cli::Cli, error::MageckError};
use std::{fs, process::Command};

fn reject(design: &str, counts: &str, config: bool, expected: &str, exact: bool) {
    let root = tempfile::tempdir().unwrap();
    let d = root.path().join("design.tsv");
    let k = root.path().join("counts.tsv");
    let prefix = root.path().join("absent").join("result");
    fs::write(&d, design).unwrap();
    fs::write(&k, counts).unwrap();
    let args = [
        "mageckrs",
        "mle",
        "-d",
        d.to_str().unwrap(),
        "-k",
        k.to_str().unwrap(),
        "-n",
        prefix.to_str().unwrap(),
    ];
    let error = mageckrs::run_with_cli(Cli::try_parse_from(args).unwrap()).unwrap_err();
    let message = match error {
        MageckError::Config(message) if config => message,
        MageckError::DesignMatrix(message) if !config => message,
        error => panic!("wrong category: {error:?}"),
    };
    if exact {
        assert_eq!(message, expected);
    } else {
        assert!(message.contains(expected), "{message}");
    }
    let mut stderr = None;
    for _ in 0..2 {
        let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
            .args(&args[1..])
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1));
        assert!(output.stdout.is_empty());
        let category = if config {
            "Configuration error: "
        } else {
            "Design matrix error: "
        };
        let text = String::from_utf8(output.stderr.clone()).unwrap();
        assert!(text.contains(&format!("{category}{expected}")), "{text}");
        if let Some(previous) = &stderr {
            assert_eq!(previous, &output.stderr);
        }
        stderr = Some(output.stderr);
        assert!(!prefix.parent().unwrap().exists());
        assert_eq!(fs::read_dir(root.path()).unwrap().count(), 2);
    }
}

#[test]
fn design_validation_precedes_counts_and_alignment() {
    let missing = "sgRNA\tGene\tS1\ns1\tG1\t10\n";
    // ADV13: duplicate Day0 is the earlier error even when count Day0 is absent.
    reject(
        "s\ti\ta\nDay0\t1\t0\nDay0\t1\tNaN\n",
        missing,
        true,
        "Line 3, column 1: duplicate sample row label 'Day0' in design matrix",
        true,
    );
    reject(
        "s\ti\ta\nDay0\t1\t0\nS1\t2\tbad\n",
        "malformed\n",
        false,
        "Line 3, column 2: intercept must evaluate to 1.0, got 2",
        true,
    );
    reject(
        "s\ti\ta\nDay0\t1\t0\n",
        missing,
        false,
        "Design matrix requires at least 2 sample rows, got 1",
        true,
    );
    reject(
        "s\ti\ta\tb\nDay0\t1\t0\t0\nS1\t1\t1\t0\n",
        missing,
        false,
        "Design matrix is not full column rank: m=2, n=3",
        true,
    );
    reject("s\ti\ta\nDay0\t1\t0\nS1\t1\t0\n", "malformed\n", false,
        "Design matrix is not full column rank: m=2, n=2; dependent column index 1, label 'a'; residual norm ", false);
    reject(
        "s\ti\ta\nDay0\t1\t0\nS1\t1\t1\n",
        missing,
        true,
        "Sample 'Day0' specified in design matrix cannot be found in count table",
        true,
    );
    // Label matching is exact: spaces remain significant at alignment.
    reject(
        "s\ti\ta\n Day0 \t1\t0\nS1\t1\t1\n",
        "sgRNA\tGene\tDay0\tS1\ns1\tG1\t10\t20\n",
        true,
        "Sample ' Day0 ' specified in design matrix cannot be found in count table",
        true,
    );
}

#[test]
fn excluded_count_column_still_requires_valid_syntax() {
    reject(
        "s\ti\ta\nDay0\t1\t0\nS1\t1\t1\n",
        "sgRNA\tGene\tDay0\tS1\tExcluded\ns1\tG1\t10\t20\tbad\n",
        true,
        "Line 2: sample 'Excluded' has non-numeric count 'bad'",
        true,
    );
}

#[test]
fn signed_zero_design_preserves_ordinary_and_cnv_output_bytes() {
    let root = tempfile::tempdir().unwrap();
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/m11a");
    let zero = root.path().join("signed-zero.tsv");
    fs::write(
        &zero,
        fs::read_to_string(fixture.join("design.tsv"))
            .unwrap()
            .replace("\t0", "\t-0.0"),
    )
    .unwrap();
    for cnv in [false, true] {
        let mut outputs = Vec::new();
        for (label, design) in [
            ("positive", fixture.join("design.tsv")),
            ("negative", zero.clone()),
        ] {
            let directory = root.path().join(format!("{cnv}-{label}"));
            fs::create_dir(&directory).unwrap();
            let mut args = vec![
                "mageckrs".to_owned(),
                "mle".into(),
                "-k".into(),
                fixture.join("counts.tsv").display().to_string(),
                "-d".into(),
                design.display().to_string(),
                "-n".into(),
                directory.join("result").display().to_string(),
            ];
            if cnv {
                args.extend([
                    "--cnv-norm".into(),
                    fixture.join("cnv.tsv").display().to_string(),
                ]);
            }
            mageckrs::run_with_cli(Cli::try_parse_from(args).unwrap()).unwrap();
            let mut files: Vec<_> = fs::read_dir(&directory)
                .unwrap()
                .filter(|entry| {
                    !entry
                        .as_ref()
                        .unwrap()
                        .file_name()
                        .to_string_lossy()
                        .ends_with("mle_diagnostics.json")
                })
                .map(|entry| {
                    let p = entry.unwrap().path();
                    (p.file_name().unwrap().to_owned(), fs::read(p).unwrap())
                })
                .collect();
            files.sort();
            assert_eq!(files.len(), if cnv { 3 } else { 2 });
            assert!(directory.join("result.mle_diagnostics.json").exists());
            outputs.push(files);
        }
        assert_eq!(outputs[0], outputs[1]);
    }
}
