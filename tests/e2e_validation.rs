//! System-level regression coverage for the implemented `count` command.
//!
//! Fixture provenance: the three inputs in `fixtures/e2e_validation/` are
//! hand-authored, repository-local minimal examples. The library has two
//! 20-nt guides; baseline contains A, A, C, and an unmapped T read, while
//! treatment contains C, C, and an unmapped T read. Therefore the expected
//! count matrix is `(g1, g2) = (2, 1)` for baseline and `(0, 2)` for treatment;
//! with `--norm-method none`, normalized values are those counts formatted to
//! four decimal places. Per-file mapped fractions are 3/4 = 0.7500 and
//! 2/3 = 0.6667, respectively.

use std::fs;
use std::path::Path;
use std::process::Command;

use tempfile::TempDir;

const LIBRARY: &str = include_str!("fixtures/e2e_validation/library.tsv");
const BASELINE_FASTQ: &str = include_str!("fixtures/e2e_validation/baseline.fastq");
const TREATMENT_FASTQ: &str = include_str!("fixtures/e2e_validation/treatment.fastq");

fn materialize_fixture(dir: &TempDir, name: &str, contents: &str) -> String {
    let path = dir.path().join(name);
    fs::write(&path, contents).expect("write deterministic fixture");
    path.to_string_lossy().into_owned()
}

fn run_count(prefix: &str, library: &str, baseline: &str, treatment: &str) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args([
            "count",
            "--list-seq",
            library,
            "--fastq",
            baseline,
            treatment,
            "--sample-label",
            "baseline,treatment",
            "--output-prefix",
            prefix,
            "--norm-method",
            "none",
        ])
        .output()
        .expect("execute mageckrs count")
}

#[test]
fn count_cli_processes_fixture_and_is_repeatable() {
    let dir = tempfile::tempdir().expect("create tempdir");
    let library = materialize_fixture(&dir, "library.tsv", LIBRARY);
    let baseline = materialize_fixture(&dir, "baseline.fastq", BASELINE_FASTQ);
    let treatment = materialize_fixture(&dir, "treatment.fastq", TREATMENT_FASTQ);
    let prefix = dir.path().join("screen").to_string_lossy().into_owned();

    let first = run_count(&prefix, &library, &baseline, &treatment);
    assert!(
        first.status.success(),
        "count command failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );

    let count_path = format!("{prefix}.count.txt");
    let normalized_path = format!("{prefix}.count.normalized.txt");
    let summary_path = format!("{prefix}.countsummary.txt");
    let first_count = fs::read_to_string(&count_path).expect("read raw count output");
    let first_normalized = fs::read_to_string(&normalized_path).expect("read normalized output");
    let first_summary = fs::read_to_string(&summary_path).expect("read summary output");

    assert_eq!(
        first_count,
        "sgRNA\tGene\tbaseline\ttreatment\n\
         g1\tGENE_A\t2\t0\n\
         g2\tGENE_B\t1\t2\n"
    );
    assert_eq!(
        first_normalized,
        "sgRNA\tGene\tbaseline\ttreatment\n\
         g1\tGENE_A\t2.0000\t0.0000\n\
         g2\tGENE_B\t1.0000\t2.0000\n"
    );

    let summary_lines: Vec<&str> = first_summary.lines().collect();
    assert_eq!(
        summary_lines[0],
        "File\tLabel\tReads\tMapped\tPercentage\tZerocounts"
    );
    assert_eq!(summary_lines.len(), 3);
    assert!(summary_lines[1].ends_with("\tbaseline\t4\t3\t0.7500\t0"));
    assert!(summary_lines[2].ends_with("\ttreatment\t3\t2\t0.6667\t1"));

    let second = run_count(&prefix, &library, &baseline, &treatment);
    assert!(
        second.status.success(),
        "repeat count command failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(fs::read_to_string(&count_path).unwrap(), first_count);
    assert_eq!(
        fs::read_to_string(&normalized_path).unwrap(),
        first_normalized
    );
    assert_eq!(fs::read_to_string(&summary_path).unwrap(), first_summary);
}

#[test]
fn count_cli_rejects_duplicate_sample_labels_before_writing_output() {
    let dir = tempfile::tempdir().expect("create tempdir");
    let library = materialize_fixture(&dir, "library.tsv", LIBRARY);
    let baseline = materialize_fixture(&dir, "baseline.fastq", BASELINE_FASTQ);
    let treatment = materialize_fixture(&dir, "treatment.fastq", TREATMENT_FASTQ);
    let prefix = dir.path().join("invalid").to_string_lossy().into_owned();

    let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args([
            "count",
            "--list-seq",
            &library,
            "--fastq",
            &baseline,
            &treatment,
            "--sample-label",
            "duplicate,duplicate",
            "--output-prefix",
            &prefix,
            "--norm-method",
            "none",
        ])
        .output()
        .expect("execute invalid mageckrs count invocation");

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("Duplicate sample label 'duplicate'"));
    for suffix in ["count.txt", "count.normalized.txt", "countsummary.txt"] {
        assert!(
            !Path::new(&format!("{prefix}.{suffix}")).exists(),
            "invalid configuration created an output file"
        );
    }
}
