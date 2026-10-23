//! Process and argument coverage for the RRA test subcommand.

use clap::Parser;
use mageckrs::cli::{Cli, Commands};
use mageckrs::config::TestConfig;
use std::fs;

fn parse_test(extra: &[&str]) -> Result<mageckrs::cli::TestArgs, clap::Error> {
    let mut args = vec![
        "mageckrs",
        "test",
        "--count-table",
        "counts.tsv",
        "--treatment-id",
        "treatment",
        "--day0-label",
        "baseline",
        "--output-prefix",
        "result",
    ];
    args.extend_from_slice(extra);
    match Cli::try_parse_from(args)?.command {
        Commands::Test(args) => Ok(args),
        _ => unreachable!(),
    }
}

fn config(extra: &[&str]) -> TestConfig {
    let dir = tempfile::tempdir().unwrap();
    let table = dir.path().join("counts.tsv");
    fs::write(&table, "sgRNA\tGene\tbaseline\ttreatment\ng\tA\t1\t2\n").unwrap();
    let mut args = parse_test(extra).unwrap();
    args.count_table = table;
    TestConfig::try_from(args).unwrap()
}

#[test]
fn parsing_rra_defaults_and_valid_u64() {
    let default = config(&[]);
    assert_eq!(default.rra_seed, 0x0123_4567_89ab_cdef);
    assert_eq!(default.permutation, 1000);
    for seed in ["0", "42", "18446744073709551615"] {
        let config = config(&["--rra-seed", seed, "--num-permutations", "17"]);
        assert_eq!(config.rra_seed, seed.parse::<u64>().unwrap());
        assert_eq!(config.permutation, 17);
    }
}

#[test]
fn parsing_rejects_invalid_seed_and_permutation_values() {
    for seed in ["abc", "-1", "18446744073709551616", "1.5"] {
        assert!(parse_test(&["--rra-seed", seed]).is_err());
    }
    assert!(parse_test(&["--num-permutations", "abc"]).is_err());
    let dir = tempfile::tempdir().unwrap();
    let table = dir.path().join("counts.tsv");
    fs::write(&table, "header").unwrap();
    let mut args = parse_test(&["--num-permutations", "0"]).unwrap();
    args.count_table = table;
    assert!(TestConfig::try_from(args)
        .unwrap_err()
        .to_string()
        .contains("greater than 0"));
}

#[test]
fn parsing_preserves_legacy_seed_and_permutation_and_rejects_conflicts() {
    let legacy = config(&["--seed", "99", "--permutation", "7"]);
    assert_eq!(legacy.rra_seed, 99);
    assert_eq!(legacy.permutation, 7);
    assert_eq!(legacy.seed, Some(99));
    assert!(parse_test(&["--seed", "1", "--rra-seed", "2"]).is_err());
    assert!(parse_test(&["--permutation", "1", "--num-permutations", "2"]).is_err());
}

#[test]
fn parsing_test_help_documents_new_flags_and_scope() {
    let error = Cli::try_parse_from(["mageckrs", "test", "--help"])
        .unwrap_err()
        .to_string();
    assert!(error.contains("--rra-seed"));
    assert!(error.contains("default: 0x0123456789abcdef"));
    assert!(error.contains("--num-permutations"));
    assert!(error.contains("default: 1000"));
    assert!(!error.contains("--fastq"));
    assert!(!error.contains("--design-matrix"));
}

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const GENE_HEADER: &str = "id\tnum\tneg|score\tneg|p-value\tneg|fdr\tneg|rank\tneg|lfc\tneg|goodsgrna\tpos|score\tpos|p-value\tpos|fdr\tpos|rank\tpos|lfc\tpos|goodsgrna";
const GUIDE_HEADER: &str = "sgRNA\tGene\tcontrol_count\ttreatment_count\tcontrol_mean\ttreat_mean\tLFC\tcontrol_var\tadj_var\tscore\tp.low\tp.high\tp.twosided\tFDR\thigh_in_treatment";

fn fixture(dir: &Path) -> PathBuf {
    for (name, contents) in [
        (
            "library.tsv",
            include_str!("fixtures/e2e_validation/library.tsv"),
        ),
        (
            "baseline.fastq",
            include_str!("fixtures/e2e_validation/baseline.fastq"),
        ),
        (
            "treatment.fastq",
            include_str!("fixtures/e2e_validation/treatment.fastq"),
        ),
    ] {
        fs::write(dir.join(name), contents).unwrap();
    }
    let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["count", "--list-seq"])
        .arg(dir.join("library.tsv"))
        .arg("--fastq")
        .arg(dir.join("baseline.fastq"))
        .arg(dir.join("treatment.fastq"))
        .args(["--sample-label", "baseline,treatment", "--output-prefix"])
        .arg(dir.join("counts"))
        .args(["--norm-method", "none"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    dir.join("counts.count.txt")
}

fn run_test(table: &Path, prefix: &Path, extra: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["test", "--count-table"])
        .arg(table)
        .args([
            "--treatment-id",
            "treatment",
            "--day0-label",
            "baseline",
            "--output-prefix",
        ])
        .arg(prefix)
        .args(extra)
        .output()
        .unwrap()
}

fn read(prefix: &Path, suffix: &str) -> String {
    fs::read_to_string(format!("{}.{suffix}", prefix.display())).unwrap()
}

fn assert_no_outputs(prefix: &Path) {
    for suffix in ["gene_summary.txt", "sgrna_summary.txt", "run_info.txt"] {
        assert!(!Path::new(&format!("{}.{suffix}", prefix.display())).exists());
    }
}

#[test]
fn count_to_test_cli_writes_real_summaries_and_ordered_metadata() {
    let dir = tempfile::tempdir().unwrap();
    let table = fixture(dir.path());
    let prefix = dir.path().join("rra");
    // Exact owner-requested test invocation, using the fixture's single control.
    let output = run_test(&table, &prefix, &[]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let gene = read(&prefix, "gene_summary.txt");
    let guide = read(&prefix, "sgrna_summary.txt");
    assert_eq!(gene.lines().next().unwrap(), GENE_HEADER);
    assert_eq!(gene.lines().count(), 3);
    assert_eq!(guide.lines().next().unwrap(), GUIDE_HEADER);
    assert_eq!(guide.lines().count(), 3);
    let guide_rows: Vec<Vec<&str>> = guide
        .lines()
        .skip(1)
        .map(|r| r.split('\t').collect())
        .collect();
    assert_eq!(guide_rows[0][0], "g1");
    assert_eq!(guide_rows[1][0], "g2");
    assert_eq!(guide_rows[0][6], "-1.5850");
    assert_eq!(guide_rows[1][6], "0.5850");
    for row in guide_rows {
        assert_eq!(row[7], "0.00000000e0");
        assert_eq!(row[8], "0.00000000e0");
        assert_eq!(&row[10..14], &["NA"; 4]);
    }
    assert_eq!(read(&prefix, "run_info.txt"), format!("key\tvalue\nmageckrs_version\t{}\ncommand\ttest\nrra_seed\t0x0123456789abcdef\nnum_genes\t2\nnum_guides\t2\ncontrol_null\tfalse\noutput_prefix\t{}\n", env!("CARGO_PKG_VERSION"), prefix.display()));
}

#[test]
fn explicit_seed_reproduces_outputs_and_legacy_flags_remain_compatible() {
    let dir = tempfile::tempdir().unwrap();
    let table = fixture(dir.path());
    let first = dir.path().join("first");
    let second = dir.path().join("second");
    assert!(run_test(
        &table,
        &first,
        &["--rra-seed", "42", "--num-permutations", "37"]
    )
    .status
    .success());
    assert!(
        run_test(&table, &second, &["--seed", "42", "--permutation", "37"])
            .status
            .success()
    );
    for suffix in ["gene_summary.txt", "sgrna_summary.txt"] {
        assert_eq!(read(&first, suffix), read(&second, suffix));
    }
    assert!(read(&first, "run_info.txt").contains("rra_seed\t0x000000000000002a\n"));
}

#[test]
fn control_file_selects_control_null_and_reports_metadata() {
    let dir = tempfile::tempdir().unwrap();
    let table = fixture(dir.path());
    let controls = dir.path().join("controls.txt");
    fs::write(&controls, "# controls\ng1\ng1\nunknown\ng2\n").unwrap();
    let prefix = dir.path().join("control");
    let output = run_test(
        &table,
        &prefix,
        &["--control-sgrna", controls.to_str().unwrap()],
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(read(&prefix, "run_info.txt").contains("control_null\ttrue\n"));
    let ordinary = dir.path().join("ordinary");
    assert!(run_test(&table, &ordinary, &[]).status.success());
    assert_ne!(
        read(&prefix, "gene_summary.txt"),
        read(&ordinary, "gene_summary.txt")
    );
}

#[test]
fn invalid_seed_mask_and_missing_control_pool_fail_cleanly() {
    let dir = tempfile::tempdir().unwrap();
    let table = fixture(dir.path());
    let controls = dir.path().join("unknown.txt");
    fs::write(&controls, "unknown\n").unwrap();
    for (i, extra) in [
        vec!["--rra-seed", "invalid"],
        vec!["--remove-zero-threshold", "100"],
        vec!["--control-sgrna", controls.to_str().unwrap()],
        vec!["--adjust-method", "unsupported"],
    ]
    .iter()
    .enumerate()
    {
        let prefix = dir.path().join(format!("invalid-{i}"));
        let output = run_test(&table, &prefix, extra);
        assert!(!output.status.success());
        let error = String::from_utf8_lossy(&output.stderr);
        assert!(!error.contains("panicked"));
        assert!(!error.is_empty());
        assert_no_outputs(&prefix);
    }
}

#[test]
fn preexisting_metadata_rejected_without_publishing_summaries() {
    let dir = tempfile::tempdir().unwrap();
    let table = fixture(dir.path());
    let prefix = dir.path().join("existing");
    let metadata = dir.path().join("existing.run_info.txt");
    fs::write(&metadata, "keep").unwrap();
    let output = run_test(&table, &prefix, &[]);
    assert!(!output.status.success());
    assert_eq!(fs::read_to_string(metadata).unwrap(), "keep");
    assert!(!dir.path().join("existing.gene_summary.txt").exists());
    assert!(!dir.path().join("existing.sgrna_summary.txt").exists());
}

#[test]
fn two_control_columns_and_invalid_table_errors_use_existing_apis() {
    let dir = tempfile::tempdir().unwrap();
    let table = dir.path().join("counts.tsv");
    fs::write(
        &table,
        "sgRNA\tGene\tbaseline\tbaseline2\ttreatment\ng1\tA\t2\t4\t1\ng2\tB\t1\t3\t5\n",
    )
    .unwrap();
    let prefix = dir.path().join("multi");
    let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["test", "--count-table"])
        .arg(&table)
        .args([
            "--treatment-id",
            "treatment",
            "--day0-label",
            "baseline,baseline2",
            "--output-prefix",
        ])
        .arg(&prefix)
        .args(["--num-permutations", "10"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let guide = read(&prefix, "sgrna_summary.txt");
    assert_eq!(
        guide.lines().nth(1).unwrap().split('\t').nth(7).unwrap(),
        "2.00000000e0"
    );
    fs::write(&table, "bad header\n").unwrap();
    let invalid = dir.path().join("bad");
    let output = run_test(&table, &invalid, &[]);
    assert!(!output.status.success());
    assert!(!String::from_utf8_lossy(&output.stderr).contains("panicked"));
    assert_no_outputs(&invalid);
}
