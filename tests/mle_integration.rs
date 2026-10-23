//! External RUN1 coverage is intentionally limited to observable CLI behavior;
//! deterministic publication fault injection remains private to `mle::run`.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/mle")
        .join(name)
}
fn temp_prefix(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("mageckrs-run1-{name}-{}", std::process::id()))
}
fn mle(prefix: &PathBuf, method: &str) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_mageckrs"));
    command
        .args(["mle", "-k"])
        .arg(fixture("count_minimal.tsv"))
        .args(["-d"])
        .arg(fixture("design_minimal.tsv"))
        .args(["-n"])
        .arg(prefix)
        .args(["--norm-method", method]);
    command
}

#[test]
fn mle_help_is_available() {
    let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["mle", "--help"])
        .output()
        .expect("run mle help");
    assert!(output.status.success());
    assert!(String::from_utf8_lossy(&output.stdout).contains("--count-table"));
}

#[test]
fn run1_oracles_01_09_10_11_publish_all_normalization_modes() {
    for method in ["none", "median", "total", "control"] {
        let prefix = temp_prefix(method);
        let mut command = mle(&prefix, method);
        if method == "control" {
            command
                .args(["--control-sgrna"])
                .arg(fixture("controls.txt"));
        }
        let output = command.output().unwrap();
        assert!(
            output.status.success(),
            "{method}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let gene = fs::read_to_string(format!("{}.gene_summary.txt", prefix.display())).unwrap();
        let sgrna = fs::read_to_string(format!("{}.sgrna_summary.txt", prefix.display())).unwrap();
        assert!(gene.starts_with("Gene\tsgRNA\ttreatment|beta"));
        assert_eq!(gene.contains("GENE_CTRL"), method != "control");
        assert!(sgrna.starts_with("Gene\tsgRNA\teff\n") && sgrna.contains("GENE_CTRL\tgc1"));
        let _ = fs::remove_file(format!("{}.gene_summary.txt", prefix.display()));
        let _ = fs::remove_file(format!("{}.sgrna_summary.txt", prefix.display()));
        let _ = fs::remove_file(format!("{}.mle_diagnostics.json", prefix.display()));
    }
}

#[test]
fn run1_oracles_02_03_and_08_are_observable() {
    for (suffix, existing) in [("gene", "gene_summary.txt"), ("sgrna", "sgrna_summary.txt")] {
        let prefix = temp_prefix(suffix);
        let target = format!("{}.{}", prefix.display(), existing);
        fs::write(&target, "ORIGINAL_USER_DATA\n").unwrap();
        let output = mle(&prefix, "none").output().unwrap();
        assert!(!output.status.success());
        assert_eq!(fs::read_to_string(&target).unwrap(), "ORIGINAL_USER_DATA\n");
        let _ = fs::remove_file(target);
    }
    for args in [
        ["--cnv-norm", "/dev/null"],
        ["--threads", "0"],
        ["--seed", "42"],
    ] {
        let output = mle(&temp_prefix("unsupported"), "none")
            .args(args)
            .output()
            .unwrap();
        assert!(!output.status.success());
    }
}

#[test]
fn run1_oracles_12_13_14_reject_before_publication() {
    let root = temp_prefix("invalid");
    fs::create_dir_all(&root).unwrap();
    let count = root.join("count.tsv");
    fs::write(&count, "sgRNA\tGene\tc0\tt0\ng1\tG\t10\t50.5\n").unwrap();
    let rank = root.join("rank.tsv");
    fs::write(&rank, "sample\tbaseline\ta\tb\nc0\t1\t0\t0\nt0\t1\t1\t2\n").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["mle", "-k"])
        .arg(&count)
        .args(["-d"])
        .arg(&rank)
        .args(["-n"])
        .arg(root.join("rank"))
        .output()
        .unwrap();
    assert!(!output.status.success());
    let valid = fixture("design_minimal.tsv");
    let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["mle", "-k"])
        .arg(&count)
        .args(["-d"])
        .arg(valid)
        .args(["-n"])
        .arg(root.join("fractional"))
        .output()
        .unwrap();
    assert!(!output.status.success());
    let missing = root.join("missing.tsv");
    fs::write(
        &missing,
        "sample\tbaseline\ttreatment\nc0\t1\t0\ntreat_MISSING\t1\t1\n",
    )
    .unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["mle", "-k"])
        .arg(&count)
        .args(["-d"])
        .arg(&missing)
        .args(["-n"])
        .arg(root.join("missing"))
        .output()
        .unwrap();
    assert!(!output.status.success());
    let _ = fs::remove_dir_all(root);
}
