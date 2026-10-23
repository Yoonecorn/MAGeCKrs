use std::{fs, path::Path, process::Command};

fn count(binary: &str, library: &Path, fastqs: &[&Path], prefix: &Path, mode: Option<&str>) {
    let mut command = Command::new(binary);
    command
        .args(["count", "--list-seq"])
        .arg(library)
        .arg("--fastq");
    command
        .args(fastqs)
        .arg("--output-prefix")
        .arg(prefix)
        .args(["--norm-method", "none"]);
    if let Some(mode) = mode {
        command.args(["--trim-5", mode]);
    }
    let output = command.env("RUST_LOG", "info").output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("guide start offsets"));
}

#[test]
fn modes_preserve_default_and_auto_matches_variable_offsets_across_samples() {
    let dir = tempfile::tempdir().unwrap();
    let library = dir.path().join("library.tsv");
    let first = dir.path().join("first.fq");
    let second = dir.path().join("second.fq");
    fs::write(
        &library,
        "id\tsequence\tgene\ng1\tACGTACGTACGTACGTACGT\tG1\n",
    )
    .unwrap();
    let record = "@r\nNNNNNACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIII\n";
    fs::write(&first, record).unwrap();
    let other = "@r\nNNNNNNNACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    fs::write(&second, format!("{record}{other}")).unwrap();
    for mode in [None, Some("none"), Some("0"), Some("5"), Some("auto")] {
        let name = mode.unwrap_or("default");
        let prefix = dir.path().join(name);
        count(
            env!("CARGO_BIN_EXE_mageckrs"),
            &library,
            &[&first, &second],
            &prefix,
            mode,
        );
        let table = fs::read_to_string(prefix.with_extension("count.txt")).unwrap();
        let metadata = fs::read_to_string(prefix.with_extension("run_info.txt")).unwrap();
        match mode {
            Some("auto") => {
                assert!(table.contains("g1\tG1\t1\t2\n"));
                assert!(metadata.contains("trim_5_offsets\t0..=40\n"));
            }
            Some("5") => {
                assert!(table.contains("g1\tG1\t1\t1\n"));
                assert!(metadata.contains("trim_5_offsets\t5\n"));
            }
            _ => {
                assert!(table.contains("g1\tG1\t0\t0\n"));
                assert!(metadata.contains("trim_5_offsets\t0\n"));
            }
        }
    }
    for suffix in ["count.txt", "count.normalized.txt", "countsummary.txt"] {
        assert_eq!(
            fs::read(dir.path().join(format!("default.{suffix}"))).unwrap(),
            fs::read(dir.path().join(format!("none.{suffix}"))).unwrap()
        );
    }
    for invalid in ["-1", "invalid", "1.5"] {
        assert!(!Command::new(env!("CARGO_BIN_EXE_mageckrs"))
            .args(["count", "--list-seq"])
            .arg(&library)
            .arg("--fastq")
            .arg(&first)
            .args(["--output-prefix", "unused", "--trim-5", invalid])
            .output()
            .unwrap()
            .status
            .success());
    }
}
