use std::{fs, path::Path, process::Command};
fn invoke(root: &Path, counts: &str, extra: &[&str], threads: usize) -> std::process::Output {
    fs::write(root.join("count.tsv"), counts).unwrap();
    fs::write(
        root.join("design.tsv"),
        "Samples\tbaseline\ttreatment\ns0\t1\t0\ns1\t1\t1\ns2\t1\t0\ns3\t1\t1\n",
    )
    .unwrap();
    Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(["mle", "-k"])
        .arg(root.join("count.tsv"))
        .arg("-d")
        .arg(root.join("design.tsv"))
        .arg("-n")
        .arg(root.join("result"))
        .args(["--norm-method", "none", "--threads", &threads.to_string()])
        .args(extra)
        .output()
        .unwrap()
}
const HEADER: &str = "sgRNA\tGene\ts0\ts1\ts2\ts3\n";
const NORMAL: &str = "n1\tNORMAL\t99\t119\t109\t129\nn2\tNORMAL\t149\t169\t159\t179\n";
const FAIL:&str="f1\tFAIL_FIRST\t0\t99999\t0\t999999\nf2\tFAIL_FIRST\t99999\t0\t999999\t0\nf3\tFAIL_SECOND\t0\t99999\t0\t999999\nf4\tFAIL_SECOND\t99999\t0\t999999\t0\n";
#[test]
fn planned_failure_publishes_only_complete_deterministic_diagnostics() {
    let mut baseline = None;
    let mut stderr = None;
    for threads in [1, 2, 4, 8] {
        let d = tempfile::tempdir().unwrap();
        let out = invoke(d.path(), &format!("{HEADER}{NORMAL}{FAIL}"), &[], threads);
        assert!(!out.status.success());
        assert!(!d.path().join("result.gene_summary.txt").exists());
        assert!(!d.path().join("result.sgrna_summary.txt").exists());
        let bytes = fs::read(d.path().join("result.mle_diagnostics.json")).unwrap();
        let text = String::from_utf8_lossy(&bytes);
        assert!(text.contains("FIT_FAILED"));
        assert!(text.contains("FAIL_SECOND"));
        assert!(text.contains("\"m\":0"));
        assert!(
            String::from_utf8_lossy(&out.stderr).contains("FAIL_FIRST")
                || String::from_utf8_lossy(&out.stderr).contains("converged eligible fit")
        );
        if let Some(b) = &baseline {
            assert_eq!(&bytes, b);
            assert_eq!(stderr.as_ref().unwrap(), &out.stderr);
        } else {
            baseline = Some(bytes);
            stderr = Some(out.stderr);
        }
    }
}
#[test]
fn control_skip_singleton_and_empty_universe() {
    let d = tempfile::tempdir().unwrap();
    let count=format!("{HEADER}{NORMAL}c1\tControl\t20\t20\t20\t20\ns1\tSKIP\t20\t20\t20\t20\ns2\tSKIP\t30\t30\t30\t30\ns3\tSKIP\t40\t40\t40\t40\n");
    let out = invoke(d.path(), &count, &["--max-sgrnas-per-gene", "2"], 2);
    assert!(out.status.success(), "{out:?}");
    let gene = fs::read_to_string(d.path().join("result.gene_summary.txt")).unwrap();
    assert!(gene.contains("NORMAL"));
    assert!(!gene.contains("Control") && !gene.contains("SKIP"));
    let diag = fs::read_to_string(d.path().join("result.mle_diagnostics.json")).unwrap();
    assert!(diag.contains("CONTROL_ONLY") && diag.contains("PREDECLARED_GUIDE_LIMIT"));
    assert!(diag.contains("\"m\":1"));
    let row: Vec<_> = gene.lines().nth(1).unwrap().split('\t').collect();
    assert_eq!(row[6], row[7]);
    let d = tempfile::tempdir().unwrap();
    let out = invoke(
        d.path(),
        &format!("{HEADER}c1\tControl\t20\t20\t20\t20\n"),
        &[],
        1,
    );
    assert!(!out.status.success());
    assert!(!d.path().join("result.gene_summary.txt").exists());
    assert!(
        fs::read_to_string(d.path().join("result.mle_diagnostics.json"))
            .unwrap()
            .contains("EMPTY_ELIGIBLE_UNIVERSE")
            || String::from_utf8_lossy(&out.stderr).contains("EMPTY_ELIGIBLE_UNIVERSE")
    );
}

#[test]
fn zero_residual_wald_is_ineligible_without_retry_or_primary_publication() {
    let d = tempfile::tempdir().unwrap();
    let out = invoke(d.path(), &format!("{HEADER}z1\tZERO\t0\t0\t0\t0\n"), &[], 1);
    assert!(!out.status.success());
    assert!(!d.path().join("result.gene_summary.txt").exists());
    let diag = fs::read_to_string(d.path().join("result.mle_diagnostics.json")).unwrap();
    assert!(diag.contains("WALD_INELIGIBLE"));
    assert!(diag.contains("INFERENCE_INVALID"));
    assert!(!diag.contains("\"relaxation\":0.5"));
    assert!(!diag.contains("1000000"));
}
