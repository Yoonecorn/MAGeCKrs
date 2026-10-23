use std::{
    fs,
    path::Path,
    process::{Command, Output},
};

fn invoke(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(args)
        .output()
        .unwrap()
}
fn reads(path: &Path, sequences: &[&str]) {
    let mut text = String::new();
    for (i, seq) in sequences.iter().enumerate() {
        text.push_str(&format!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
    }
    fs::write(path, text).unwrap();
}
fn assert_no_output(dir: &Path) {
    assert!(fs::read_dir(dir).unwrap().all(|p| {
        let name = p.unwrap().file_name();
        let name = name.to_string_lossy();
        !name.starts_with("out") && !name.starts_with(".mageckrs")
    }));
}

#[test]
fn explicit_file_and_sample_scopes_preserve_read_qc_and_group_boundaries() {
    let d = tempfile::tempdir().unwrap();
    let lib = d.path().join("lib");
    fs::write(&lib, "id\tsequence\tgene\ng1\tACGT\tA\ng2\tTGCA\tB\n").unwrap();
    let r1 = d.path().join("r1.fq");
    let r2 = d.path().join("r2.fq");
    reads(&r1, &["AAACGT", "aaACGT", "CCACGT", "AATGCA"]);
    reads(&r2, &["AAACGT", "GGACGT"]);
    let grouped = format!("{},{}", r1.display(), r2.display());
    for (scope, expected, new, duplicate) in [("file", 4, 2, 0), ("sample", 3, 1, 1)] {
        let prefix = d.path().join(format!("out-{scope}")).display().to_string();
        let output = invoke(&[
            "count",
            "-l",
            lib.to_str().unwrap(),
            "-f",
            &grouped,
            r1.to_str().unwrap(),
            "-n",
            &prefix,
            "--sgrna-len",
            "4",
            "--trim-5",
            "2",
            "--norm-method",
            "none",
            "--umi-start",
            "0",
            "--umi-length",
            "2",
            "--umi-dedup-scope",
            scope,
            "--umi-max-keys",
            "4",
        ]);
        assert!(output.status.success(), "{output:?}");
        assert_eq!(
            fs::read_to_string(format!("{prefix}.count.txt")).unwrap(),
            format!("sgRNA\tGene\tsample1\tsample2\ng1\tA\t{expected}\t2\ng2\tB\t1\t1\n")
        );
        let qc = fs::read_to_string(format!("{prefix}.countsummary.txt")).unwrap();
        assert!(
            qc.contains(&format!("{}\tsample1\t4\t4\t1.0000\t0", r1.display())),
            "{qc}"
        );
        assert!(
            qc.contains(&format!("{}\tsample1\t2\t2\t1.0000\t1", r2.display())),
            "{qc}"
        );
        let umi = fs::read_to_string(format!("{prefix}.umi_summary.txt")).unwrap();
        assert!(umi.contains(&format!("{}\tsample1\t4\t3\t1", r1.display())));
        assert!(umi.contains(&format!("{}\tsample1\t2\t{new}\t{duplicate}", r2.display())));
        let meta = fs::read_to_string(format!("{prefix}.run_info.txt")).unwrap();
        assert!(meta.contains(&format!("umi_dedup_scope\t{scope}\n")));
        assert!(meta.contains("umi_cap_is_process_memory_bound\tfalse\n"));
        let norm = fs::read_to_string(format!("{prefix}.count.normalized.txt")).unwrap();
        assert!(norm.contains(&format!("g1\tA\t{expected}.0000\t2.0000")));
    }
}

#[test]
fn cap_is_per_active_scope_and_duplicates_are_legal_at_capacity() {
    let d = tempfile::tempdir().unwrap();
    let lib = d.path().join("lib");
    fs::write(&lib, "id\tsequence\tgene\ng1\tACGT\tA\n").unwrap();
    let a = d.path().join("a");
    let b = d.path().join("b");
    reads(&a, &["AAACGT", "AAACGT"]);
    reads(&b, &["CCACGT"]);
    let group = format!("{},{}", a.display(), b.display());
    for (scope, success) in [("file", true), ("sample", false)] {
        let prefix = d.path().join("out").display().to_string();
        let output = invoke(&[
            "count",
            "-l",
            lib.to_str().unwrap(),
            "-f",
            &group,
            "-n",
            &prefix,
            "--sgrna-len",
            "4",
            "--trim-5",
            "2",
            "--norm-method",
            "none",
            "--umi-start",
            "0",
            "--umi-length",
            "2",
            "--umi-dedup-scope",
            scope,
            "--umi-max-keys",
            "1",
        ]);
        assert_eq!(output.status.success(), success, "{output:?}");
        if success {
            assert!(fs::read_to_string(format!("{prefix}.count.txt"))
                .unwrap()
                .contains("g1\tA\t2\n"));
            for suffix in [
                "count.txt",
                "count.normalized.txt",
                "countsummary.txt",
                "run_info.txt",
                "umi_summary.txt",
            ] {
                fs::remove_file(format!("{prefix}.{suffix}")).unwrap();
            }
        } else {
            let err = String::from_utf8(output.stderr).unwrap();
            assert!(
                err.contains("cap 1 exceeded") && err.contains("record 1"),
                "{err}"
            );
            assert_no_output(d.path());
        }
    }
}

#[test]
fn packed_and_byte_lengths_normalize_case_and_keep_exact_distinctions() {
    for length in [1, 32, 33, 64] {
        let d = tempfile::tempdir().unwrap();
        let lib = d.path().join("lib");
        fs::write(&lib, "id\tsequence\tgene\ng1\tACGT\tA\n").unwrap();
        let input = d.path().join("r");
        let upper = format!("ACGT{}", "A".repeat(length));
        let lower = format!("ACGT{}", "a".repeat(length));
        let distinct = format!("ACGT{}C", "A".repeat(length - 1));
        reads(&input, &[&upper, &lower, &distinct]);
        let prefix = d.path().join("out").display().to_string();
        let length = length.to_string();
        let out = invoke(&[
            "count",
            "-l",
            lib.to_str().unwrap(),
            "-f",
            input.to_str().unwrap(),
            "-n",
            &prefix,
            "--sgrna-len",
            "4",
            "--norm-method",
            "none",
            "--umi-start",
            "4",
            "--umi-length",
            &length,
            "--umi-dedup-scope",
            "sample",
            "--umi-max-keys",
            "2",
        ]);
        assert!(out.status.success(), "{out:?}");
        assert!(fs::read_to_string(format!("{prefix}.count.txt"))
            .unwrap()
            .contains("g1\tA\t2\n"));
    }
}

#[test]
fn mapped_invalid_short_or_overlapping_umi_fails_without_publication() {
    for (seq, trim, start, length, needle) in [
        ("NNACGT", "2", "0", "2", "invalid UMI"),
        ("ACGT", "none", "4", "2", "too short"),
        ("ACGTAA", "none", "0", "2", "overlaps"),
        ("GGACGTAA", "auto", "2", "2", "overlaps"),
    ] {
        let d = tempfile::tempdir().unwrap();
        let lib = d.path().join("lib");
        fs::write(&lib, "id\tsequence\tgene\ng1\tACGT\tA\n").unwrap();
        let r = d.path().join("r");
        reads(&r, &[seq]);
        let prefix = d.path().join("out").display().to_string();
        let o = invoke(&[
            "count",
            "-l",
            lib.to_str().unwrap(),
            "-f",
            r.to_str().unwrap(),
            "-n",
            &prefix,
            "--sgrna-len",
            "4",
            "--trim-5",
            trim,
            "--norm-method",
            "none",
            "--umi-start",
            start,
            "--umi-length",
            length,
            "--umi-dedup-scope",
            "file",
            "--umi-max-keys",
            "10",
            "--unmapped-to-file",
        ]);
        assert_eq!(o.status.code(), Some(1));
        let e = String::from_utf8(o.stderr).unwrap();
        assert!(e.contains(needle) && e.contains("record 1"), "{e}");
        assert_no_output(d.path());
    }
}

#[test]
fn incomplete_scope_zero_cap_and_overflow_rejected_before_io() {
    let max = usize::MAX.to_string();
    for extra in [
        vec![
            "--umi-start",
            "0",
            "--umi-length",
            "2",
            "--umi-max-keys",
            "2",
        ],
        vec![
            "--umi-start",
            "0",
            "--umi-length",
            "0",
            "--umi-dedup-scope",
            "file",
            "--umi-max-keys",
            "1",
        ],
        vec![
            "--umi-start",
            "0",
            "--umi-length",
            "2",
            "--umi-dedup-scope",
            "file",
            "--umi-max-keys",
            "0",
        ],
        vec![
            "--umi-start",
            &max,
            "--umi-length",
            "2",
            "--umi-dedup-scope",
            "sample",
            "--umi-max-keys",
            "1",
        ],
    ] {
        let mut args = vec!["count", "-l", "missing", "-f", "missing", "-n", "unused"];
        args.extend(extra);
        let o = invoke(&args);
        assert_eq!(o.status.code(), Some(1));
        let e = String::from_utf8(o.stderr).unwrap();
        assert!(
            e.contains("UMI") && !e.contains("file does not exist"),
            "{e}"
        );
    }
}

#[test]
fn auto_disjoint_suffix_works_and_unmapped_invalid_umi_is_not_a_molecule() {
    let d = tempfile::tempdir().unwrap();
    let l = d.path().join("lib");
    fs::write(&l, "id\tsequence\tgene\ng\tACGT\tA\n").unwrap();
    let r = d.path().join("r");
    reads(&r, &["GGACGTaa", "GGACGTAA", "NNNNNNNN"]);
    let p = d.path().join("out").display().to_string();
    let o = invoke(&[
        "count",
        "-l",
        l.to_str().unwrap(),
        "-f",
        r.to_str().unwrap(),
        "-n",
        &p,
        "--sgrna-len",
        "4",
        "--trim-5",
        "auto",
        "--test-mismatch",
        "1",
        "--norm-method",
        "none",
        "--umi-start",
        "6",
        "--umi-length",
        "2",
        "--umi-dedup-scope",
        "sample",
        "--umi-max-keys",
        "1",
    ]);
    assert!(o.status.success(), "{o:?}");
    assert!(fs::read_to_string(format!("{p}.count.txt"))
        .unwrap()
        .contains("g\tA\t1\n"));
}

#[test]
fn sample_totals_are_order_independent_but_file_attribution_is_explicit() {
    let d = tempfile::tempdir().unwrap();
    let lib = d.path().join("lib");
    fs::write(&lib, "id\tsequence\tgene\ng1\tACGT\tA\n").unwrap();
    let a = d.path().join("a");
    let b = d.path().join("b");
    reads(&a, &["AAACGT", "CCACGT"]);
    reads(&b, &["AAACGT", "GGACGT"]);
    let mut results = Vec::new();
    let mut metadata = Vec::new();
    for (i, group) in [
        format!("{},{}", a.display(), b.display()),
        format!("{},{}", b.display(), a.display()),
    ]
    .iter()
    .enumerate()
    {
        let p = d.path().join(format!("out{i}")).display().to_string();
        let o = invoke(&[
            "count",
            "-l",
            lib.to_str().unwrap(),
            "-f",
            group,
            "-n",
            &p,
            "--sgrna-len",
            "4",
            "--trim-5",
            "2",
            "--norm-method",
            "none",
            "--umi-start",
            "0",
            "--umi-length",
            "2",
            "--umi-dedup-scope",
            "sample",
            "--umi-max-keys",
            "3",
        ]);
        assert!(o.status.success(), "{o:?}");
        results.push(fs::read(format!("{p}.count.txt")).unwrap());
        metadata.push(fs::read_to_string(format!("{p}.umi_summary.txt")).unwrap());
    }
    assert_eq!(results[0], results[1]);
    assert!(String::from_utf8(results[0].clone())
        .unwrap()
        .contains("g1\tA\t3\n"));
    assert!(metadata[0].contains(&format!("{}\tsample1\t2\t2\t0", a.display())));
    assert!(metadata[1].contains(&format!("{}\tsample1\t2\t1\t1", a.display())));
}

#[test]
fn unique_mismatch_is_deduplicated_while_ambiguous_reads_have_no_umi_observation() {
    let d = tempfile::tempdir().unwrap();
    let l = d.path().join("lib");
    fs::write(
        &l,
        "id\tsequence\tgene\na\tAAAA\tA\nb\tAAAC\tB\nc\tCCCC\tC\n",
    )
    .unwrap();
    let r = d.path().join("r");
    reads(&r, &["AACCCC", "AACCCT", "NNAAAG"]);
    let p = d.path().join("out").display().to_string();
    let o = invoke(&[
        "count",
        "-l",
        l.to_str().unwrap(),
        "-f",
        r.to_str().unwrap(),
        "-n",
        &p,
        "--sgrna-len",
        "4",
        "--trim-5",
        "2",
        "--test-mismatch",
        "1",
        "--norm-method",
        "none",
        "--umi-start",
        "0",
        "--umi-length",
        "2",
        "--umi-dedup-scope",
        "sample",
        "--umi-max-keys",
        "1",
        "--unmapped-to-file",
    ]);
    assert!(o.status.success(), "{o:?}");
    assert!(fs::read_to_string(format!("{p}.count.txt"))
        .unwrap()
        .contains("c\tC\t1\n"));
    assert_eq!(
        fs::read_to_string(format!("{p}.unmapped.fastq")).unwrap(),
        "@r2\nNNAAAG\n+\nIIIIII\n"
    );
}

#[test]
fn empty_physical_files_create_read_and_molecule_zero_rows() {
    let d = tempfile::tempdir().unwrap();
    let l = d.path().join("lib");
    fs::write(&l, "id\tsequence\tgene\ng\tACGT\tA\n").unwrap();
    let r = d.path().join("r");
    fs::write(&r, "").unwrap();
    let group = format!("{},{}", r.display(), r.display());
    for scope in ["file", "sample"] {
        let p = d.path().join(scope).display().to_string();
        let o = invoke(&[
            "count",
            "-l",
            l.to_str().unwrap(),
            "-f",
            &group,
            "-n",
            &p,
            "--sgrna-len",
            "4",
            "--norm-method",
            "none",
            "--umi-start",
            "4",
            "--umi-length",
            "1",
            "--umi-dedup-scope",
            scope,
            "--umi-max-keys",
            "1",
        ]);
        assert!(o.status.success(), "{o:?}");
        let m = fs::read_to_string(format!("{p}.umi_summary.txt")).unwrap();
        assert_eq!(m.lines().count(), 3);
        assert!(m.lines().skip(1).all(|l| l.ends_with("\t0\t0\t0")));
    }
}

#[test]
fn mixed_plain_gzip_replicates_apply_the_selected_scope() {
    let d = tempfile::tempdir().unwrap();
    let lib = d.path().join("lib");
    fs::write(&lib, "id\tsequence\tgene\ng1\tACGT\tA\n").unwrap();
    let gzip = d.path().join("r.gz");
    fs::write(&gzip, include_bytes!("fixtures/campaign_a/umi.fastq.gz")).unwrap();
    let plain = d.path().join("r.fq");
    reads(&plain, &["AAACGT"]);
    let group = format!("{},{}", gzip.display(), plain.display());
    for (scope, molecules) in [("file", 3), ("sample", 2)] {
        let p = d.path().join(scope).display().to_string();
        let o = invoke(&[
            "count",
            "-l",
            lib.to_str().unwrap(),
            "-f",
            &group,
            "-n",
            &p,
            "--sgrna-len",
            "4",
            "--trim-5",
            "2",
            "--norm-method",
            "none",
            "--umi-start",
            "0",
            "--umi-length",
            "2",
            "--umi-dedup-scope",
            scope,
            "--umi-max-keys",
            "2",
        ]);
        assert!(o.status.success(), "{o:?}");
        assert!(fs::read_to_string(format!("{p}.count.txt"))
            .unwrap()
            .contains(&format!("g1\tA\t{molecules}\n")));
    }
}
