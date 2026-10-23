use std::{
    fs,
    process::{Command, Output},
};
fn invoke(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(args)
        .output()
        .unwrap()
}

#[test]
fn unmapped_preserves_ids_sequence_quality_and_group_file_record_order() {
    let d = tempfile::tempdir().unwrap();
    let l = d.path().join("lib");
    fs::write(
        &l,
        "id\tsequence\tgene\na\tAAAA\tA\nb\tAAAC\tB\nc\tCCCC\tC\n",
    )
    .unwrap();
    let a = d.path().join("a");
    let b = d.path().join("b");
    let c = d.path().join("c");
    fs::write(&a,"@exact description\nAAAA\n+original\n!#$%\n@ambiguous\nAAAG\n+\nABCD\n@short\nTT\n+\n!!\n@unique mismatch\nCCCT\n+\nHHHH\n").unwrap();
    fs::write(&b, "@no match\nGGGG\n+\nJKLM\n").unwrap();
    fs::write(&c, "@later sample\nTTTT\n+\nNOPQ\n").unwrap();
    let p = d.path().join("out").display().to_string();
    let group = format!("{},{}", a.display(), b.display());
    let o = invoke(&[
        "count",
        "-l",
        l.to_str().unwrap(),
        "-f",
        &group,
        c.to_str().unwrap(),
        "-n",
        &p,
        "--sgrna-len",
        "4",
        "--test-mismatch",
        "1",
        "--norm-method",
        "none",
        "--unmapped-to-file",
    ]);
    assert!(o.status.success(), "{o:?}");
    assert_eq!(fs::read_to_string(format!("{p}.unmapped.fastq")).unwrap(),"@ambiguous\nAAAG\n+\nABCD\n@short\nTT\n+\n!!\n@no match\nGGGG\n+\nJKLM\n@later sample\nTTTT\n+\nNOPQ\n");
    assert_eq!(
        fs::read_to_string(format!("{p}.count.txt")).unwrap(),
        "sgRNA\tGene\tsample1\tsample2\na\tA\t1\t0\nb\tB\t0\t0\nc\tC\t1\t0\n"
    );
}

#[test]
fn existing_output_or_dangling_link_is_never_overwritten() {
    let d = tempfile::tempdir().unwrap();
    let l = d.path().join("lib");
    fs::write(&l, "id\tsequence\tgene\na\tAAAA\tA\n").unwrap();
    let f = d.path().join("r");
    fs::write(&f, "").unwrap();
    let p = d.path().join("out").display().to_string();
    for (suffix, link) in [
        ("unmapped.fastq", false),
        ("count.txt", false),
        ("umi_summary.txt", true),
    ] {
        let dest = format!("{p}.{suffix}");
        if link {
            std::os::unix::fs::symlink("missing", &dest).unwrap();
        } else {
            fs::write(&dest, "sentinel").unwrap();
        }
        let o = invoke(&[
            "count",
            "-l",
            l.to_str().unwrap(),
            "-f",
            f.to_str().unwrap(),
            "-n",
            &p,
            "--sgrna-len",
            "4",
            "--norm-method",
            "none",
            "--unmapped-to-file",
            "--umi-start",
            "4",
            "--umi-length",
            "1",
            "--umi-dedup-scope",
            "file",
            "--umi-max-keys",
            "1",
        ]);
        assert_eq!(o.status.code(), Some(1));
        assert!(String::from_utf8(o.stderr)
            .unwrap()
            .contains("Output already exists"));
        if link {
            assert_eq!(
                fs::read_link(&dest).unwrap(),
                std::path::PathBuf::from("missing")
            );
        } else {
            assert_eq!(fs::read_to_string(&dest).unwrap(), "sentinel");
        }
        fs::remove_file(dest).unwrap();
        assert_eq!(fs::read_dir(d.path()).unwrap().count(), 2);
    }
}

#[test]
fn empty_and_gzip_inputs_stream_correctly() {
    let d = tempfile::tempdir().unwrap();
    let l = d.path().join("lib");
    fs::write(&l, "id\tsequence\tgene\na\tAAAA\tA\n").unwrap();
    let empty = d.path().join("empty");
    fs::write(&empty, "").unwrap();
    let gzip = d.path().join("r.gz");
    fs::write(
        &gzip,
        include_bytes!("fixtures/campaign_a/unmapped.fastq.gz"),
    )
    .unwrap();
    let group = format!("{},{}", empty.display(), gzip.display());
    let p = d.path().join("out").display().to_string();
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
        "--unmapped-to-file",
    ]);
    assert!(o.status.success(), "{o:?}");
    assert_eq!(
        fs::read_to_string(format!("{p}.unmapped.fastq")).unwrap(),
        "@unmapped gzip\nGGGG\n+\n!#$%\n"
    );
    assert!(fs::read_to_string(format!("{p}.countsummary.txt"))
        .unwrap()
        .contains(&format!("{}\tsample1\t0\t0\t0.0000\t1", empty.display())));
}

#[test]
fn malformed_or_fasta_input_and_missing_output_parent_leave_no_success() {
    for contents in [
        "@good\nGGGG\n+\nIIII\n@bad\nAAAA\n+\nII\n",
        ">fasta\nAAAA\n",
    ] {
        let d = tempfile::tempdir().unwrap();
        let l = d.path().join("lib");
        fs::write(&l, "id\tsequence\tgene\na\tAAAA\tA\n").unwrap();
        let f = d.path().join("r");
        fs::write(&f, contents).unwrap();
        let p = d.path().join("out").display().to_string();
        let o = invoke(&[
            "count",
            "-l",
            l.to_str().unwrap(),
            "-f",
            f.to_str().unwrap(),
            "-n",
            &p,
            "--sgrna-len",
            "4",
            "--norm-method",
            "none",
            "--unmapped-to-file",
        ]);
        assert_eq!(o.status.code(), Some(1));
        assert_eq!(fs::read_dir(d.path()).unwrap().count(), 2);
    }
    let d = tempfile::tempdir().unwrap();
    let l = d.path().join("lib");
    let f = d.path().join("r");
    fs::write(&l, "id\tsequence\tgene\na\tAAAA\tA\n").unwrap();
    fs::write(&f, "").unwrap();
    let p = d.path().join("missing/out").display().to_string();
    assert_eq!(
        invoke(&[
            "count",
            "-l",
            l.to_str().unwrap(),
            "-f",
            f.to_str().unwrap(),
            "-n",
            &p,
            "--sgrna-len",
            "4",
            "--unmapped-to-file"
        ])
        .status
        .code(),
        Some(1)
    );
    assert_eq!(fs::read_dir(d.path()).unwrap().count(), 2);
}
