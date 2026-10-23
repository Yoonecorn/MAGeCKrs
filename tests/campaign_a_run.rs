use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
};
fn invoke(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(args)
        .output()
        .unwrap()
}
fn success(args: &[&str]) {
    let o = invoke(args);
    assert!(o.status.success(), "{args:?}: {o:?}");
}
fn fixture(d: &Path) -> (PathBuf, Vec<String>) {
    let guides = ["AAAA", "CCCC", "GGGG", "TTTT", "ACGT", "TGCA"];
    let lib = d.join("library.tsv");
    let mut library = "id\tsequence\tgene\n".to_string();
    for (i, s) in guides.iter().enumerate() {
        library.push_str(&format!("g{}\t{s}\tGENE{}\n", i + 1, i / 2));
    }
    fs::write(&lib, library).unwrap();
    let mut paths = Vec::new();
    for (sample, counts) in [
        [5, 6, 4, 7, 5, 6],
        [5, 6, 4, 7, 5, 5],
        [11, 10, 9, 13, 12, 10],
        [2, 3, 30, 28, 10, 13],
    ]
    .iter()
    .enumerate()
    {
        let p = d.join(format!("sample-{sample}.fq"));
        let mut fq = String::new();
        for (g, &count) in counts.iter().enumerate() {
            for r in 0..count {
                fq.push_str(&format!("@{sample}-{g}-{r}\n{}AA\n+\nIIIIII\n", guides[g]));
            }
        }
        fs::write(&p, fq).unwrap();
        paths.push(p.display().to_string());
    }
    (
        lib,
        vec![
            format!("{},{}", paths[0], paths[1]),
            paths[2].clone(),
            paths[3].clone(),
        ],
    )
}
fn metadata(p: &str) -> String {
    fs::read_to_string(format!("{p}.run_info.txt")).unwrap()
}
fn no_success(d: &Path) {
    assert!(fs::read_dir(d).unwrap().all(|p| {
        let n = p.unwrap().file_name();
        let n = n.to_string_lossy();
        !n.starts_with("out") && !n.starts_with(".mageckrs")
    }));
}

#[test]
fn combined_matches_separate_count_then_test_for_all_normalizations_roles_and_seed_choices() {
    for norm in ["none", "median", "total", "control"] {
        for seed in [None, Some("0"), Some("18446744073709551615")] {
            let d = tempfile::tempdir().unwrap();
            let (lib, groups) = fixture(d.path());
            let controls = d.path().join("controls.txt");
            fs::write(&controls, "g5\ng6\n").unwrap();
            let sep = d.path().join("separate").display().to_string();
            let sep_test = d.path().join("separate-test").display().to_string();
            let combined = d.path().join("combined").display().to_string();
            let mut common = vec![
                "-l",
                lib.to_str().unwrap(),
                "-f",
                &groups[0],
                &groups[1],
                &groups[2],
                "--sample-label",
                "C1,C2,T",
                "--sgrna-len",
                "4",
                "--norm-method",
                norm,
            ];
            if norm == "control" {
                common.extend(["--control-sgrna", controls.to_str().unwrap()]);
            }
            let mut count = vec!["count"];
            count.extend(&common);
            count.extend(["-n", &sep]);
            success(&count);
            let count_path = format!("{sep}.count.txt");
            // Index roles exercise exactly the same existing resolver as label roles.
            let (control, treatment) = if seed.is_some() {
                ("0,1", "2")
            } else {
                ("C1,C2", "T")
            };
            let mut test = vec![
                "test",
                "-k",
                &count_path,
                "-c",
                control,
                "-t",
                treatment,
                "-n",
                &sep_test,
                "--norm-method",
                norm,
                "--permutation",
                "19",
            ];
            if norm == "control" {
                test.extend(["--control-sgrna", controls.to_str().unwrap()]);
            }
            if let Some(seed) = seed {
                test.extend(["--seed", seed]);
            }
            success(&test);
            let mut run = vec!["run"];
            run.extend(&common);
            run.extend([
                "-c",
                control,
                "-t",
                treatment,
                "-n",
                &combined,
                "--permutation",
                "19",
            ]);
            if let Some(seed) = seed {
                run.extend(["--seed", seed]);
            }
            success(&run);
            for suffix in ["count.txt", "count.normalized.txt", "countsummary.txt"] {
                assert_eq!(
                    fs::read(format!("{combined}.{suffix}")).unwrap(),
                    fs::read(format!("{sep}.{suffix}")).unwrap(),
                    "{norm} {seed:?} {suffix}"
                );
            }
            for suffix in ["gene_summary.txt", "sgrna_summary.txt"] {
                assert_eq!(
                    fs::read(format!("{combined}.{suffix}")).unwrap(),
                    fs::read(format!("{sep_test}.{suffix}")).unwrap(),
                    "{norm} {seed:?} {suffix}"
                );
            }
            assert_eq!(
                fs::read(format!("{combined}.count_run_info.txt")).unwrap(),
                fs::read(format!("{sep}.run_info.txt")).unwrap()
            );
            // Pin the intentionally different output_prefix; this is not full-file identity.
            let expected = metadata(&sep_test).replace(
                &format!("output_prefix\t{sep_test}\n"),
                &format!("output_prefix\t{combined}\n"),
            );
            assert_eq!(metadata(&combined), expected);
            assert!(!metadata(&combined).contains(".mageckrs-stage"));
        }
    }
}

#[test]
fn combined_umi_scopes_and_unmapped_match_separate_outputs() {
    for scope in ["file", "sample"] {
        let d = tempfile::tempdir().unwrap();
        let (lib, groups) = fixture(d.path());
        let sep = d.path().join("separate").display().to_string();
        let sep_test = d.path().join("separate-test").display().to_string();
        let runprefix = d.path().join("combined").display().to_string();
        let common = [
            "-l",
            lib.to_str().unwrap(),
            "-f",
            &groups[0],
            &groups[1],
            &groups[2],
            "--sample-label",
            "C1,C2,T",
            "--sgrna-len",
            "4",
            "--norm-method",
            "none",
            "--umi-start",
            "4",
            "--umi-length",
            "2",
            "--umi-dedup-scope",
            scope,
            "--umi-max-keys",
            "6",
            "--unmapped-to-file",
        ];
        let mut count = vec!["count"];
        count.extend(common);
        count.extend(["-n", &sep]);
        success(&count);
        let raw = format!("{sep}.count.txt");
        success(&[
            "test",
            "-k",
            &raw,
            "-c",
            "C1,C2",
            "-t",
            "T",
            "-n",
            &sep_test,
            "--norm-method",
            "none",
            "--permutation",
            "7",
        ]);
        let mut run = vec!["run"];
        run.extend(common);
        run.extend([
            "-c",
            "C1,C2",
            "-t",
            "T",
            "-n",
            &runprefix,
            "--permutation",
            "7",
        ]);
        success(&run);
        for suffix in [
            "count.txt",
            "count.normalized.txt",
            "countsummary.txt",
            "umi_summary.txt",
            "unmapped.fastq",
        ] {
            assert_eq!(
                fs::read(format!("{runprefix}.{suffix}")).unwrap(),
                fs::read(format!("{sep}.{suffix}")).unwrap()
            );
        }
        for suffix in ["gene_summary.txt", "sgrna_summary.txt"] {
            assert_eq!(
                fs::read(format!("{runprefix}.{suffix}")).unwrap(),
                fs::read(format!("{sep_test}.{suffix}")).unwrap()
            );
        }
        assert_eq!(
            fs::read(format!("{runprefix}.count_run_info.txt")).unwrap(),
            fs::read(format!("{sep}.run_info.txt")).unwrap()
        );
    }
}

#[test]
fn count_stage_test_stage_and_control_grammar_failures_publish_nothing() {
    let d = tempfile::tempdir().unwrap();
    let (lib, groups) = fixture(d.path());
    let p = d.path().join("out").display().to_string();
    let base = [
        "run",
        "-l",
        lib.to_str().unwrap(),
        "-f",
        &groups[0],
        &groups[1],
        &groups[2],
        "--sample-label",
        "C1,C2,T",
        "--sgrna-len",
        "4",
        "-n",
        &p,
        "--norm-method",
        "none",
        "--permutation",
        "7",
    ];
    let mut role = base.to_vec();
    role.extend(["-t", "missing", "-c", "C1,C2"]);
    let o = invoke(&role);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8(o.stderr).unwrap().contains("missing"));
    no_success(d.path());
    let mut umi = base.to_vec();
    umi.extend([
        "-t",
        "T",
        "-c",
        "C1,C2",
        "--umi-start",
        "0",
        "--umi-length",
        "2",
        "--umi-dedup-scope",
        "sample",
        "--umi-max-keys",
        "6",
    ]);
    let o = invoke(&umi);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8(o.stderr).unwrap().contains("overlaps"));
    no_success(d.path());
    let control = d.path().join("sequence-controls");
    fs::write(&control, "ACGT\nTGCA\n").unwrap();
    let mut grammar = base.to_vec();
    let pos = grammar.iter().position(|a| *a == "none").unwrap();
    grammar[pos] = "control";
    grammar.extend([
        "-t",
        "T",
        "-c",
        "C1,C2",
        "--control-sgrna",
        control.to_str().unwrap(),
    ]);
    let o = invoke(&grammar);
    assert_eq!(o.status.code(), Some(1));
    let error = String::from_utf8(o.stderr).unwrap();
    assert_eq!(error, "Configuration error: Every guide has at least one zero (or subnormal) count; cannot compute median-ratio size factors\n");
    no_success(d.path());
    fs::write(&groups[2], "@bad\nAAAA\n+\nI\n").unwrap();
    let mut bad = base.to_vec();
    bad.extend(["-t", "T", "-c", "C1,C2"]);
    let o = invoke(&bad);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8(o.stderr).unwrap().contains("FASTQ"));
    no_success(d.path());
}

#[test]
fn every_combined_destination_is_protected_including_dangling_links() {
    let d = tempfile::tempdir().unwrap();
    let (lib, groups) = fixture(d.path());
    let p = d.path().join("out").display().to_string();
    for suffix in [
        "count.txt",
        "count.normalized.txt",
        "countsummary.txt",
        "count_run_info.txt",
        "gene_summary.txt",
        "sgrna_summary.txt",
        "run_info.txt",
    ] {
        let dest = format!("{p}.{suffix}");
        std::os::unix::fs::symlink("missing", &dest).unwrap();
        let o = invoke(&[
            "run",
            "-l",
            lib.to_str().unwrap(),
            "-f",
            &groups[0],
            &groups[1],
            &groups[2],
            "--sample-label",
            "C1,C2,T",
            "--sgrna-len",
            "4",
            "-n",
            &p,
            "-t",
            "T",
            "-c",
            "C1,C2",
        ]);
        assert_eq!(o.status.code(), Some(1));
        assert!(String::from_utf8(o.stderr)
            .unwrap()
            .contains("Output already exists"));
        assert_eq!(fs::read_link(&dest).unwrap(), PathBuf::from("missing"));
        fs::remove_file(dest).unwrap();
        no_success(d.path());
    }
}

#[test]
fn documented_example_runs_end_to_end_and_repeats_exactly() {
    let d = tempfile::tempdir().unwrap();
    let e = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/campaign-a");
    let l = e.join("library.tsv");
    let paths = (0..4)
        .map(|i| e.join(format!("sample-{i}.fastq")).display().to_string())
        .collect::<Vec<_>>();
    let group = format!("{},{}", paths[0], paths[1]);
    for trim in ["auto", "0"] {
        let sep = d.path().join(format!("count-{trim}")).display().to_string();
        let test = d.path().join(format!("test-{trim}")).display().to_string();
        let common = [
            "-l",
            l.to_str().unwrap(),
            "-f",
            &group,
            &paths[2],
            &paths[3],
            "--sample-label",
            "C1,C2,T",
            "--sgrna-len",
            "4",
            "--trim-5",
            trim,
            "--norm-method",
            "none",
        ];
        let mut count = vec!["count"];
        count.extend(common);
        count.extend(["-n", &sep]);
        success(&count);
        let raw = format!("{sep}.count.txt");
        success(&[
            "test",
            "-k",
            &raw,
            "-c",
            "C1,C2",
            "-t",
            "T",
            "-n",
            &test,
            "--norm-method",
            "none",
            "--permutation",
            "19",
            "--seed",
            "17",
        ]);
        let mut earlier = None;
        for repeat in 0..2 {
            let p = d
                .path()
                .join(format!("combined-{trim}-{repeat}"))
                .display()
                .to_string();
            let mut run = vec!["run"];
            run.extend(common);
            run.extend([
                "-n",
                &p,
                "-c",
                "C1,C2",
                "-t",
                "T",
                "--permutation",
                "19",
                "--seed",
                "17",
            ]);
            success(&run);
            let mut files = Vec::new();
            for suffix in [
                "count.txt",
                "count.normalized.txt",
                "countsummary.txt",
                "gene_summary.txt",
                "sgrna_summary.txt",
            ] {
                let bytes = fs::read(format!("{p}.{suffix}")).unwrap();
                let separate = if suffix.contains("summary.txt") && !suffix.contains("countsummary")
                {
                    &test
                } else {
                    &sep
                };
                assert_eq!(bytes, fs::read(format!("{separate}.{suffix}")).unwrap());
                files.push(bytes);
            }
            if let Some(prior) = earlier {
                assert_eq!(files, prior);
            }
            earlier = Some(files);
        }
    }
    let p = d.path().join("molecules").display().to_string();
    success(&[
        "count",
        "-l",
        l.to_str().unwrap(),
        "-f",
        &group,
        &paths[2],
        &paths[3],
        "--sample-label",
        "C1,C2,T",
        "--sgrna-len",
        "4",
        "--norm-method",
        "none",
        "--umi-start",
        "4",
        "--umi-length",
        "2",
        "--umi-dedup-scope",
        "sample",
        "--umi-max-keys",
        "6",
        "--unmapped-to-file",
        "-n",
        &p,
    ]);
    assert_eq!(fs::read(format!("{p}.unmapped.fastq")).unwrap(), b"");
    assert!(fs::read_to_string(format!("{p}.count.txt"))
        .unwrap()
        .lines()
        .skip(1)
        .all(|l| l.ends_with("\t1\t1\t1")));
}
