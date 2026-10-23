use clap::Parser;
use mageckrs::cli::{Cli, Commands, CountArgs, NormMethod};
use mageckrs::error::MageckError;
use mageckrs::run_with_cli;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

/// Helper to write text content to a file.
fn create_file(path: &Path, content: &str) {
    fs::write(path, content).unwrap_or_else(|e| {
        panic!("Failed to write to file '{}': {e}", path.display());
    });
}

/// Helper to create a standard 2-guide reference library in TSV format.
fn create_standard_library(path: &Path) {
    let content = "id\tsequence\tgene\n\
        g1\tAAAAAAAAAAAAAAAAAAAA\tGENE_A\n\
        g2\tCCCCCCCCCCCCCCCCCCCC\tGENE_B\n";
    create_file(path, content);
}

/// Helper to create a 3-guide reference library in TSV format.
fn create_three_guide_library(path: &Path) {
    let content = "id\tsequence\tgene\n\
        g1\tAAAAAAAAAAAAAAAAAAAA\tGENE_A\n\
        g2\tCCCCCCCCCCCCCCCCCCCC\tGENE_B\n\
        g3\tGGGGGGGGGGGGGGGGGGGG\tGENE_C\n";
    create_file(path, content);
}

/// Helper to construct default `CountArgs` for integration testing.
fn default_count_args(lib: &Path, fastq: Vec<String>, prefix: &str) -> CountArgs {
    CountArgs {
        list_seq: lib.to_path_buf(),
        fastq,
        output_prefix: prefix.to_string(),
        sample_label: None,
        norm_method: NormMethod::None,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    }
}

/// Helper to verify that none of the three count pipeline output files exist on disk.
fn assert_no_output_files(prefix: &str) {
    let count_file = format!("{prefix}.count.txt");
    let norm_file = format!("{prefix}.count.normalized.txt");
    let summary_file = format!("{prefix}.countsummary.txt");

    assert!(
        !Path::new(&count_file).exists(),
        "Expected raw count file '{count_file}' NOT to exist on error path"
    );
    assert!(
        !Path::new(&norm_file).exists(),
        "Expected normalized count file '{norm_file}' NOT to exist on error path"
    );
    assert!(
        !Path::new(&summary_file).exists(),
        "Expected count summary file '{summary_file}' NOT to exist on error path"
    );
}

// -----------------------------------------------------------------------------
// 1. Single sample group, single FASTQ — headers and content
// -----------------------------------------------------------------------------
#[test]
fn test_count_single_sample_single_fastq_headers_and_content() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    // FASTQ: 2 reads for g1, 1 read for g2, 1 unmapped read
    let fq_path = dir.path().join("sample1.fastq");
    create_file(
        &fq_path,
        "@read1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @read2\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @read3\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @read4\nTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let mut args = default_count_args(
        &lib_path,
        vec![fq_path.to_string_lossy().to_string()],
        &prefix,
    );
    args.norm_method = NormMethod::None;

    let cli = Cli {
        command: Commands::Count(args),
    };
    run_with_cli(cli).expect("run_with_cli count single sample");

    let count_file = format!("{prefix}.count.txt");
    let norm_file = format!("{prefix}.count.normalized.txt");
    let summary_file = format!("{prefix}.countsummary.txt");

    assert!(Path::new(&count_file).exists());
    assert!(Path::new(&norm_file).exists());
    assert!(Path::new(&summary_file).exists());

    // Verify .count.txt header and content
    let count_content = fs::read_to_string(&count_file).expect("read count.txt");
    assert_eq!(
        count_content,
        "sgRNA\tGene\tsample1\ng1\tGENE_A\t2\ng2\tGENE_B\t1\n"
    );

    // Verify .count.normalized.txt header and content (None norm preserves raw values to 4 decimal places)
    let norm_content = fs::read_to_string(&norm_file).expect("read count.normalized.txt");
    assert_eq!(
        norm_content,
        "sgRNA\tGene\tsample1\ng1\tGENE_A\t2.0000\ng2\tGENE_B\t1.0000\n"
    );

    // Verify .countsummary.txt header and content
    let summary_content = fs::read_to_string(&summary_file).expect("read countsummary.txt");
    let fq_display = fq_path.to_string_lossy();
    assert_eq!(
        summary_content,
        format!(
            "File\tLabel\tReads\tMapped\tPercentage\tZerocounts\n\
             {fq_display}\tsample1\t4\t3\t0.7500\t0\n"
        )
    );
}

// -----------------------------------------------------------------------------
// 2. Multiple sample groups — verify one column per group in count table
// -----------------------------------------------------------------------------
#[test]
fn test_count_multiple_sample_groups_columns() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    // Sample group 1 (ctrl): 2 reads g1, 0 reads g2
    let fq1 = dir.path().join("ctrl.fastq");
    create_file(
        &fq1,
        "@r1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r2\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    // Sample group 2 (treat): 1 read g1, 3 reads g2
    let fq2 = dir.path().join("treat.fastq");
    create_file(
        &fq2,
        "@r1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r2\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r3\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r4\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let mut args = default_count_args(
        &lib_path,
        vec![
            fq1.to_string_lossy().to_string(),
            fq2.to_string_lossy().to_string(),
        ],
        &prefix,
    );
    args.sample_label = Some(vec!["ctrl".to_string(), "treat".to_string()]);
    args.norm_method = NormMethod::None;

    let cli = Cli {
        command: Commands::Count(args),
    };
    run_with_cli(cli).expect("run_with_cli count multiple sample groups");

    let count_file = format!("{prefix}.count.txt");
    let norm_file = format!("{prefix}.count.normalized.txt");
    let summary_file = format!("{prefix}.countsummary.txt");

    let count_content = fs::read_to_string(&count_file).expect("read count.txt");
    assert_eq!(
        count_content,
        "sgRNA\tGene\tctrl\ttreat\ng1\tGENE_A\t2\t1\ng2\tGENE_B\t0\t3\n"
    );

    let norm_content = fs::read_to_string(&norm_file).expect("read count.normalized.txt");
    assert_eq!(
        norm_content,
        "sgRNA\tGene\tctrl\ttreat\ng1\tGENE_A\t2.0000\t1.0000\ng2\tGENE_B\t0.0000\t3.0000\n"
    );

    let summary_content = fs::read_to_string(&summary_file).expect("read countsummary.txt");
    let f1 = fq1.to_string_lossy();
    let f2 = fq2.to_string_lossy();
    assert_eq!(
        summary_content,
        format!(
            "File\tLabel\tReads\tMapped\tPercentage\tZerocounts\n\
             {f1}\tctrl\t2\t2\t1.0000\t1\n\
             {f2}\ttreat\t4\t4\t1.0000\t0\n"
        )
    );
}

// -----------------------------------------------------------------------------
// 3. Technical replicates (comma-separated) — summed counts, per-file summary
// -----------------------------------------------------------------------------
#[test]
fn test_count_technical_replicates_summed_and_summarized() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    // Replicate 1: 2 reads g1, 0 reads g2
    let rep1 = dir.path().join("rep1.fastq");
    create_file(
        &rep1,
        "@r1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r2\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    // Replicate 2: 3 reads g1, 4 reads g2
    let rep2 = dir.path().join("rep2.fastq");
    create_file(
        &rep2,
        "@r1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r2\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r3\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r4\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r5\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r6\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r7\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let comma_group = format!("{},{}", rep1.to_string_lossy(), rep2.to_string_lossy());

    let mut args = default_count_args(&lib_path, vec![comma_group], &prefix);
    args.sample_label = Some(vec!["merged_sample".to_string()]);
    args.norm_method = NormMethod::None;

    let cli = Cli {
        command: Commands::Count(args),
    };
    run_with_cli(cli).expect("run_with_cli technical replicates");

    let count_file = format!("{prefix}.count.txt");
    let summary_file = format!("{prefix}.countsummary.txt");

    // Single sample column in count table; counts are summed (2+3=5, 0+4=4)
    let count_content = fs::read_to_string(&count_file).expect("read count.txt");
    assert_eq!(
        count_content,
        "sgRNA\tGene\tmerged_sample\ng1\tGENE_A\t5\ng2\tGENE_B\t4\n"
    );

    // One summary row per physical file
    let summary_content = fs::read_to_string(&summary_file).expect("read countsummary.txt");
    let r1 = rep1.to_string_lossy();
    let r2 = rep2.to_string_lossy();
    assert_eq!(
        summary_content,
        format!(
            "File\tLabel\tReads\tMapped\tPercentage\tZerocounts\n\
             {r1}\tmerged_sample\t2\t2\t1.0000\t1\n\
             {r2}\tmerged_sample\t7\t7\t1.0000\t0\n"
        )
    );
}

// -----------------------------------------------------------------------------
// 4. Median normalization — four decimal places
// -----------------------------------------------------------------------------
#[test]
fn test_count_median_normalization_precision() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_three_guide_library(&lib_path);

    // Sample 1: 10 of g1, 20 of g2, 30 of g3
    let mut fq1_content = String::new();
    for i in 0..10 {
        fq1_content.push_str(&format!(
            "@g1_{i}\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    for i in 0..20 {
        fq1_content.push_str(&format!(
            "@g2_{i}\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    for i in 0..30 {
        fq1_content.push_str(&format!(
            "@g3_{i}\nGGGGGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    let fq1 = dir.path().join("sample1.fastq");
    create_file(&fq1, &fq1_content);

    // Sample 2: 20 of g1, 40 of g2, 60 of g3 (exactly 2x sample 1)
    let mut fq2_content = String::new();
    for i in 0..20 {
        fq2_content.push_str(&format!(
            "@g1_{i}\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    for i in 0..40 {
        fq2_content.push_str(&format!(
            "@g2_{i}\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    for i in 0..60 {
        fq2_content.push_str(&format!(
            "@g3_{i}\nGGGGGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    let fq2 = dir.path().join("sample2.fastq");
    create_file(&fq2, &fq2_content);

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let mut args = default_count_args(
        &lib_path,
        vec![
            fq1.to_string_lossy().to_string(),
            fq2.to_string_lossy().to_string(),
        ],
        &prefix,
    );
    args.sample_label = Some(vec!["s1".to_string(), "s2".to_string()]);
    args.norm_method = NormMethod::Median;

    let cli = Cli {
        command: Commands::Count(args),
    };
    run_with_cli(cli).expect("run_with_cli median normalization");

    let norm_file = format!("{prefix}.count.normalized.txt");
    let norm_content = fs::read_to_string(&norm_file).expect("read count.normalized.txt");

    // Theoretical normalized counts:
    // gm_g1 = sqrt(10*20) = 14.1421356...
    // gm_g2 = sqrt(20*40) = 28.2842712...
    // gm_g3 = sqrt(30*60) = 42.4264068...
    // Ratios s1: 1/sqrt(2) ≈ 0.70710678. Median sf1 = 1/sqrt(2).
    // Ratios s2: sqrt(2) ≈ 1.41421356. Median sf2 = sqrt(2).
    // Normalized counts: 10 / (1/sqrt(2)) = 14.1421; 20 / (1/sqrt(2)) = 28.2843; 30 / (1/sqrt(2)) = 42.4264.
    // Sample 2 normalized counts are identical.
    let expected = "sgRNA\tGene\ts1\ts2\n\
                    g1\tGENE_A\t14.1421\t14.1421\n\
                    g2\tGENE_B\t28.2843\t28.2843\n\
                    g3\tGENE_C\t42.4264\t42.4264\n";
    assert_eq!(norm_content, expected);

    // Verify each value has exactly four decimal places
    for line in norm_content.lines().skip(1) {
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts.len(), 4);
        for &val in &parts[2..] {
            let dot_pos = val.find('.').expect("expected decimal dot");
            assert_eq!(
                val[dot_pos + 1..].len(),
                4,
                "expected 4 decimal places: {val}"
            );
        }
    }
}

// -----------------------------------------------------------------------------
// 5. Total normalization — four decimal places
// -----------------------------------------------------------------------------
#[test]
fn test_count_total_normalization_precision() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    // Sample 1: 10 g1, 30 g2 (total = 40)
    let mut fq1_content = String::new();
    for i in 0..10 {
        fq1_content.push_str(&format!(
            "@g1_{i}\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    for i in 0..30 {
        fq1_content.push_str(&format!(
            "@g2_{i}\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    let fq1 = dir.path().join("sample1.fastq");
    create_file(&fq1, &fq1_content);

    // Sample 2: 20 g1, 60 g2 (total = 80)
    let mut fq2_content = String::new();
    for i in 0..20 {
        fq2_content.push_str(&format!(
            "@g1_{i}\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    for i in 0..60 {
        fq2_content.push_str(&format!(
            "@g2_{i}\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    let fq2 = dir.path().join("sample2.fastq");
    create_file(&fq2, &fq2_content);

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let mut args = default_count_args(
        &lib_path,
        vec![
            fq1.to_string_lossy().to_string(),
            fq2.to_string_lossy().to_string(),
        ],
        &prefix,
    );
    args.sample_label = Some(vec!["s1".to_string(), "s2".to_string()]);
    args.norm_method = NormMethod::Total;

    let cli = Cli {
        command: Commands::Count(args),
    };
    run_with_cli(cli).expect("run_with_cli total normalization");

    let norm_file = format!("{prefix}.count.normalized.txt");
    let norm_content = fs::read_to_string(&norm_file).expect("read count.normalized.txt");

    // Mean total = (40 + 80) / 2 = 60.
    // Size factors: s1 = 40/60 = 2/3; s2 = 80/60 = 4/3.
    // Normalized counts:
    // s1: 10 / (2/3) = 15.0000; 30 / (2/3) = 45.0000.
    // s2: 20 / (4/3) = 15.0000; 60 / (4/3) = 45.0000.
    let expected = "sgRNA\tGene\ts1\ts2\n\
                    g1\tGENE_A\t15.0000\t15.0000\n\
                    g2\tGENE_B\t45.0000\t45.0000\n";
    assert_eq!(norm_content, expected);

    // Verify each value has exactly four decimal places
    for line in norm_content.lines().skip(1) {
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts.len(), 4);
        for &val in &parts[2..] {
            let dot_pos = val.find('.').expect("expected decimal dot");
            assert_eq!(
                val[dot_pos + 1..].len(),
                4,
                "expected 4 decimal places: {val}"
            );
        }
    }
}

// -----------------------------------------------------------------------------
// 6. None normalization — raw values preserved in normalized output
// -----------------------------------------------------------------------------
#[test]
fn test_count_none_normalization_preserves_raw() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    // 7 reads of g1, 13 reads of g2
    let mut fq_content = String::new();
    for i in 0..7 {
        fq_content.push_str(&format!(
            "@g1_{i}\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    for i in 0..13 {
        fq_content.push_str(&format!(
            "@g2_{i}\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    let fq = dir.path().join("sample.fastq");
    create_file(&fq, &fq_content);

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let mut args = default_count_args(&lib_path, vec![fq.to_string_lossy().to_string()], &prefix);
    args.norm_method = NormMethod::None;
    args.sample_label = Some(vec!["raw_sample".to_string()]);

    let cli = Cli {
        command: Commands::Count(args),
    };
    run_with_cli(cli).expect("run_with_cli none normalization");

    let count_file = format!("{prefix}.count.txt");
    let norm_file = format!("{prefix}.count.normalized.txt");

    let count_content = fs::read_to_string(&count_file).expect("read count.txt");
    assert_eq!(
        count_content,
        "sgRNA\tGene\traw_sample\ng1\tGENE_A\t7\ng2\tGENE_B\t13\n"
    );

    let norm_content = fs::read_to_string(&norm_file).expect("read count.normalized.txt");
    assert_eq!(
        norm_content,
        "sgRNA\tGene\traw_sample\ng1\tGENE_A\t7.0000\ng2\tGENE_B\t13.0000\n"
    );
}

// -----------------------------------------------------------------------------
// 7. Control normalization with valid control file — run succeeds
// -----------------------------------------------------------------------------
#[test]
fn test_count_control_normalization_with_valid_control_file() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    // Control file pointing to g1
    let ctrl_path = dir.path().join("control_guides.txt");
    create_file(&ctrl_path, "g1\n");

    let mut fq_content = String::new();
    for i in 0..5 {
        fq_content.push_str(&format!(
            "@g1_{i}\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    for i in 0..10 {
        fq_content.push_str(&format!(
            "@g2_{i}\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    let fq = dir.path().join("sample.fastq");
    create_file(&fq, &fq_content);

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let mut args = default_count_args(&lib_path, vec![fq.to_string_lossy().to_string()], &prefix);
    args.norm_method = NormMethod::Control;
    args.control_sgrna = Some(ctrl_path);
    args.sample_label = Some(vec!["ctrl_norm_sample".to_string()]);

    let cli = Cli {
        command: Commands::Count(args),
    };
    let res = run_with_cli(cli);
    assert!(res.is_ok(), "Expected control normalization to succeed");

    let count_file = format!("{prefix}.count.txt");
    let norm_file = format!("{prefix}.count.normalized.txt");
    let summary_file = format!("{prefix}.countsummary.txt");

    assert!(Path::new(&count_file).exists());
    assert!(Path::new(&norm_file).exists());
    assert!(Path::new(&summary_file).exists());
}

// -----------------------------------------------------------------------------
// 8. Control normalization without control file — MageckError::Config before output
// -----------------------------------------------------------------------------
#[test]
fn test_count_control_norm_without_control_file_errors_before_output() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    let fq = dir.path().join("sample.fastq");
    create_file(&fq, "@r1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n");

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let mut args = default_count_args(&lib_path, vec![fq.to_string_lossy().to_string()], &prefix);
    args.norm_method = NormMethod::Control;
    args.control_sgrna = None;

    let cli = Cli {
        command: Commands::Count(args),
    };
    let res = run_with_cli(cli);

    match res {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("requires --control-sgrna")
                    || msg.contains("requires a control-sgrna file"),
                "Unexpected error message: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got {other:?}"),
    }

    assert_no_output_files(&prefix);
}

// -----------------------------------------------------------------------------
// 9. Missing library file — MageckError returned before output
// -----------------------------------------------------------------------------
#[test]
fn test_count_missing_library_file_errors_before_output() {
    let dir = TempDir::new().expect("tempdir");
    let missing_lib = dir.path().join("nonexistent_lib.tsv");

    let fq = dir.path().join("sample.fastq");
    create_file(&fq, "@r1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n");

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let args = default_count_args(
        &missing_lib,
        vec![fq.to_string_lossy().to_string()],
        &prefix,
    );

    let cli = Cli {
        command: Commands::Count(args),
    };
    let res = run_with_cli(cli);

    match res {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("file does not exist") || msg.contains("does not exist"),
                "Unexpected error message: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got {other:?}"),
    }

    assert_no_output_files(&prefix);
}

// -----------------------------------------------------------------------------
// 10. Missing FASTQ file — MageckError returned before output
// -----------------------------------------------------------------------------
#[test]
fn test_count_missing_fastq_file_errors_before_output() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    let missing_fq = dir.path().join("nonexistent.fastq");

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let args = default_count_args(
        &lib_path,
        vec![missing_fq.to_string_lossy().to_string()],
        &prefix,
    );

    let cli = Cli {
        command: Commands::Count(args),
    };
    let res = run_with_cli(cli);

    match res {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("file does not exist") || msg.contains("does not exist"),
                "Unexpected error message: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got {other:?}"),
    }

    assert_no_output_files(&prefix);
}

// -----------------------------------------------------------------------------
// 11. Empty FASTQ — zero summary row written correctly
// -----------------------------------------------------------------------------
#[test]
fn test_count_empty_fastq_zero_summary_row() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    // Completely empty FASTQ file
    let empty_fq = dir.path().join("empty.fastq");
    create_file(&empty_fq, "");

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let mut args = default_count_args(
        &lib_path,
        vec![empty_fq.to_string_lossy().to_string()],
        &prefix,
    );
    args.norm_method = NormMethod::None;
    args.sample_label = Some(vec!["empty_sample".to_string()]);

    let cli = Cli {
        command: Commands::Count(args),
    };
    run_with_cli(cli).expect("run_with_cli empty fastq");

    let count_file = format!("{prefix}.count.txt");
    let summary_file = format!("{prefix}.countsummary.txt");

    let count_content = fs::read_to_string(&count_file).expect("read count.txt");
    assert_eq!(
        count_content,
        "sgRNA\tGene\tempty_sample\ng1\tGENE_A\t0\ng2\tGENE_B\t0\n"
    );

    let summary_content = fs::read_to_string(&summary_file).expect("read countsummary.txt");
    let fq_display = empty_fq.to_string_lossy();
    // Reads=0, Mapped=0, Percentage=0.0000, Zerocounts=2 (all 2 guides in library are zero)
    assert_eq!(
        summary_content,
        format!(
            "File\tLabel\tReads\tMapped\tPercentage\tZerocounts\n\
             {fq_display}\tempty_sample\t0\t0\t0.0000\t2\n"
        )
    );
}

// -----------------------------------------------------------------------------
// 12. One-mismatch mode enabled — matched reads are counted
// -----------------------------------------------------------------------------
#[test]
fn test_count_one_mismatch_mode_counts_matches() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    // Read has 1 mismatch against g1 (first base is T instead of A)
    let fq = dir.path().join("mismatch.fastq");
    create_file(
        &fq,
        "@read_1mm\nTAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    // Run 1: test_mismatch = 0 (exact match only). Read must NOT be mapped.
    let prefix0 = dir.path().join("out_exact").to_string_lossy().to_string();
    let mut args0 = default_count_args(&lib_path, vec![fq.to_string_lossy().to_string()], &prefix0);
    args0.test_mismatch = 0;
    args0.norm_method = NormMethod::None;
    run_with_cli(Cli {
        command: Commands::Count(args0),
    })
    .expect("run exact match");

    let count_exact = fs::read_to_string(format!("{prefix0}.count.txt")).expect("read count");
    assert_eq!(
        count_exact,
        "sgRNA\tGene\tsample1\ng1\tGENE_A\t0\ng2\tGENE_B\t0\n"
    );

    let summary_exact =
        fs::read_to_string(format!("{prefix0}.countsummary.txt")).expect("read sum");
    let fq_display = fq.to_string_lossy();
    assert!(
        summary_exact.contains(&format!("{fq_display}\tsample1\t1\t0\t0.0000\t2")),
        "Expected 0 mapped reads with test_mismatch=0: {summary_exact}"
    );

    // Run 2: test_mismatch = 1. Read MUST be mapped.
    let prefix1 = dir.path().join("out_1mm").to_string_lossy().to_string();
    let mut args1 = default_count_args(&lib_path, vec![fq.to_string_lossy().to_string()], &prefix1);
    args1.test_mismatch = 1;
    args1.norm_method = NormMethod::None;
    run_with_cli(Cli {
        command: Commands::Count(args1),
    })
    .expect("run 1-mismatch");

    let count_1mm = fs::read_to_string(format!("{prefix1}.count.txt")).expect("read count");
    assert_eq!(
        count_1mm,
        "sgRNA\tGene\tsample1\ng1\tGENE_A\t1\ng2\tGENE_B\t0\n"
    );

    let summary_1mm = fs::read_to_string(format!("{prefix1}.countsummary.txt")).expect("read sum");
    assert!(
        summary_1mm.contains(&format!("{fq_display}\tsample1\t1\t1\t1.0000\t1")),
        "Expected 1 mapped read with test_mismatch=1: {summary_1mm}"
    );
}

// -----------------------------------------------------------------------------
// 13. Verify exact six-column countsummary.txt header
// -----------------------------------------------------------------------------
#[test]
fn test_countsummary_exact_six_column_header() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    let fq = dir.path().join("sample.fastq");
    create_file(&fq, "@r1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n");

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let args = default_count_args(&lib_path, vec![fq.to_string_lossy().to_string()], &prefix);

    run_with_cli(Cli {
        command: Commands::Count(args),
    })
    .expect("run count");

    let summary_path = format!("{prefix}.countsummary.txt");
    let content = fs::read_to_string(&summary_path).expect("read summary");
    let first_line = content.lines().next().expect("first line");

    assert_eq!(
        first_line, "File\tLabel\tReads\tMapped\tPercentage\tZerocounts",
        "Header must match exact six-column contract"
    );

    let cols: Vec<&str> = first_line.split('\t').collect();
    assert_eq!(
        cols.len(),
        6,
        "Expected exactly 6 columns in countsummary header"
    );
    assert_eq!(
        cols,
        vec![
            "File",
            "Label",
            "Reads",
            "Mapped",
            "Percentage",
            "Zerocounts"
        ]
    );
}

// -----------------------------------------------------------------------------
// 14. Verify Percentage is mapped fraction [0,1] with exactly four decimal places
// -----------------------------------------------------------------------------
#[test]
fn test_countsummary_percentage_mapped_fraction_scale_and_precision() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    // 4 reads: 1 mapped to g1, 3 unmapped
    let fq = dir.path().join("frac.fastq");
    create_file(
        &fq,
        "@r1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r2\nTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r3\nGGGGGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r4\nNNNNNNNNNNNNNNNNNNNN\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    let prefix = dir.path().join("out").to_string_lossy().to_string();
    let args = default_count_args(&lib_path, vec![fq.to_string_lossy().to_string()], &prefix);

    run_with_cli(Cli {
        command: Commands::Count(args),
    })
    .expect("run count");

    let summary_path = format!("{prefix}.countsummary.txt");
    let content = fs::read_to_string(&summary_path).expect("read summary");
    let data_line = content.lines().nth(1).expect("data line");
    let cols: Vec<&str> = data_line.split('\t').collect();
    assert_eq!(cols.len(), 6);

    let reads: usize = cols[2].parse().expect("parse Reads");
    let mapped: usize = cols[3].parse().expect("parse Mapped");
    let pct_str = cols[4];

    assert_eq!(reads, 4);
    assert_eq!(mapped, 1);
    // Fraction: 1 / 4 = 0.2500 (NOT 25.0000 or 25%)
    assert_eq!(
        pct_str, "0.2500",
        "Percentage column must be mapped fraction in [0,1] formatted to {{:.4}}"
    );

    let pct_f64: f64 = pct_str.parse().expect("parse Percentage as f64");
    assert!(
        (0.0..=1.0).contains(&pct_f64),
        "Percentage must be in [0.0, 1.0], got {pct_f64}"
    );
    let dot_pos = pct_str.find('.').expect("expected decimal point");
    assert_eq!(pct_str[dot_pos + 1..].len(), 4);
}

// -----------------------------------------------------------------------------
// 15. Verify no output files are created on error paths
// -----------------------------------------------------------------------------
#[test]
fn test_count_no_output_files_on_error_paths() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    let fq_path = dir.path().join("sample.fastq");
    create_file(
        &fq_path,
        "@r1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    // Subtest A: Empty output prefix
    {
        let mut args =
            default_count_args(&lib_path, vec![fq_path.to_string_lossy().to_string()], "");
        args.output_prefix = "".to_string();
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files("");
    }

    // Subtest B: Missing library file
    {
        let missing_lib = dir.path().join("missing.tsv");
        let prefix = dir.path().join("err_lib").to_string_lossy().to_string();
        let args = default_count_args(
            &missing_lib,
            vec![fq_path.to_string_lossy().to_string()],
            &prefix,
        );
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files(&prefix);
    }

    // Subtest C: Missing FASTQ file
    {
        let missing_fq = dir.path().join("missing.fastq");
        let prefix = dir.path().join("err_fq").to_string_lossy().to_string();
        let args = default_count_args(
            &lib_path,
            vec![missing_fq.to_string_lossy().to_string()],
            &prefix,
        );
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files(&prefix);
    }

    // Subtest D: Control normalization without control file
    {
        let prefix = dir.path().join("err_ctrl").to_string_lossy().to_string();
        let mut args = default_count_args(
            &lib_path,
            vec![fq_path.to_string_lossy().to_string()],
            &prefix,
        );
        args.norm_method = NormMethod::Control;
        args.control_sgrna = None;
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files(&prefix);
    }

    // Subtest E: Nonexistent control file
    {
        let prefix = dir
            .path()
            .join("err_ctrl_file")
            .to_string_lossy()
            .to_string();
        let mut args = default_count_args(
            &lib_path,
            vec![fq_path.to_string_lossy().to_string()],
            &prefix,
        );
        args.norm_method = NormMethod::Control;
        args.control_sgrna = Some(dir.path().join("nonexistent_ctrl.txt"));
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files(&prefix);
    }

    // Subtest F: Invalid sgRNA length (0)
    {
        let prefix = dir.path().join("err_len_0").to_string_lossy().to_string();
        let mut args = default_count_args(
            &lib_path,
            vec![fq_path.to_string_lossy().to_string()],
            &prefix,
        );
        args.sgrna_len = 0;
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files(&prefix);
    }

    // Subtest G: Invalid sgRNA length (> 32)
    {
        let prefix = dir.path().join("err_len_33").to_string_lossy().to_string();
        let mut args = default_count_args(
            &lib_path,
            vec![fq_path.to_string_lossy().to_string()],
            &prefix,
        );
        args.sgrna_len = 33;
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files(&prefix);
    }

    // Subtest H: Invalid test_mismatch (> 1)
    {
        let prefix = dir.path().join("err_mm").to_string_lossy().to_string();
        let mut args = default_count_args(
            &lib_path,
            vec![fq_path.to_string_lossy().to_string()],
            &prefix,
        );
        args.test_mismatch = 2;
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files(&prefix);
    }

    // Subtest I: Duplicate sample labels
    {
        let prefix = dir
            .path()
            .join("err_dup_label")
            .to_string_lossy()
            .to_string();
        let mut args = default_count_args(
            &lib_path,
            vec![
                fq_path.to_string_lossy().to_string(),
                fq_path.to_string_lossy().to_string(),
            ],
            &prefix,
        );
        args.sample_label = Some(vec!["dup".to_string(), "dup".to_string()]);
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files(&prefix);
    }

    // Subtest J: Label count does not match group count
    {
        let prefix = dir
            .path()
            .join("err_label_count")
            .to_string_lossy()
            .to_string();
        let mut args = default_count_args(
            &lib_path,
            vec![fq_path.to_string_lossy().to_string()],
            &prefix,
        );
        args.sample_label = Some(vec!["s1".to_string(), "s2".to_string()]);
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files(&prefix);
    }

    // Subtest K: Empty path in comma-separated replicates
    {
        let prefix = dir.path().join("err_comma").to_string_lossy().to_string();
        let args = default_count_args(&lib_path, vec![format!("{},", fq_path.display())], &prefix);
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files(&prefix);
    }

    // Subtest L: Truncated/corrupt FASTQ record
    {
        let corrupt_fq = dir.path().join("corrupt.fastq");
        create_file(&corrupt_fq, "@r1\nACGT\n+\n"); // incomplete 4th line
        let prefix = dir.path().join("err_corrupt").to_string_lossy().to_string();
        let args = default_count_args(
            &lib_path,
            vec![corrupt_fq.to_string_lossy().to_string()],
            &prefix,
        );
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Fastq(_))));
        assert_no_output_files(&prefix);
    }

    // Subtest M: Median normalization failure when all counts are zero
    {
        let empty_fq = dir.path().join("all_zeros.fastq");
        create_file(&empty_fq, "");
        let prefix = dir
            .path()
            .join("err_zero_median")
            .to_string_lossy()
            .to_string();
        let mut args = default_count_args(
            &lib_path,
            vec![empty_fq.to_string_lossy().to_string()],
            &prefix,
        );
        args.norm_method = NormMethod::Median;
        let res = run_with_cli(Cli {
            command: Commands::Count(args),
        });
        assert!(matches!(res, Err(MageckError::Config(_))));
        assert_no_output_files(&prefix);
    }
}

// -----------------------------------------------------------------------------
// 16. CLI end-to-end parse and dispatch integration
// -----------------------------------------------------------------------------
#[test]
fn test_count_cli_parse_from_end_to_end() {
    let dir = TempDir::new().expect("tempdir");
    let lib_path = dir.path().join("library.tsv");
    create_standard_library(&lib_path);

    let fq = dir.path().join("sample.fastq");
    create_file(
        &fq,
        "@r1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r2\nCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    let prefix = dir.path().join("cli_out").to_string_lossy().to_string();

    let cli = Cli::try_parse_from([
        "mageckrs",
        "count",
        "-l",
        lib_path.to_str().expect("lib path str"),
        "-f",
        fq.to_str().expect("fq path str"),
        "-n",
        &prefix,
        "--sample-label",
        "cli_sample",
        "--norm-method",
        "none",
    ])
    .expect("Cli::try_parse_from");

    run_with_cli(cli).expect("run_with_cli from parsed CLI args");

    let count_file = format!("{prefix}.count.txt");
    let norm_file = format!("{prefix}.count.normalized.txt");
    let summary_file = format!("{prefix}.countsummary.txt");

    assert!(Path::new(&count_file).exists());
    assert!(Path::new(&norm_file).exists());
    assert!(Path::new(&summary_file).exists());

    let count_content = fs::read_to_string(&count_file).expect("read count");
    assert_eq!(
        count_content,
        "sgRNA\tGene\tcli_sample\ng1\tGENE_A\t1\ng2\tGENE_B\t1\n"
    );
}
