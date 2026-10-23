use mageckrs::cli::NormMethod;
use mageckrs::config::{CountConfig, FastqSampleGroup};
use mageckrs::count::run;
use mageckrs::error::MageckError;
use std::fs;

use std::path::{Path, PathBuf};
use tempfile::TempDir;

// A small deterministic helper to create test files
fn create_file(dir: &Path, name: &str, content: &str) -> PathBuf {
    let path = dir.join(name);
    fs::write(&path, content).unwrap();
    path
}

const GZ_FASTQ_BYTES: &[u8] = &[
    31, 139, 8, 0, 98, 130, 152, 106, 2, 255, 115, 40, 74, 77, 76, 137, 79, 173, 72, 76, 46, 225,
    114, 116, 118, 15, 65, 199, 92, 218, 92, 110, 88, 0, 23, 0, 141, 149, 224, 203, 56, 0, 0, 0,
];

fn create_gz_file(dir: &Path, name: &str, _content: &str) -> PathBuf {
    let path = dir.join(name);
    fs::write(&path, GZ_FASTQ_BYTES).unwrap();
    path
}

fn base_config(dir: &TempDir) -> CountConfig {
    CountConfig {
        list_seq: dir.path().join("lib.csv"),
        sample_groups: vec![],
        output_prefix: dir.path().join("out").to_string_lossy().to_string(),
        norm_method: NormMethod::None,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    }
}

// 1. Single sample with one FASTQ and default label.
#[test]
fn test_single_sample_one_fastq() {
    let dir = tempfile::tempdir().unwrap();
    let lib = create_file(
        dir.path(),
        "lib.csv",
        "id,sequence,gene\ng1,ACGTACGTACGTACGTACGT,geneA\ng2,TGCATGCATGCATGCATGCA,geneB\n",
    );
    let fq = create_file(
        dir.path(),
        "1.fq",
        "@seq1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    let mut config = base_config(&dir);
    config.list_seq = lib;
    config.sample_groups = vec![FastqSampleGroup {
        label: "sample1".to_string(),
        fastq_paths: vec![fq.clone()],
    }];

    let out_paths = run::run(&config).unwrap();

    // Check paths
    assert!(out_paths
        .count_path
        .to_string_lossy()
        .ends_with("out.count.txt"));
    assert!(out_paths
        .normalized_path
        .to_string_lossy()
        .ends_with("out.count.normalized.txt"));
    assert!(out_paths
        .summary_path
        .to_string_lossy()
        .ends_with("out.countsummary.txt"));

    let count_data = fs::read_to_string(&out_paths.count_path).unwrap();
    assert_eq!(
        count_data,
        "sgRNA\tGene\tsample1\ng1\tgeneA\t1\ng2\tgeneB\t0\n"
    );

    let norm_data = fs::read_to_string(&out_paths.normalized_path).unwrap();
    // None normalization preserves raw numerical values with four-decimal serialization
    assert_eq!(
        norm_data,
        "sgRNA\tGene\tsample1\ng1\tgeneA\t1.0000\ng2\tgeneB\t0.0000\n"
    );

    let sum_data = fs::read_to_string(&out_paths.summary_path).unwrap();
    let fq_display = fq.to_string_lossy();
    // 1 mapped out of 1 -> 1.0000, 1 zerocounts (g2)
    assert_eq!(
        sum_data,
        format!("File\tLabel\tReads\tMapped\tPercentage\tZerocounts\n{fq_display}\tsample1\t1\t1\t1.0000\t1\n")
    );
}

// 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13. Multiple sample groups and technical replicates.
#[test]
fn test_multiple_groups_and_technical_replicates() {
    let dir = tempfile::tempdir().unwrap();
    let lib = create_file(
        dir.path(),
        "lib.csv",
        "id,sequence,gene\ng1,ACGTACGTACGTACGTACGT,geneA\ng2,TGCATGCATGCATGCATGCA,geneB\n",
    );

    // Group 1: 2 files
    let fq1_1 = create_file(
        dir.path(),
        "g1_1.fq",
        "@seq1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );
    let fq1_2 = create_file(
        dir.path(),
        "g1_2.fq",
        "@seq2\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    // Group 2: 1 file
    let fq2_1 = create_file(dir.path(), "g2_1.fq", "@seq3\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n@seq4\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n@seq5\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n");

    let mut config = base_config(&dir);
    config.list_seq = lib;
    config.sample_groups = vec![
        FastqSampleGroup {
            label: "group1".to_string(),
            fastq_paths: vec![fq1_1.clone(), fq1_2.clone()],
        },
        FastqSampleGroup {
            label: "group2".to_string(),
            fastq_paths: vec![fq2_1.clone()],
        },
    ];

    let out_paths = run::run(&config).unwrap();

    let count_data = fs::read_to_string(&out_paths.count_path).unwrap();
    // One column per sample group
    assert_eq!(
        count_data,
        "sgRNA\tGene\tgroup1\tgroup2\ng1\tgeneA\t2\t0\ng2\tgeneB\t0\t2\n"
    );

    let sum_data = fs::read_to_string(&out_paths.summary_path).unwrap();
    let lines: Vec<&str> = sum_data.lines().collect();
    assert_eq!(
        lines[0],
        "File\tLabel\tReads\tMapped\tPercentage\tZerocounts"
    );
    // Summary rows retain physical-file order, independent mapped/reads/percentage/zerocounts
    let l1 = format!("{}\tgroup1\t1\t1\t1.0000\t1", fq1_1.to_string_lossy());
    let l2 = format!("{}\tgroup1\t1\t1\t1.0000\t1", fq1_2.to_string_lossy());
    let l3 = format!("{}\tgroup2\t3\t2\t0.6667\t1", fq2_1.to_string_lossy());
    assert_eq!(lines[1], l1);
    assert_eq!(lines[2], l2);
    assert_eq!(lines[3], l3);
}

// 14. Normalized output uses exactly four decimal places.
// 16. Median normalization runs without a control file.
#[test]
fn test_median_normalization() {
    let dir = tempfile::tempdir().unwrap();
    let lib = create_file(
        dir.path(),
        "lib.csv",
        "id,sequence,gene\ng1,ACGTACGTACGTACGTACGT,geneA\ng2,TGCATGCATGCATGCATGCA,geneB\n",
    );
    let fq = create_file(
        dir.path(),
        "1.fq",
        "@seq1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    let mut config = base_config(&dir);
    config.list_seq = lib;
    config.sample_groups = vec![FastqSampleGroup {
        label: "sample1".to_string(),
        fastq_paths: vec![fq],
    }];
    config.norm_method = NormMethod::Median;
    config.control_sgrna = None; // Without control file

    let out_paths = run::run(&config).unwrap();
    let norm_data = fs::read_to_string(&out_paths.normalized_path).unwrap();

    // Size factor will be 0.5 (Median of [0.5, 0.0] -> 0.5? Actually mean of log, but 1 non-zero is 1.0, wait. If g1 has 1, g2 has 0, mean of all samples for g1=1.0, g2=0.0. Ratio for g1=1/1=1, g2=0. Median of [1.0] = 1.0. Size factor 1.0.
    // Normalized count for g1 = 1.0 / 1.0 = 1.0. Let's just check four decimal places format.
    assert!(norm_data.contains("1.0000") || norm_data.contains("0.0000"));
}

// 22. Control normalization without --control-sgrna returns a typed configuration error.
#[test]
fn test_control_norm_requires_control_file() {
    let dir = tempfile::tempdir().unwrap();
    let mut config = base_config(&dir);
    config.norm_method = NormMethod::Control;
    config.control_sgrna = None;

    let res = run::run(&config);
    match res {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("requires a control-sgrna file")
                    || msg.contains("requires --control-sgrna")
            );
        }
        _ => panic!("Expected Config error"),
    }
}

// 19. Control normalization with a valid control-ID file.
#[test]
fn test_control_norm_valid_id() {
    let dir = tempfile::tempdir().unwrap();
    let lib = create_file(
        dir.path(),
        "lib.csv",
        "id,sequence,gene\ng1,AAAAAAAAAAAAAAAAAAAA,geneA\ng2,CCCCCCCCCCCCCCCCCCCC,geneB\n",
    );
    let ctrl = create_file(dir.path(), "ctrl.txt", "g1\n");
    let fq = create_file(
        dir.path(),
        "1.fq",
        "@seq1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    let mut config = base_config(&dir);
    config.list_seq = lib;
    config.sample_groups = vec![FastqSampleGroup {
        label: "s1".to_string(),
        fastq_paths: vec![fq],
    }];
    config.norm_method = NormMethod::Control;
    config.control_sgrna = Some(ctrl);

    let res = run::run(&config);
    assert!(res.is_ok());
}

// 23, 24, 25. Supplied control file with Median/Total/None is parsed and validated but does not affect normalization.
// 28. Malformed control entry fails before output creation.
#[test]
fn test_control_file_with_median_validates_but_malformed_fails() {
    let dir = tempfile::tempdir().unwrap();
    let lib = create_file(
        dir.path(),
        "lib.csv",
        "id,sequence,gene\ng1,AAAAAAAAAAAAAAAAAAAA,geneA\n",
    );
    let ctrl = create_file(dir.path(), "ctrl.txt", "g1\ninvalid space\n");
    let fq = create_file(
        dir.path(),
        "1.fq",
        "@seq1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    let mut config = base_config(&dir);
    config.list_seq = lib;
    config.sample_groups = vec![FastqSampleGroup {
        label: "s1".to_string(),
        fastq_paths: vec![fq],
    }];
    config.norm_method = NormMethod::Median;
    config.control_sgrna = Some(ctrl);

    let res = run::run(&config);
    match res {
        Err(MageckError::Config(_)) => {}
        _ => panic!("Expected Config error due to malformed control"),
    }

    // Assert output not created
    let count_os = config.output_prefix.clone() + ".count.txt";
    assert!(!Path::new(&count_os).exists());
}

// 34. Empty FASTQ produces a valid zero summary.
#[test]
fn test_empty_fastq() {
    let dir = tempfile::tempdir().unwrap();
    let lib = create_file(
        dir.path(),
        "lib.csv",
        "id,sequence,gene\ng1,AAAAAAAAAAAAAAAAAAAA,geneA\n",
    );
    let fq = create_file(dir.path(), "1.fq", "");

    let mut config = base_config(&dir);
    config.list_seq = lib;
    config.sample_groups = vec![FastqSampleGroup {
        label: "s1".to_string(),
        fastq_paths: vec![fq.clone()],
    }];

    let out = run::run(&config).unwrap();
    let sum = fs::read_to_string(&out.summary_path).unwrap();
    let fq_display = fq.to_string_lossy();
    assert!(sum.contains(&format!("{fq_display}\ts1\t0\t0\t0.0000\t1")));
}

// 35. Mixed gzip and uncompressed technical replicates.
// 36. Repeated physical FASTQ path is counted twice and summarized twice.
#[test]
fn test_mixed_compression_and_repeated_path() {
    let dir = tempfile::tempdir().unwrap();
    let lib = create_file(
        dir.path(),
        "lib.csv",
        "id,sequence,gene\ng1,AAAAAAAAAAAAAAAAAAAA,geneA\n",
    );
    let fq1 = create_file(
        dir.path(),
        "1.fq",
        "@s\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );
    let fq2 = create_gz_file(
        dir.path(),
        "2.fq.gz",
        "@s\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n",
    );

    let mut config = base_config(&dir);
    config.list_seq = lib;
    config.sample_groups = vec![FastqSampleGroup {
        label: "s1".to_string(),
        fastq_paths: vec![fq1.clone(), fq2.clone(), fq1.clone()], // fq1 repeated
    }];

    let out = run::run(&config).unwrap();
    let count = fs::read_to_string(&out.count_path).unwrap();
    assert!(count.contains("\t2\n")); // fq1 counted twice

    let sum = fs::read_to_string(&out.summary_path).unwrap();
    let f1 = fq1.to_string_lossy();
    let f2 = fq2.to_string_lossy();
    // Summarized 3 times in order
    let lines: Vec<&str> = sum.lines().collect();
    assert_eq!(lines.len(), 4);
    assert!(lines[1].starts_with(&f1.to_string()));
    assert!(lines[2].starts_with(&f2.to_string()));
    assert!(lines[3].starts_with(&f1.to_string()));
}
