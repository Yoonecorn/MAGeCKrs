use std::fs;
use std::path::Path;
use tempfile::TempDir;

use mageckrs::count::matcher::SgRnaRecord;
use mageckrs::count::reader::{PhysicalFastqSummary, ReadSummary};
use mageckrs::count::writer::{
    write_count_outputs, write_count_summary, write_normalized_counts, write_raw_counts,
};
use mageckrs::data::count_matrix::CountMatrix;
use mageckrs::error::MageckError;

fn make_test_records(items: &[(&str, &str, &str)]) -> Vec<SgRnaRecord> {
    items
        .iter()
        .map(|(id, gene, seq)| SgRnaRecord {
            id: id.to_string(),
            gene: gene.to_string(),
            sequence: seq.to_string(),
        })
        .collect()
}

fn read_file(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()))
}

// 1. Exact raw header and row order
#[test]
fn test_1_exact_raw_header_and_row_order() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.count.txt");

    let records = make_test_records(&[
        ("g0", "GENE_A", "ACGTACGTACGTACGTACGT"),
        ("g1", "GENE_B", "TGCATGCATGCATGCATGCA"),
        ("g2", "GENE_C", "GGGGCCCCAAAATTTTGGGG"),
    ]);
    let labels = vec!["sample1"];
    let matrix = CountMatrix::from_vec(vec![10.0, 20.0, 30.0], 1, 3).unwrap();

    write_raw_counts(&out_file, &matrix, &records, &labels).expect("write raw counts");

    let content = read_file(&out_file);
    let lines: Vec<&str> = content.split('\n').collect();

    assert_eq!(lines[0], "sgRNA\tGene\tsample1");
    assert_eq!(lines[1], "g0\tGENE_A\t10");
    assert_eq!(lines[2], "g1\tGENE_B\t20");
    assert_eq!(lines[3], "g2\tGENE_C\t30");
    assert_eq!(lines[4], ""); // Trailing empty string from final newline
    assert_eq!(lines.len(), 5);
}

// 2. Multiple samples preserve label and matrix order
#[test]
fn test_2_multiple_samples_preserve_label_and_matrix_order() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.count.txt");

    let records = make_test_records(&[
        ("g0", "GENE_A", "ACGTACGTACGTACGTACGT"),
        ("g1", "GENE_B", "TGCATGCATGCATGCATGCA"),
    ]);
    let labels = vec!["sample_B", "sample_A", "sample_C"];
    // 3 samples x 2 guides:
    // sample 0 (sample_B): g0=1, g1=2
    // sample 1 (sample_A): g0=3, g1=4
    // sample 2 (sample_C): g0=5, g1=6
    let matrix = CountMatrix::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 3, 2).unwrap();

    write_raw_counts(&out_file, &matrix, &records, &labels).expect("write raw counts");

    let content = read_file(&out_file);
    let lines: Vec<&str> = content.split('\n').collect();

    assert_eq!(lines[0], "sgRNA\tGene\tsample_B\tsample_A\tsample_C");
    assert_eq!(lines[1], "g0\tGENE_A\t1\t3\t5");
    assert_eq!(lines[2], "g1\tGENE_B\t2\t4\t6");
}

// 3. Raw values serialize as unpadded integers
#[test]
fn test_3_raw_values_serialize_as_unpadded_integers() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.count.txt");

    let records = make_test_records(&[
        ("g0", "GENE_A", "ACGTACGTACGTACGTACGT"),
        ("g1", "GENE_B", "TGCATGCATGCATGCATGCA"),
        ("g2", "GENE_C", "GGGGCCCCAAAATTTTGGGG"),
        ("g3", "GENE_D", "AAAACCCCGGGGTTTTCCCC"),
    ]);
    let labels = vec!["s1"];
    let matrix = CountMatrix::from_vec(vec![0.0, 1.0, 42.0, 100000.0], 1, 4).unwrap();

    write_raw_counts(&out_file, &matrix, &records, &labels).unwrap();

    let content = read_file(&out_file);
    let lines: Vec<&str> = content.lines().collect();

    assert_eq!(lines[1], "g0\tGENE_A\t0");
    assert_eq!(lines[2], "g1\tGENE_B\t1");
    assert_eq!(lines[3], "g2\tGENE_C\t42");
    assert_eq!(lines[4], "g3\tGENE_D\t100000");
}

// 4. Raw writer rejects fractional values
#[test]
fn test_4_raw_writer_rejects_fractional_values() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.count.txt");

    let records = make_test_records(&[("g0", "GENE_A", "ACGTACGTACGTACGTACGT")]);
    let labels = vec!["s1"];

    let matrix = CountMatrix::from_vec(vec![1.5], 1, 1).unwrap();
    let res = write_raw_counts(&out_file, &matrix, &records, &labels);
    assert!(matches!(res, Err(MageckError::Config(_))));
    assert!(!out_file.exists());

    let matrix2 = CountMatrix::from_vec(vec![0.0001], 1, 1).unwrap();
    let res2 = write_raw_counts(&out_file, &matrix2, &records, &labels);
    assert!(matches!(res2, Err(MageckError::Config(_))));
    assert!(!out_file.exists());
}

// 5. Raw writer rejects negative values
#[test]
fn test_5_raw_writer_rejects_negative_values() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.count.txt");

    let records = make_test_records(&[("g0", "GENE_A", "ACGTACGTACGTACGTACGT")]);
    let labels = vec!["s1"];

    let matrix = CountMatrix::from_vec(vec![-1.0], 1, 1).unwrap();
    let res = write_raw_counts(&out_file, &matrix, &records, &labels);
    assert!(matches!(res, Err(MageckError::Config(_))));
    assert!(!out_file.exists());

    let matrix2 = CountMatrix::from_vec(vec![-5.5], 1, 1).unwrap();
    let res2 = write_raw_counts(&out_file, &matrix2, &records, &labels);
    assert!(matches!(res2, Err(MageckError::Config(_))));
    assert!(!out_file.exists());
}

// 6. Raw writer rejects NaN and infinity
#[test]
fn test_6_raw_writer_rejects_nan_and_infinity() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.count.txt");

    let records = make_test_records(&[("g0", "GENE_A", "ACGTACGTACGTACGTACGT")]);
    let labels = vec!["s1"];

    let matrix_nan = CountMatrix::from_vec(vec![f64::NAN], 1, 1).unwrap();
    let res_nan = write_raw_counts(&out_file, &matrix_nan, &records, &labels);
    assert!(matches!(res_nan, Err(MageckError::Config(_))));
    assert!(!out_file.exists());

    let matrix_inf = CountMatrix::from_vec(vec![f64::INFINITY], 1, 1).unwrap();
    let res_inf = write_raw_counts(&out_file, &matrix_inf, &records, &labels);
    assert!(matches!(res_inf, Err(MageckError::Config(_))));
    assert!(!out_file.exists());

    let matrix_neginf = CountMatrix::from_vec(vec![f64::NEG_INFINITY], 1, 1).unwrap();
    let res_neginf = write_raw_counts(&out_file, &matrix_neginf, &records, &labels);
    assert!(matches!(res_neginf, Err(MageckError::Config(_))));
    assert!(!out_file.exists());
}

// 7. Normalized values contain exactly four digits after the decimal point
#[test]
fn test_7_normalized_values_contain_exactly_four_digits() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.count.normalized.txt");

    let records = make_test_records(&[
        ("g0", "GENE_A", "ACGTACGTACGTACGTACGT"),
        ("g1", "GENE_B", "TGCATGCATGCATGCATGCA"),
    ]);
    let labels = vec!["s1"];
    let matrix = CountMatrix::from_vec(vec![1.234567, 42.1], 1, 2).unwrap();

    write_normalized_counts(&out_file, &matrix, &records, &labels).unwrap();

    let content = read_file(&out_file);
    let lines: Vec<&str> = content.lines().collect();

    assert_eq!(lines[1], "g0\tGENE_A\t1.2346");
    assert_eq!(lines[2], "g1\tGENE_B\t42.1000");
}

// 8. Normalized output preserves trailing zeroes
#[test]
fn test_8_normalized_output_preserves_trailing_zeroes() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.count.normalized.txt");

    let records = make_test_records(&[
        ("g0", "GENE_A", "ACGTACGTACGTACGTACGT"),
        ("g1", "GENE_B", "TGCATGCATGCATGCATGCA"),
        ("g2", "GENE_C", "GGGGCCCCAAAATTTTGGGG"),
    ]);
    let labels = vec!["s1"];
    let matrix = CountMatrix::from_vec(vec![0.0, 1.0, 2.5], 1, 3).unwrap();

    write_normalized_counts(&out_file, &matrix, &records, &labels).unwrap();

    let content = read_file(&out_file);
    let lines: Vec<&str> = content.lines().collect();

    assert_eq!(lines[1], "g0\tGENE_A\t0.0000");
    assert_eq!(lines[2], "g1\tGENE_B\t1.0000");
    assert_eq!(lines[3], "g2\tGENE_C\t2.5000");
}

// 9. Normalized serialization does not mutate source values
#[test]
fn test_9_normalized_serialization_does_not_mutate_source_values() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.count.normalized.txt");

    let records = make_test_records(&[
        ("g0", "GENE_A", "ACGTACGTACGTACGTACGT"),
        ("g1", "GENE_B", "TGCATGCATGCATGCATGCA"),
    ]);
    let labels = vec!["s1"];
    let initial_data = vec![1.234567890123, 9.876543210987];
    let matrix = CountMatrix::from_vec(initial_data.clone(), 1, 2).unwrap();

    write_normalized_counts(&out_file, &matrix, &records, &labels).unwrap();

    // Verify matrix in memory retains full precision and was not truncated or modified
    assert_eq!(matrix.as_slice(), initial_data.as_slice());
    assert_eq!(matrix.entry(0, 0), 1.234567890123);
    assert_eq!(matrix.entry(0, 1), 9.876543210987);
}

// 10. Normalized writer rejects NaN, infinity, and negative values
#[test]
fn test_10_normalized_writer_rejects_nan_infinity_and_negative_values() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.count.normalized.txt");

    let records = make_test_records(&[("g0", "GENE_A", "ACGTACGTACGTACGTACGT")]);
    let labels = vec!["s1"];

    let mat_nan = CountMatrix::from_vec(vec![f64::NAN], 1, 1).unwrap();
    assert!(matches!(
        write_normalized_counts(&out_file, &mat_nan, &records, &labels),
        Err(MageckError::Config(_))
    ));
    assert!(!out_file.exists());

    let mat_inf = CountMatrix::from_vec(vec![f64::INFINITY], 1, 1).unwrap();
    assert!(matches!(
        write_normalized_counts(&out_file, &mat_inf, &records, &labels),
        Err(MageckError::Config(_))
    ));
    assert!(!out_file.exists());

    let mat_neg = CountMatrix::from_vec(vec![-0.001], 1, 1).unwrap();
    assert!(matches!(
        write_normalized_counts(&out_file, &mat_neg, &records, &labels),
        Err(MageckError::Config(_))
    ));
    assert!(!out_file.exists());
}

// 11. Exact summary header and field order
#[test]
fn test_11_exact_summary_header_and_field_order() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.countsummary.txt");

    let raw_matrix = CountMatrix::from_vec(vec![10.0], 1, 1).unwrap();
    let fq_paths = vec!["sample1.fastq.gz"];
    let labels = vec!["s1"];
    let summaries = vec![ReadSummary {
        total_reads: 100,
        mapped_reads: 80,
        exact_matches: 70,
        one_mismatch_matches: 10,
        ambiguous_reads: 5,
        unmapped_reads: 10,
        short_reads: 5,
    }];

    let _records = make_test_records(&[("g", "G", "S")]);
    {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&summaries;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_summary(&out_file, &phys_summaries)
    }
    .unwrap();

    let content = read_file(&out_file);
    let lines: Vec<&str> = content.lines().collect();

    assert_eq!(
        lines[0],
        "File\tLabel\tReads\tMapped\tPercentage\tZerocounts"
    );
    assert_eq!(lines[1], "sample1.fastq.gz\ts1\t100\t80\t0.8000\t0");
}

// 12. Mapped equals exact plus one-mismatch counts
#[test]
fn test_12_mapped_equals_exact_plus_one_mismatch_counts() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.countsummary.txt");

    let raw_matrix = CountMatrix::from_vec(vec![1.0], 1, 1).unwrap();
    let fq_paths = vec!["data.fq"];
    let labels = vec!["sample_1"];

    // Valid invariant: mapped = 100 + 25 = 125
    let valid_summary = vec![ReadSummary {
        total_reads: 200,
        mapped_reads: 125,
        exact_matches: 100,
        one_mismatch_matches: 25,
        ambiguous_reads: 10,
        unmapped_reads: 50,
        short_reads: 15,
    }];
    let _records = make_test_records(&[("g", "G", "S")]);
    {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&valid_summary;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_summary(&out_file, &phys_summaries)
    }
    .unwrap();
    let content = read_file(&out_file);
    assert!(content.contains("data.fq\tsample_1\t200\t125\t0.6250\t0"));

    // Invalid invariant: mapped (120) != exact (100) + one_mismatch (25)
    let invalid_summary = vec![ReadSummary {
        total_reads: 200,
        mapped_reads: 120,
        exact_matches: 100,
        one_mismatch_matches: 25,
        ambiguous_reads: 10,
        unmapped_reads: 55,
        short_reads: 15,
    }];
    let res = {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&invalid_summary;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_summary(&out_file, &phys_summaries)
    };
    assert!(matches!(res, Err(MageckError::Config(_))));
}

// 13. Ambiguous, unmapped, and short reads do not enter Mapped
#[test]
fn test_13_ambiguous_unmapped_and_short_reads_do_not_enter_mapped() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.countsummary.txt");

    let raw_matrix = CountMatrix::from_vec(vec![1.0], 1, 1).unwrap();
    let fq_paths = vec!["reads.fq"];
    let labels = vec!["s1"];

    let summary = vec![ReadSummary {
        total_reads: 100,
        mapped_reads: 50, // exact(40) + mismatch(10)
        exact_matches: 40,
        one_mismatch_matches: 10,
        ambiguous_reads: 20,
        unmapped_reads: 20,
        short_reads: 10,
    }];

    let _records = make_test_records(&[("g", "G", "S")]);
    {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&summary;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_summary(&out_file, &phys_summaries)
    }
    .unwrap();

    let content = read_file(&out_file);
    let line = content.lines().nth(1).unwrap();
    let fields: Vec<&str> = line.split('\t').collect();

    assert_eq!(fields[2], "100"); // Reads
    assert_eq!(fields[3], "50"); // Mapped = 50, not 100
    assert_eq!(fields[4], "0.5000"); // Percentage = 0.5000 (mapped fraction)
}

// 14. Zero-read percentage is finite and zero
#[test]
fn test_14_zero_read_percentage_is_finite_and_zero() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.countsummary.txt");

    let raw_matrix = CountMatrix::from_vec(vec![0.0], 1, 1).unwrap();
    let fq_paths = vec!["empty.fastq"];
    let labels = vec!["empty_sample"];
    let summary = vec![ReadSummary {
        total_reads: 0,
        mapped_reads: 0,
        exact_matches: 0,
        one_mismatch_matches: 0,
        ambiguous_reads: 0,
        unmapped_reads: 0,
        short_reads: 0,
    }];

    let _records = make_test_records(&[("g", "G", "S")]);
    {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&summary;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_summary(&out_file, &phys_summaries)
    }
    .unwrap();

    let content = read_file(&out_file);
    let line = content.lines().nth(1).unwrap();
    let fields: Vec<&str> = line.split('\t').collect();

    assert_eq!(fields[2], "0"); // Reads
    assert_eq!(fields[3], "0"); // Mapped
    assert_eq!(fields[4], "0.0000"); // Percentage
    let pct: f64 = fields[4].parse().expect("parse percentage float");
    assert!(pct.is_finite());
    assert_eq!(pct, 0.0);
}

// 14b. Percentage column serializes mapped fraction in [0, 1] with exactly 4 decimal places
#[test]
fn test_14b_percentage_serializes_mapped_fraction_scale() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.countsummary.txt");

    let raw_matrix = CountMatrix::new(6, 1);
    let fq_paths = vec![
        "zero_reads.fq",
        "zero_mapped.fq",
        "quarter_mapped.fq",
        "half_mapped.fq",
        "all_mapped.fq",
        "ref_mapped.fq",
    ];
    let labels = vec![
        "s_zero_reads",
        "s_zero_mapped",
        "s_quarter",
        "s_half",
        "s_all",
        "s_ref",
    ];
    let summaries = vec![
        // 1. zero reads gives 0.0000
        ReadSummary {
            total_reads: 0,
            mapped_reads: 0,
            ..Default::default()
        },
        // 2. zero mapped out of nonzero reads gives 0.0000
        ReadSummary {
            total_reads: 100,
            mapped_reads: 0,
            unmapped_reads: 100,
            ..Default::default()
        },
        // 3. one-quarter mapped gives 0.2500
        ReadSummary {
            total_reads: 4,
            mapped_reads: 1,
            exact_matches: 1,
            unmapped_reads: 3,
            ..Default::default()
        },
        // 4. one-half mapped gives 0.5000
        ReadSummary {
            total_reads: 100,
            mapped_reads: 50,
            exact_matches: 50,
            unmapped_reads: 50,
            ..Default::default()
        },
        // 5. all reads mapped gives 1.0000
        ReadSummary {
            total_reads: 100,
            mapped_reads: 100,
            exact_matches: 100,
            ..Default::default()
        },
        // 6. reference values: 12,052,478 / 36,736,318 serialize as 0.3281
        ReadSummary {
            total_reads: 36_736_318,
            mapped_reads: 12_052_478,
            exact_matches: 12_052_478,
            unmapped_reads: 36_736_318 - 12_052_478,
            ..Default::default()
        },
    ];

    let _records = make_test_records(&[("g", "G", "S")]);
    {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&summaries;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_summary(&out_file, &phys_summaries)
    }
    .unwrap();

    let content = read_file(&out_file);
    let lines: Vec<&str> = content.lines().collect();

    assert_eq!(
        lines[0],
        "File\tLabel\tReads\tMapped\tPercentage\tZerocounts"
    );
    assert_eq!(lines[1], "zero_reads.fq\ts_zero_reads\t0\t0\t0.0000\t1");
    assert_eq!(lines[2], "zero_mapped.fq\ts_zero_mapped\t100\t0\t0.0000\t1");
    assert_eq!(lines[3], "quarter_mapped.fq\ts_quarter\t4\t1\t0.2500\t1");
    assert_eq!(lines[4], "half_mapped.fq\ts_half\t100\t50\t0.5000\t1");
    assert_eq!(lines[5], "all_mapped.fq\ts_all\t100\t100\t1.0000\t1");
    assert_eq!(
        lines[6],
        "ref_mapped.fq\ts_ref\t36736318\t12052478\t0.3281\t1"
    );
}

// 15. Zerocounts is computed per sample from the raw matrix
#[test]
fn test_15_zerocounts_computed_per_sample_from_raw_matrix() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.countsummary.txt");

    // 3 samples x 4 guides:
    // s0: [0, 1, 0, 5] -> 2 zerocounts
    // s1: [1, 2, 3, 4] -> 0 zerocounts
    // s2: [0, 0, 0, 0] -> 4 zerocounts
    let raw_matrix = CountMatrix::from_vec(
        vec![
            0.0, 1.0, 0.0, 5.0, // s0
            1.0, 2.0, 3.0, 4.0, // s1
            0.0, 0.0, 0.0, 0.0, // s2
        ],
        3,
        4,
    )
    .unwrap();

    let fq_paths = vec!["s0.fq", "s1.fq", "s2.fq"];
    let labels = vec!["s0", "s1", "s2"];
    let summaries = vec![
        ReadSummary {
            total_reads: 10,
            mapped_reads: 6,
            exact_matches: 6,
            unmapped_reads: 4,
            ..Default::default()
        },
        ReadSummary {
            total_reads: 20,
            mapped_reads: 10,
            exact_matches: 10,
            unmapped_reads: 10,
            ..Default::default()
        },
        ReadSummary {
            total_reads: 5,
            mapped_reads: 0,
            exact_matches: 0,
            unmapped_reads: 5,
            ..Default::default()
        },
    ];

    let _records = make_test_records(&[
        ("g", "G", "S"),
        ("g", "G", "S"),
        ("g", "G", "S"),
        ("g", "G", "S"),
    ]);
    {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&summaries;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_summary(&out_file, &phys_summaries)
    }
    .unwrap();

    let content = read_file(&out_file);
    let lines: Vec<&str> = content.lines().collect();

    assert_eq!(lines[1].split('\t').nth(5).unwrap(), "2");
    assert_eq!(lines[2].split('\t').nth(5).unwrap(), "0");
    assert_eq!(lines[3].split('\t').nth(5).unwrap(), "4");
}

// 16. Empty library writes the header and no data rows if CountMatrix supports zero guides
#[test]
fn test_16_empty_library_writes_header_and_no_data_rows() {
    let temp_dir = TempDir::new().unwrap();
    let raw_file = temp_dir.path().join("out.count.txt");
    let norm_file = temp_dir.path().join("out.count.normalized.txt");

    let matrix = CountMatrix::new(2, 0); // 2 samples, 0 guides
    let records: Vec<SgRnaRecord> = Vec::new();
    let labels = vec!["s1", "s2"];

    write_raw_counts(&raw_file, &matrix, &records, &labels).expect("write empty raw");
    let raw_content = read_file(&raw_file);
    assert_eq!(raw_content, "sgRNA\tGene\ts1\ts2\n");

    write_normalized_counts(&norm_file, &matrix, &records, &labels).expect("write empty norm");
    let norm_content = read_file(&norm_file);
    assert_eq!(norm_content, "sgRNA\tGene\ts1\ts2\n");
}

// 17. Dimension mismatches return errors before output is written where practical
#[test]
fn test_17_dimension_mismatches_return_errors_before_output_is_written() {
    let temp_dir = TempDir::new().unwrap();
    let matrix = CountMatrix::from_vec(vec![1.0, 2.0, 3.0, 4.0], 2, 2).unwrap();
    let records = make_test_records(&[
        ("g0", "GENE_A", "ACGTACGTACGTACGTACGT"),
        ("g1", "GENE_B", "TGCATGCATGCATGCATGCA"),
    ]);
    let fq_paths = ["s1.fq", "s2.fq"];
    let labels = vec!["s1", "s2"];
    let summaries = [
        ReadSummary {
            total_reads: 10,
            mapped_reads: 5,
            exact_matches: 5,
            unmapped_reads: 5,
            ..Default::default()
        },
        ReadSummary {
            total_reads: 20,
            mapped_reads: 15,
            exact_matches: 15,
            unmapped_reads: 5,
            ..Default::default()
        },
    ];

    // Combined writer dimension mismatch: raw [2, 2] vs norm [2, 3]
    let norm_matrix = CountMatrix::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 2, 3).unwrap();
    let prefix = temp_dir.path().join("combined_test");
    let res_combined = {
        let mut phys_summaries = Vec::new();
        for i in 0..summaries.len() {
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(fq_paths[i]),
                label: labels[i].to_string(),
                summary: summaries[i].clone(),
                zerocounts: matrix.sample_slice(i).iter().filter(|&&c| c == 0.0).count(),
            });
        }
        write_count_outputs(
            &prefix,
            &matrix,
            &norm_matrix,
            &records,
            &labels,
            &phys_summaries,
        )
    };
    assert!(matches!(res_combined, Err(MageckError::Config(_))));
    assert!(!temp_dir.path().join("combined_test.count.txt").exists());
}

// 18. Paths and labels are emitted in supplied order
#[test]
fn test_18_paths_and_labels_are_emitted_in_supplied_order() {
    let temp_dir = TempDir::new().unwrap();
    let out_file = temp_dir.path().join("out.countsummary.txt");

    let matrix = CountMatrix::from_vec(vec![1.0, 2.0, 3.0], 3, 1).unwrap();
    let fq_paths = vec!["/path/z.fastq", "/path/a.fastq", "/path/m.fastq"];
    let labels = vec!["sample_Z", "sample_A", "sample_M"];
    let summaries = vec![
        ReadSummary {
            total_reads: 10,
            mapped_reads: 5,
            exact_matches: 5,
            unmapped_reads: 5,
            ..Default::default()
        },
        ReadSummary {
            total_reads: 20,
            mapped_reads: 10,
            exact_matches: 10,
            unmapped_reads: 10,
            ..Default::default()
        },
        ReadSummary {
            total_reads: 30,
            mapped_reads: 15,
            exact_matches: 15,
            unmapped_reads: 15,
            ..Default::default()
        },
    ];

    let _records = make_test_records(&[("g", "G", "S")]);
    {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&summaries;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &matrix.sample_slice(i).iter().filter(|&&c| c == 0.0).count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_summary(&out_file, &phys_summaries)
    }
    .unwrap();

    let content = read_file(&out_file);
    let lines: Vec<&str> = content.lines().collect();

    assert!(lines[1].starts_with("/path/z.fastq\tsample_Z\t"));
    assert!(lines[2].starts_with("/path/a.fastq\tsample_A\t"));
    assert!(lines[3].starts_with("/path/m.fastq\tsample_M\t"));
}

// 19. All output uses tabs rather than commas
#[test]
fn test_19_all_output_uses_tabs_rather_than_commas() {
    let temp_dir = TempDir::new().unwrap();
    let prefix = temp_dir.path().join("tab_test");

    let records = make_test_records(&[
        ("g0", "GENE_A", "ACGTACGTACGTACGTACGT"),
        ("g1", "GENE_B", "TGCATGCATGCATGCATGCA"),
    ]);
    let labels = vec!["s1", "s2"];
    let raw_matrix = CountMatrix::from_vec(vec![1.0, 2.0, 3.0, 4.0], 2, 2).unwrap();
    let norm_matrix = CountMatrix::from_vec(vec![1.5, 2.5, 3.5, 4.5], 2, 2).unwrap();
    let fq_paths = vec!["s1.fq", "s2.fq"];
    let summaries = vec![
        ReadSummary {
            total_reads: 10,
            mapped_reads: 3,
            exact_matches: 3,
            unmapped_reads: 7,
            ..Default::default()
        },
        ReadSummary {
            total_reads: 10,
            mapped_reads: 7,
            exact_matches: 7,
            unmapped_reads: 3,
            ..Default::default()
        },
    ];

    let paths = {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&summaries;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_outputs(
            &prefix,
            &raw_matrix,
            &norm_matrix,
            &records,
            &labels,
            &phys_summaries,
        )
    }
    .unwrap();

    for path in [
        &paths.count_path,
        &paths.normalized_path,
        &paths.summary_path,
    ] {
        let content = read_file(path);
        for line in content.lines() {
            assert!(line.contains('\t'), "Line missing tabs: {line}");
            assert!(!line.contains(','), "Line contains comma delimiter: {line}");
        }
    }
}

// 20. Files use newline characters and end with a final newline
#[test]
fn test_20_files_use_newline_characters_and_end_with_final_newline() {
    let temp_dir = TempDir::new().unwrap();
    let prefix = temp_dir.path().join("newline_test");

    let records = make_test_records(&[("g0", "GENE_A", "ACGTACGTACGTACGTACGT")]);
    let labels = vec!["s1"];
    let raw_matrix = CountMatrix::from_vec(vec![42.0], 1, 1).unwrap();
    let norm_matrix = CountMatrix::from_vec(vec![42.0], 1, 1).unwrap();
    let fq_paths = vec!["s1.fq"];
    let summaries = vec![ReadSummary {
        total_reads: 100,
        mapped_reads: 42,
        exact_matches: 42,
        unmapped_reads: 58,
        ..Default::default()
    }];

    let paths = {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&summaries;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_outputs(
            &prefix,
            &raw_matrix,
            &norm_matrix,
            &records,
            &labels,
            &phys_summaries,
        )
    }
    .unwrap();

    for path in [
        &paths.count_path,
        &paths.normalized_path,
        &paths.summary_path,
    ] {
        let raw_bytes = fs::read(path).unwrap();
        // File must not be empty
        assert!(!raw_bytes.is_empty());
        // Must end with Unix newline '\n' (0x0A)
        assert_eq!(*raw_bytes.last().unwrap(), b'\n');
        // Must not contain Windows carriage returns '\r' (0x0D)
        assert!(!raw_bytes.contains(&b'\r'));
    }
}

// 21. Writers reject invalid TSV control characters
#[test]
fn test_21_writers_reject_invalid_tsv_control_characters() {
    let temp_dir = TempDir::new().unwrap();
    let prefix = temp_dir.path().join("tsv_test");

    let mut records = make_test_records(&[("g0", "GENE_A", "ACGTACGTACGTACGTACGT")]);
    let labels = vec!["s1\tsneaky"];
    let raw_matrix = CountMatrix::from_vec(vec![42.0], 1, 1).unwrap();
    let norm_matrix = CountMatrix::from_vec(vec![42.0], 1, 1).unwrap();
    let fq_paths = vec!["s1.fq"];
    let summaries = vec![ReadSummary {
        total_reads: 100,
        mapped_reads: 42,
        exact_matches: 42,
        unmapped_reads: 58,
        ..Default::default()
    }];

    // Label with tab
    let res = {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&summaries;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_outputs(
            &prefix,
            &raw_matrix,
            &norm_matrix,
            &records,
            &labels,
            &phys_summaries,
        )
    };
    assert!(matches!(res, Err(MageckError::Config(_))));

    // Gene with newline
    records[0].gene = "GENE\nA".to_string();
    let valid_labels = vec!["s1"];
    let res = {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&summaries;
        let paths_ref = &&fq_paths;
        let labels_ref = &&valid_labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_outputs(
            &prefix,
            &raw_matrix,
            &norm_matrix,
            &records,
            &valid_labels,
            &phys_summaries,
        )
    };
    assert!(matches!(res, Err(MageckError::Config(_))));

    // FASTQ path with carriage return
    records[0].gene = "GENE_A".to_string();
    let bad_paths = vec!["s1\r.fq"];
    let res = {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&summaries;
        let paths_ref = &&bad_paths;
        let labels_ref = &&valid_labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_outputs(
            &prefix,
            &raw_matrix,
            &norm_matrix,
            &records,
            &valid_labels,
            &phys_summaries,
        )
    };
    assert!(matches!(res, Err(MageckError::Config(_))));
}

// 22. Partial output prevention on numerical or invariant errors
#[test]
fn test_22_partial_output_prevention_on_numerical_errors() {
    let temp_dir = TempDir::new().unwrap();
    let prefix = temp_dir.path().join("partial_test");

    let records = make_test_records(&[("g0", "GENE_A", "ACGTACGTACGTACGTACGT")]);
    let labels = vec!["s1"];
    let fq_paths = vec!["s1.fq"];

    // Valid inputs
    let raw_matrix = CountMatrix::from_vec(vec![42.0], 1, 1).unwrap();
    let valid_norm = CountMatrix::from_vec(vec![42.0], 1, 1).unwrap();
    let valid_summary = vec![ReadSummary {
        total_reads: 100,
        mapped_reads: 42,
        exact_matches: 42,
        unmapped_reads: 58,
        ..Default::default()
    }];

    // Test 1: Invalid normalized matrix (NaN)
    let invalid_norm = CountMatrix::from_vec(vec![f64::NAN], 1, 1).unwrap();
    let res1 = {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&valid_summary;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_outputs(
            &prefix,
            &raw_matrix,
            &invalid_norm,
            &records,
            &labels,
            &phys_summaries,
        )
    };
    assert!(matches!(res1, Err(MageckError::Config(_))));
    assert!(!temp_dir.path().join("partial_test.count.txt").exists());
    assert!(!temp_dir
        .path()
        .join("partial_test.count.normalized.txt")
        .exists());
    assert!(!temp_dir
        .path()
        .join("partial_test.countsummary.txt")
        .exists());

    // Test 2: Invalid summary (invariant violation)
    let invalid_summary = vec![ReadSummary {
        total_reads: 100,
        mapped_reads: 42,
        exact_matches: 40, // mismatch!
        unmapped_reads: 58,
        ..Default::default()
    }];
    let res2 = {
        let mut phys_summaries = Vec::new();
        let summaries_ref = &&invalid_summary;
        let paths_ref = &&fq_paths;
        let labels_ref = &&labels;
        let len = std::cmp::min(
            summaries_ref.len(),
            std::cmp::min(paths_ref.len(), labels_ref.len()),
        );
        for i in 0..len {
            let zerocounts = &raw_matrix
                .sample_slice(i)
                .iter()
                .filter(|&&c| c == 0.0)
                .count();
            phys_summaries.push(PhysicalFastqSummary {
                path: std::path::PathBuf::from(paths_ref[i]),
                label: labels_ref[i].to_string(),
                summary: summaries_ref[i].clone(),
                zerocounts: *zerocounts,
            });
        }
        write_count_outputs(
            &prefix,
            &raw_matrix,
            &valid_norm,
            &records,
            &labels,
            &phys_summaries,
        )
    };
    assert!(matches!(res2, Err(MageckError::Config(_))));
    assert!(!temp_dir.path().join("partial_test.count.txt").exists());
    assert!(!temp_dir
        .path()
        .join("partial_test.count.normalized.txt")
        .exists());
    assert!(!temp_dir
        .path()
        .join("partial_test.countsummary.txt")
        .exists());
}
