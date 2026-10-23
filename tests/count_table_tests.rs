use std::io::{Cursor, Write};
use std::path::PathBuf;
use tempfile::NamedTempFile;

use mageckrs::data::{count_matrix::CountMatrix, ParsedCountTable};
use mageckrs::error::MageckError;
use mageckrs::mageck_test::{
    parse_count_table, parse_count_table_file, resolve_sample_indices, CountTable,
};

const VALID_TSV_4COL: &str = "\
sgRNA\tGene\tctrl\ttreat\n\
g1\tGENE_A\t10\t25\n\
g2\tGENE_B\t0\t50\n\
g3\tGENE_A\t100\t15\n";

#[test]
fn test_01_valid_count_table_label_based_selection() {
    let table = parse_count_table(Cursor::new(VALID_TSV_4COL), &["treat"], &["ctrl"])
        .expect("valid table must parse");

    assert_eq!(table.treatment_indices(), &[1]);
    assert_eq!(table.control_indices(), &[0]);
    assert_eq!(table.sample_labels(), &["ctrl", "treat"]);
    assert_eq!(table.n_guides(), 3);
    assert_eq!(table.n_samples(), 2);
}

#[test]
fn test_02_valid_count_table_index_based_selection() {
    let table = parse_count_table(Cursor::new(VALID_TSV_4COL), &["1"], &["0"])
        .expect("valid table with numeric indices must parse");

    assert_eq!(table.treatment_indices(), &[1]);
    assert_eq!(table.control_indices(), &[0]);
}

#[test]
fn test_03_exact_label_precedence_for_numeric_sample_label() {
    // Column index 0 has label "ctrl", column index 1 has literal label "0", column index 2 has label "treat"
    let tsv = "\
sgRNA\tGene\tctrl\t0\ttreat\n\
g1\tGENE_A\t10\t20\t30\n\
g2\tGENE_B\t40\t50\t60\n";

    // "0" matches sample label "0" at column index 1, NOT 0-based column index 0 ("ctrl")
    let table = parse_count_table(Cursor::new(tsv), &["0"], &["ctrl"])
        .expect("numeric label precedence must succeed");

    assert_eq!(table.treatment_indices(), &[1]);
    assert_eq!(table.control_indices(), &[0]);
}

#[test]
fn test_04_mixed_label_and_index_selection() {
    let tsv = "\
sgRNA\tGene\tc1\tc2\tt1\tt2\n\
g1\tGENE_A\t1\t2\t3\t4\n";

    // Treatment: label "t1" (idx 2) + index "3" (idx 3)
    // Control: index "0" (idx 0) + label "c2" (idx 1)
    let table = parse_count_table(Cursor::new(tsv), &["t1", "3"], &["0", "c2"])
        .expect("mixed label and index selection must succeed");

    assert_eq!(table.treatment_indices(), &[2, 3]);
    assert_eq!(table.control_indices(), &[0, 1]);
}

#[test]
fn test_05_preservation_of_sample_selection_order() {
    let tsv = "\
sgRNA\tGene\ts0\ts1\ts2\ts3\n\
g1\tGENE_A\t1\t2\t3\t4\n";

    let table = parse_count_table(Cursor::new(tsv), &["s3", "s1"], &["s2", "s0"])
        .expect("must preserve selection order");

    assert_eq!(table.treatment_indices(), &[3, 1]);
    assert_eq!(table.control_indices(), &[2, 0]);
}

#[test]
fn test_06_preservation_of_guide_row_order() {
    let tsv = "\
sgRNA\tGene\tctrl\ttreat\n\
guide_z\tGENE_Z\t1\t2\n\
guide_a\tGENE_A\t3\t4\n\
guide_m\tGENE_M\t5\t6\n\
guide_b\tGENE_B\t7\t8\n";

    let table = parse_count_table(Cursor::new(tsv), &["treat"], &["ctrl"])
        .expect("must preserve guide row order");

    assert_eq!(
        table.sgrna_ids(),
        &["guide_z", "guide_a", "guide_m", "guide_b"]
    );
    assert_eq!(table.genes(), &["GENE_Z", "GENE_A", "GENE_M", "GENE_B"]);
}

#[test]
fn test_07_correct_countmatrix_dimensions_and_values() {
    let tsv = "\
sgRNA\tGene\ts1\ts2\n\
g1\tG1\t10.5\t20.0\n\
g2\tG2\t30.0\t40.5\n\
g3\tG3\t50.0\t60.0\n";

    let table = parse_count_table(Cursor::new(tsv), &["s2"], &["s1"])
        .expect("must parse dimensions and fractional counts correctly");

    let matrix = table.counts();
    assert_eq!(matrix.n_samples(), 2);
    assert_eq!(matrix.n_guides(), 3);

    // Sample 0 (s1)
    assert_eq!(matrix.sample_slice(0), &[10.5, 30.0, 50.0]);
    // Sample 1 (s2)
    assert_eq!(matrix.sample_slice(1), &[20.0, 40.5, 60.0]);

    assert_eq!(matrix.entry(0, 0), 10.5);
    assert_eq!(matrix.entry(0, 1), 30.0);
    assert_eq!(matrix.entry(0, 2), 50.0);
    assert_eq!(matrix.entry(1, 0), 20.0);
    assert_eq!(matrix.entry(1, 1), 40.5);
    assert_eq!(matrix.entry(1, 2), 60.0);
}

#[test]
fn test_08_treatment_token_not_found() {
    let err = parse_count_table(Cursor::new(VALID_TSV_4COL), &["nonexistent"], &["ctrl"])
        .expect_err("unknown treatment token must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("nonexistent"),
                "error message must mention token: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_09_control_token_not_found() {
    let err = parse_count_table(Cursor::new(VALID_TSV_4COL), &["treat"], &["nonexistent"])
        .expect_err("unknown control token must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("nonexistent"),
                "error message must mention token: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_10_out_of_range_index() {
    // 2 samples: valid indices 0 and 1
    let err = parse_count_table(Cursor::new(VALID_TSV_4COL), &["2"], &["0"])
        .expect_err("out of range index 2 must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("out of range"),
                "error must mention out of range: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_11_duplicate_treatment_selection() {
    // Label and label
    let err = parse_count_table(Cursor::new(VALID_TSV_4COL), &["treat", "treat"], &["ctrl"])
        .expect_err("duplicate treatment selection must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("Duplicate treatment"),
                "error must identify duplicate: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }

    // Label and index resolving to same column
    let err2 = parse_count_table(Cursor::new(VALID_TSV_4COL), &["treat", "1"], &["ctrl"])
        .expect_err("label and index resolving to same treatment must fail");

    match err2 {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("Duplicate treatment"),
                "error must identify duplicate: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_12_duplicate_control_selection() {
    let err = parse_count_table(Cursor::new(VALID_TSV_4COL), &["treat"], &["ctrl", "ctrl"])
        .expect_err("duplicate control selection must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("Duplicate control"),
                "error must identify duplicate: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_13_overlap_between_treatment_and_control_selections() {
    // Same label in both
    let err = parse_count_table(Cursor::new(VALID_TSV_4COL), &["ctrl"], &["ctrl"])
        .expect_err("same column in treatment and control must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("both treatment and control"),
                "error must identify overlap: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }

    // Label and index resolving to same column
    let err2 = parse_count_table(Cursor::new(VALID_TSV_4COL), &["0"], &["ctrl"])
        .expect_err("index 0 and label ctrl overlap must fail");

    match err2 {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("both treatment and control"),
                "error must identify overlap: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_14_duplicate_sample_header() {
    let tsv = "\
sgRNA\tGene\tsample1\tsample1\n\
g1\tGENE_A\t10\t20\n";

    let err = ParsedCountTable::from_reader(Cursor::new(tsv))
        .expect_err("duplicate sample header must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("duplicate") || msg.contains("Duplicate"),
                "must identify duplicate sample header: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_15_duplicate_sgrna_identifier() {
    let tsv = "\
sgRNA\tGene\tctrl\ttreat\n\
g1\tGENE_A\t10\t20\n\
g2\tGENE_B\t30\t40\n\
g1\tGENE_C\t50\t60\n";

    let err =
        ParsedCountTable::from_reader(Cursor::new(tsv)).expect_err("duplicate sgRNA ID must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("duplicate sgRNA identifier 'g1'"),
                "must report duplicate guide: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_16_missing_sgrna_header() {
    // Missing sgRNA column header (has "guide" instead)
    let tsv = "\
guide\tGene\tctrl\ttreat\n\
g1\tGENE_A\t10\t20\n";

    let err = ParsedCountTable::from_reader(Cursor::new(tsv))
        .expect_err("missing sgRNA header must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("sgRNA"),
                "error must mention missing sgRNA: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }

    // Lowercase sgrna header rejected
    let tsv_lower = "\
sgrna\tGene\tctrl\ttreat\n\
g1\tGENE_A\t10\t20\n";

    let err_lower = ParsedCountTable::from_reader(Cursor::new(tsv_lower))
        .expect_err("lowercase sgrna header must fail");

    match err_lower {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("sgRNA"),
                "error must mention expected sgRNA: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_17_missing_gene_header() {
    let tsv = "\
sgRNA\tTarget\tctrl\ttreat\n\
g1\tGENE_A\t10\t20\n";

    let err =
        ParsedCountTable::from_reader(Cursor::new(tsv)).expect_err("missing Gene header must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("Gene"),
                "error must mention missing Gene: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }

    // Lowercase gene header rejected
    let tsv_lower = "\
sgRNA\tgene\tctrl\ttreat\n\
g1\tGENE_A\t10\t20\n";

    let err_lower = ParsedCountTable::from_reader(Cursor::new(tsv_lower))
        .expect_err("lowercase gene header must fail");

    match err_lower {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("Gene"),
                "error must mention expected Gene: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_18_no_sample_columns() {
    let tsv = "\
sgRNA\tGene\n\
g1\tGENE_A\n";

    let err =
        ParsedCountTable::from_reader(Cursor::new(tsv)).expect_err("no sample columns must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("sample column"),
                "error must mention lack of sample columns: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_19_empty_data_body() {
    let tsv = "sgRNA\tGene\tctrl\ttreat\n";

    let err =
        ParsedCountTable::from_reader(Cursor::new(tsv)).expect_err("empty data body must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("no data rows"),
                "error must mention no data rows: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_20_empty_sgrna_identifier() {
    let tsv = "\
sgRNA\tGene\tctrl\ttreat\n\
\tGENE_A\t10\t20\n";

    let err = ParsedCountTable::from_reader(Cursor::new(tsv))
        .expect_err("empty sgRNA identifier must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("sgRNA identifier is empty"),
                "error must identify empty sgRNA: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_21_empty_gene_identifier() {
    let tsv = "\
sgRNA\tGene\tctrl\ttreat\n\
g1\t\t10\t20\n";

    let err = ParsedCountTable::from_reader(Cursor::new(tsv))
        .expect_err("empty gene identifier must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("gene identifier is empty"),
                "error must identify empty gene: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_22_missing_field() {
    let tsv = "\
sgRNA\tGene\tctrl\ttreat\n\
g1\tGENE_A\t10\n";

    let err = ParsedCountTable::from_reader(Cursor::new(tsv))
        .expect_err("missing field on data row must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("field count mismatch"),
                "error must identify field mismatch: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_23_extra_field() {
    let tsv = "\
sgRNA\tGene\tctrl\ttreat\n\
g1\tGENE_A\t10\t20\t999\n";

    let err = ParsedCountTable::from_reader(Cursor::new(tsv))
        .expect_err("extra field on data row must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("field count mismatch"),
                "error must identify field mismatch: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_24_non_numeric_count() {
    let tsv = "\
sgRNA\tGene\tctrl\ttreat\n\
g1\tGENE_A\t10\tabc\n";

    let err = ParsedCountTable::from_reader(Cursor::new(tsv))
        .expect_err("non-numeric count value must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("non-numeric count 'abc'"),
                "error must identify non-numeric count: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_25_negative_count() {
    let tsv = "\
sgRNA\tGene\tctrl\ttreat\n\
g1\tGENE_A\t10\t-5\n";

    let err = ParsedCountTable::from_reader(Cursor::new(tsv))
        .expect_err("negative count value must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("negative count '-5'"),
                "error must identify negative count: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_26_nan_count() {
    let tsv = "\
sgRNA\tGene\tctrl\ttreat\n\
g1\tGENE_A\t10\tNaN\n";

    let err =
        ParsedCountTable::from_reader(Cursor::new(tsv)).expect_err("NaN count value must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("non-finite count 'NaN'"),
                "error must identify NaN count: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_27_positive_infinity_count() {
    let tsv = "\
sgRNA\tGene\tctrl\ttreat\n\
g1\tGENE_A\t10\tinf\n";

    let err = ParsedCountTable::from_reader(Cursor::new(tsv))
        .expect_err("positive infinity count value must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("non-finite count 'inf'"),
                "error must identify infinite count: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_28_negative_infinity_count() {
    let tsv = "\
sgRNA\tGene\tctrl\ttreat\n\
g1\tGENE_A\t-inf\t20\n";

    let err = parse_count_table(Cursor::new(tsv), &["treat"], &["ctrl"])
        .expect_err("negative infinity count value must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("non-finite count '-inf'"),
                "error must identify infinite count: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_29_crlf_input() {
    let crlf_tsv = "sgRNA\tGene\tctrl\ttreat\r\ng1\tGENE_A\t10\t25\r\ng2\tGENE_B\t0\t50\r\n";

    let table = ParsedCountTable::from_reader(Cursor::new(crlf_tsv))
        .expect("CRLF input must parse successfully");

    assert_eq!(table.n_guides(), 2);
    assert_eq!(table.sample_labels(), &["ctrl", "treat"]);
    assert_eq!(table.counts().sample_slice(1), &[25.0, 50.0]);
}

#[test]
fn test_30_mandatory_final_newline_and_missing_final_newline() {
    let with_newline = "sgRNA\tGene\tctrl\ttreat\ng1\tGENE_A\t10\t20\n";
    let without_newline = "sgRNA\tGene\tctrl\ttreat\ng1\tGENE_A\t10\t20";

    let t1 = ParsedCountTable::from_reader(Cursor::new(with_newline))
        .expect("with final newline must parse");
    let t2 = ParsedCountTable::from_reader(Cursor::new(without_newline))
        .expect("without final newline must parse");

    assert_eq!(t1, t2);
}

#[test]
fn test_31_selection_token_with_whitespace_rejected() {
    let err = parse_count_table(Cursor::new(VALID_TSV_4COL), &[" treat "], &["ctrl"])
        .expect_err("token with whitespace must be rejected");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("leading or trailing whitespace"),
                "error must mention whitespace: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_32_selection_token_case_sensitivity() {
    let tsv = "\
sgRNA\tGene\tCtrl\tTreat\n\
g1\tGENE_A\t10\t20\n";

    // Lowercase token should not match uppercase label
    let err = parse_count_table(Cursor::new(tsv), &["treat"], &["Ctrl"])
        .expect_err("case-sensitive token must not match differently cased label");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("treat"),
                "error must identify unmatched token: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_33_empty_selection_token_rejected() {
    let err = parse_count_table(Cursor::new(VALID_TSV_4COL), &[""], &["ctrl"])
        .expect_err("empty treatment token must fail");

    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("must not be empty"),
                "error must state token cannot be empty: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_34_parse_count_table_file_tempfile() {
    let mut tmp = NamedTempFile::new().expect("create temp file");
    tmp.write_all(VALID_TSV_4COL.as_bytes())
        .expect("write TSV to temp file");
    tmp.flush().expect("flush temp file");

    let table = parse_count_table_file(tmp.path(), &["treat"], &["ctrl"])
        .expect("parsing valid temp file must succeed");

    assert_eq!(table.n_guides(), 3);
    assert_eq!(table.n_samples(), 2);
    assert_eq!(table.treatment_indices(), &[1]);
    assert_eq!(table.control_indices(), &[0]);
}

#[test]
fn test_35_parse_count_table_file_not_found() {
    let missing_path = PathBuf::from("/non/existent/path/count_table.txt");
    let err = parse_count_table_file(&missing_path, &["treat"], &["ctrl"])
        .expect_err("missing file must return error");

    match err {
        MageckError::Io(io_err) => {
            assert_eq!(io_err.kind(), std::io::ErrorKind::NotFound);
            assert!(
                io_err.to_string().contains(missing_path.to_str().unwrap()),
                "error message must contain file path: {io_err}"
            );
        }
        other => panic!("expected MageckError::Io, got {other:?}"),
    }
}

#[test]
fn test_36_count_table_accessors_and_into_parts() {
    let table = parse_count_table(Cursor::new(VALID_TSV_4COL), &["treat"], &["ctrl"])
        .expect("must parse valid table");

    assert_eq!(table.sgrna_ids(), &["g1", "g2", "g3"]);
    assert_eq!(table.sgrnas(), &["g1", "g2", "g3"]);
    assert_eq!(table.genes(), &["GENE_A", "GENE_B", "GENE_A"]);
    assert_eq!(table.sample_labels(), &["ctrl", "treat"]);
    assert_eq!(table.treatment_indices(), &[1]);
    assert_eq!(table.control_indices(), &[0]);
    assert_eq!(table.n_guides(), 3);
    assert_eq!(table.n_samples(), 2);
    assert!(!table.is_empty());
    assert_eq!(table.sample_index("ctrl"), Some(0));
    assert_eq!(table.sample_index("treat"), Some(1));
    assert_eq!(table.sample_index("missing"), None);

    let (sgrnas, genes, labels, matrix, treat_idx, ctrl_idx) = table.clone().into_parts();
    assert_eq!(sgrnas, vec!["g1", "g2", "g3"]);
    assert_eq!(genes, vec!["GENE_A", "GENE_B", "GENE_A"]);
    assert_eq!(labels, vec!["ctrl", "treat"]);
    assert_eq!(treat_idx, vec![1]);
    assert_eq!(ctrl_idx, vec![0]);
    assert_eq!(matrix.n_guides(), 3);

    let counts = table.into_counts();
    assert_eq!(counts.n_samples(), 2);
    assert_eq!(counts.n_guides(), 3);
}

#[test]
fn test_37_resolve_sample_indices_direct() {
    let labels = vec![
        "ctrl_1".to_string(),
        "ctrl_2".to_string(),
        "trt_1".to_string(),
        "trt_2".to_string(),
    ];

    let (treat, ctrl) = resolve_sample_indices(&labels, &["trt_2", "2"], &["0", "ctrl_2"])
        .expect("resolution must succeed");

    assert_eq!(treat, vec![3, 2]);
    assert_eq!(ctrl, vec![0, 1]);
}

#[test]
fn test_38_count_table_constructors() {
    let matrix = CountMatrix::from_slice(&[1.0, 2.0, 3.0, 4.0], 2, 2).expect("valid matrix");
    let table = CountTable::new(
        vec!["g1".into(), "g2".into()],
        vec!["G1".into(), "G2".into()],
        vec!["ctrl".into(), "treat".into()],
        matrix.clone(),
        &["treat"],
        &["ctrl"],
    )
    .expect("CountTable::new must succeed");

    assert_eq!(table.n_guides(), 2);
    assert_eq!(table.n_samples(), 2);
    assert_eq!(table.treatment_indices(), &[1]);
    assert_eq!(table.control_indices(), &[0]);

    // from_parts with dimension mismatch
    let bad = CountTable::from_parts(
        vec!["g1".into()], // only 1 guide, but matrix has 2
        vec!["G1".into()],
        vec!["ctrl".into(), "treat".into()],
        matrix,
        vec![1],
        vec![0],
    );
    assert!(bad.is_err());
}

#[test]
fn test_39_parsed_count_table_direct_api() {
    let table = ParsedCountTable::from_reader(Cursor::new(VALID_TSV_4COL))
        .expect("ParsedCountTable::from_reader must parse valid TSV");

    assert_eq!(table.sgrna_ids(), &["g1", "g2", "g3"]);
    assert_eq!(table.genes(), &["GENE_A", "GENE_B", "GENE_A"]);
    assert_eq!(table.sample_labels(), &["ctrl", "treat"]);
    assert_eq!(table.n_guides(), 3);
    assert_eq!(table.n_samples(), 2);
    assert!(!table.is_empty());
    assert_eq!(table.sample_index("ctrl"), Some(0));
    assert_eq!(table.sample_index("treat"), Some(1));
    assert_eq!(table.sample_index("unknown"), None);

    let (sgrna_ids, genes, sample_labels, counts) = table.into_parts();
    assert_eq!(sgrna_ids.len(), 3);
    assert_eq!(genes.len(), 3);
    assert_eq!(sample_labels.len(), 2);
    assert_eq!(counts.n_guides(), 3);
    assert_eq!(counts.n_samples(), 2);
}
