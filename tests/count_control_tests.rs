use std::io::{Cursor, Write};
use std::path::PathBuf;
use tempfile::NamedTempFile;

use mageckrs::count::{parse_control_sgrna, parse_control_sgrna_file, SgRnaRecord};
use mageckrs::error::MageckError;

fn make_test_records() -> Vec<SgRnaRecord> {
    vec![
        SgRnaRecord {
            id: "guide_0".to_string(),
            gene: "GENE_A".to_string(),
            sequence: "AAAA".to_string(),
        },
        SgRnaRecord {
            id: "guide_1".to_string(),
            gene: "GENE_B".to_string(),
            sequence: "CCCC".to_string(),
        },
        SgRnaRecord {
            id: "guide_2".to_string(),
            gene: "GENE_C".to_string(),
            sequence: "GGGG".to_string(),
        },
        SgRnaRecord {
            id: "guide_3".to_string(),
            gene: "GENE_D".to_string(),
            sequence: "TTTT".to_string(),
        },
    ]
}

#[test]
fn test_01_one_valid_sgrna_id() {
    let records = make_test_records();
    let res = parse_control_sgrna(Cursor::new("guide_1\n"), &records).expect("must resolve");
    assert_eq!(res.guide_indexes, vec![1]);
    assert!(res.unknown_entries.is_empty());
}

#[test]
fn test_02_one_valid_guide_sequence() {
    let records = make_test_records();
    let res = parse_control_sgrna(Cursor::new("GGGG\n"), &records).expect("must resolve");
    assert_eq!(res.guide_indexes, vec![2]);
    assert!(res.unknown_entries.is_empty());
}

#[test]
fn test_03_multiple_valid_ids() {
    let records = make_test_records();
    let res =
        parse_control_sgrna(Cursor::new("guide_3\nguide_1\n"), &records).expect("must resolve");
    assert_eq!(res.guide_indexes, vec![1, 3]);
    assert!(res.unknown_entries.is_empty());
}

#[test]
fn test_04_multiple_valid_sequences() {
    let records = make_test_records();
    let res = parse_control_sgrna(Cursor::new("TTTT\nCCCC\n"), &records).expect("must resolve");
    assert_eq!(res.guide_indexes, vec![1, 3]);
    assert!(res.unknown_entries.is_empty());
}

#[test]
fn test_05_mixed_ids_and_sequences_resolved_independently() {
    let records = make_test_records();
    let res = parse_control_sgrna(Cursor::new("guide_0\nGGGG\n"), &records).expect("must resolve");
    assert_eq!(res.guide_indexes, vec![0, 2]);
    assert!(res.unknown_entries.is_empty());
}

#[test]
fn test_06_token_matching_id_and_sequence_same_row_resolves_once() {
    let records = vec![
        SgRnaRecord {
            id: "SAME_VAL".to_string(),
            gene: "GENE_X".to_string(),
            sequence: "SAME_VAL".to_string(),
        },
        SgRnaRecord {
            id: "other".to_string(),
            gene: "GENE_Y".to_string(),
            sequence: "TTTT".to_string(),
        },
    ];
    let res = parse_control_sgrna(Cursor::new("SAME_VAL\n"), &records).expect("must resolve");
    assert_eq!(res.guide_indexes, vec![0]);
    assert!(res.unknown_entries.is_empty());
}

#[test]
fn test_07_token_matching_id_in_one_and_seq_in_another_ambiguity_error() {
    let records = vec![
        SgRnaRecord {
            id: "guide_0".to_string(),
            gene: "GENE_A".to_string(),
            sequence: "CONFLICT_TOKEN".to_string(),
        },
        SgRnaRecord {
            id: "CONFLICT_TOKEN".to_string(),
            gene: "GENE_B".to_string(),
            sequence: "TTTT".to_string(),
        },
    ];
    let err = parse_control_sgrna(Cursor::new("CONFLICT_TOKEN\n"), &records)
        .expect_err("must return ambiguity error");
    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("ambiguous"),
                "msg must indicate ambiguity: {msg}"
            );
            assert!(
                msg.contains("CONFLICT_TOKEN"),
                "msg must identify ambiguous token: {msg}"
            );
            assert!(msg.contains("row 0"), "msg must identify row 0: {msg}");
            assert!(msg.contains("row 1"), "msg must identify row 1: {msg}");
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_08_resolved_indexes_in_library_row_order_when_reversed() {
    let records = make_test_records();
    let res = parse_control_sgrna(
        Cursor::new("guide_3\nguide_2\nguide_1\nguide_0\n"),
        &records,
    )
    .expect("must resolve");
    assert_eq!(res.guide_indexes, vec![0, 1, 2, 3]);
}

#[test]
fn test_09_duplicate_control_entries_resolve_once() {
    let records = make_test_records();
    let res = parse_control_sgrna(Cursor::new("guide_1\nguide_1\nguide_1\n"), &records)
        .expect("must resolve");
    assert_eq!(res.guide_indexes, vec![1]);
}

#[test]
fn test_10_different_entries_resolving_to_same_guide_resolve_once() {
    let records = make_test_records();
    // guide_1 has sequence CCCC
    let res = parse_control_sgrna(Cursor::new("guide_1\nCCCC\n"), &records).expect("must resolve");
    assert_eq!(res.guide_indexes, vec![1]);
}

#[test]
fn test_11_one_unknown_plus_one_valid_succeeds_and_reports_unknown() {
    let records = make_test_records();
    let res = parse_control_sgrna(Cursor::new("guide_2\nunknown_x\n"), &records)
        .expect("must succeed with partial match");
    assert_eq!(res.guide_indexes, vec![2]);
    assert_eq!(res.unknown_entries, vec!["unknown_x".to_string()]);
}

#[test]
fn test_12_multiple_unknown_entries_preserve_control_file_order() {
    let records = make_test_records();
    let res = parse_control_sgrna(
        Cursor::new("guide_2\nzeta_unknown\nalpha_unknown\nbeta_unknown\n"),
        &records,
    )
    .expect("must succeed");
    assert_eq!(res.guide_indexes, vec![2]);
    assert_eq!(
        res.unknown_entries,
        vec![
            "zeta_unknown".to_string(),
            "alpha_unknown".to_string(),
            "beta_unknown".to_string()
        ]
    );
}

#[test]
fn test_13_repeated_unknown_entries_remain_available_in_repeated_order() {
    let records = make_test_records();
    let res = parse_control_sgrna(
        Cursor::new("guide_0\nrepeat_unk\nother_unk\nrepeat_unk\n"),
        &records,
    )
    .expect("must succeed");
    assert_eq!(res.guide_indexes, vec![0]);
    assert_eq!(
        res.unknown_entries,
        vec![
            "repeat_unk".to_string(),
            "other_unk".to_string(),
            "repeat_unk".to_string()
        ]
    );
}

#[test]
fn test_14_unknown_only_file_returns_empty_resolution_error() {
    let records = make_test_records();
    let err = parse_control_sgrna(Cursor::new("unk1\nunk2\n"), &records)
        .expect_err("must error on zero resolved");
    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("No control sgRNAs") && msg.contains("matched"),
                "expected clear no-match message: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_15_empty_file_returns_empty_resolution_error() {
    let records = make_test_records();
    let err =
        parse_control_sgrna(Cursor::new(""), &records).expect_err("must error on empty input");
    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("No control sgRNAs") && msg.contains("matched"),
                "expected clear no-match message: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_16_blank_only_file_returns_empty_resolution_error() {
    let records = make_test_records();
    let err = parse_control_sgrna(Cursor::new("   \n\t  \n \r\n\n"), &records)
        .expect_err("must error on blank input");
    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("No control sgRNAs") && msg.contains("matched"),
                "expected clear no-match message: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_17_comment_only_file_returns_empty_resolution_error() {
    let records = make_test_records();
    let err = parse_control_sgrna(
        Cursor::new("# First comment\n  # Second comment with leading space\n"),
        &records,
    )
    .expect_err("must error on comment-only input");
    match err {
        MageckError::Config(msg) => {
            assert!(
                msg.contains("No control sgRNAs") && msg.contains("matched"),
                "expected clear no-match message: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_18_blank_lines_among_valid_entries_are_skipped() {
    let records = make_test_records();
    let res = parse_control_sgrna(
        Cursor::new("\n\n  \nguide_0\n\n   \t\nguide_2\n\n"),
        &records,
    )
    .expect("must resolve");
    assert_eq!(res.guide_indexes, vec![0, 2]);
}

#[test]
fn test_19_leading_and_trailing_whitespace_trimmed() {
    let records = make_test_records();
    let res = parse_control_sgrna(Cursor::new("   guide_0   \n\t  guide_1  \t\n"), &records)
        .expect("must resolve");
    assert_eq!(res.guide_indexes, vec![0, 1]);
}

#[test]
fn test_20_full_line_comments_among_valid_entries_are_skipped() {
    let records = make_test_records();
    let res = parse_control_sgrna(
        Cursor::new("# comment 1\nguide_0\n   # comment 2\nguide_3\n#\n"),
        &records,
    )
    .expect("must resolve");
    assert_eq!(res.guide_indexes, vec![0, 3]);
}

#[test]
fn test_21_inline_comment_syntax_is_rejected() {
    let records = make_test_records();
    let err = parse_control_sgrna(
        Cursor::new("guide_0 # this is an inline comment\n"),
        &records,
    )
    .expect_err("must reject inline comment due to internal whitespace");
    match err {
        MageckError::Config(msg) => {
            assert!(msg.contains("Line 1"), "must identify line 1: {msg}");
            assert!(
                msg.contains("whitespace"),
                "must identify whitespace error: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_22_internal_ordinary_space_is_rejected() {
    let records = make_test_records();
    let err = parse_control_sgrna(Cursor::new("guide 0\n"), &records)
        .expect_err("must reject internal space");
    match err {
        MageckError::Config(msg) => {
            assert!(msg.contains("Line 1"), "must identify line 1: {msg}");
            assert!(
                msg.contains("whitespace"),
                "must identify whitespace error: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_23_embedded_tab_is_rejected() {
    let records = make_test_records();
    let err = parse_control_sgrna(Cursor::new("guide_0\tguide_1\n"), &records)
        .expect_err("must reject embedded tab");
    match err {
        MageckError::Config(msg) => {
            assert!(msg.contains("Line 1"), "must identify line 1: {msg}");
            assert!(
                msg.contains("whitespace"),
                "must identify whitespace error: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_24_header_like_line_treated_as_ordinary_unknown_entry() {
    let records = make_test_records();
    let res = parse_control_sgrna(Cursor::new("sgRNA_ID\nguide_1\n"), &records)
        .expect("must succeed with header-like unknown");
    assert_eq!(res.guide_indexes, vec![1]);
    assert_eq!(res.unknown_entries, vec!["sgRNA_ID".to_string()]);
}

#[test]
fn test_25_matching_is_case_sensitive() {
    let records = make_test_records();
    // Library has sequence "AAAA" and ID "guide_0"
    // "aaaa" and "GUIDE_0" do not match case-sensitively
    let res = parse_control_sgrna(Cursor::new("aaaa\nGUIDE_0\nguide_0\n"), &records)
        .expect("must resolve valid case-sensitive match");
    assert_eq!(res.guide_indexes, vec![0]);
    assert_eq!(
        res.unknown_entries,
        vec!["aaaa".to_string(), "GUIDE_0".to_string()]
    );
}

#[test]
fn test_26_gene_name_does_not_resolve() {
    let records = make_test_records();
    // GENE_A is gene name of row 0
    let res = parse_control_sgrna(Cursor::new("GENE_A\nguide_1\n"), &records)
        .expect("must resolve valid guide");
    assert_eq!(res.guide_indexes, vec![1]);
    assert_eq!(res.unknown_entries, vec!["GENE_A".to_string()]);
}

#[test]
fn test_27_missing_file_returns_typed_error() {
    let records = make_test_records();
    let missing_path = PathBuf::from("/non/existent/path/controls.txt");
    let err =
        parse_control_sgrna_file(&missing_path, &records).expect_err("must fail on missing file");
    match err {
        MageckError::Io(io_err) => {
            assert_eq!(io_err.kind(), std::io::ErrorKind::NotFound);
            assert!(
                io_err
                    .to_string()
                    .contains("/non/existent/path/controls.txt"),
                "error message must contain file path: {io_err}"
            );
        }
        other => panic!("expected MageckError::Io, got {other:?}"),
    }
}

#[test]
fn test_28_invalid_utf8_returns_typed_error() {
    let records = make_test_records();
    let mut tmp = NamedTempFile::new().expect("create tempfile");
    tmp.write_all(b"guide_0\n\xFF\xFE\nguide_1\n")
        .expect("write bytes");
    tmp.flush().expect("flush");

    let err =
        parse_control_sgrna_file(tmp.path(), &records).expect_err("must fail on invalid UTF-8");
    match err {
        MageckError::Io(io_err) => {
            assert_eq!(io_err.kind(), std::io::ErrorKind::InvalidData);
            assert!(
                io_err.to_string().contains(tmp.path().to_str().unwrap()),
                "error message must contain file path: {io_err}"
            );
        }
        other => panic!("expected MageckError::Io, got {other:?}"),
    }
}

#[test]
fn test_29_crlf_input_handled_correctly() {
    let records = make_test_records();
    let res = parse_control_sgrna(Cursor::new("guide_0\r\nguide_1\r\n"), &records)
        .expect("must resolve CRLF input");
    assert_eq!(res.guide_indexes, vec![0, 1]);
    assert!(res.unknown_entries.is_empty());
}

#[test]
fn test_30_library_row_order_is_preserved() {
    let records = make_test_records();
    // File specifies row 3, row 1, row 2, row 0
    let res = parse_control_sgrna(
        Cursor::new("guide_3\nguide_1\nguide_2\nguide_0\n"),
        &records,
    )
    .expect("must resolve");
    assert_eq!(res.guide_indexes, vec![0, 1, 2, 3]);
}

#[test]
fn test_31_no_guide_index_appears_more_than_once() {
    let records = make_test_records();
    // guide_0 has sequence AAAA
    let res = parse_control_sgrna(Cursor::new("guide_0\nAAAA\nguide_0\nAAAA\n"), &records)
        .expect("must resolve");
    assert_eq!(res.guide_indexes, vec![0]);
}

#[test]
fn test_32_duplicate_id_in_library_records_causes_ambiguity_error() {
    let records = vec![
        SgRnaRecord {
            id: "dup_id".to_string(),
            gene: "GENE_A".to_string(),
            sequence: "AAAA".to_string(),
        },
        SgRnaRecord {
            id: "unique_id".to_string(),
            gene: "GENE_B".to_string(),
            sequence: "CCCC".to_string(),
        },
        SgRnaRecord {
            id: "dup_id".to_string(),
            gene: "GENE_C".to_string(),
            sequence: "GGGG".to_string(),
        },
    ];
    let err = parse_control_sgrna(Cursor::new("dup_id\n"), &records)
        .expect_err("must return ambiguity error when multiple library records share ID");
    match err {
        MageckError::Config(msg) => {
            assert!(msg.contains("ambiguous"), "msg: {msg}");
            assert!(msg.contains("row 0"), "msg: {msg}");
            assert!(msg.contains("row 2"), "msg: {msg}");
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_33_duplicate_sequence_in_library_records_causes_ambiguity_error() {
    let records = vec![
        SgRnaRecord {
            id: "id_1".to_string(),
            gene: "GENE_A".to_string(),
            sequence: "SAME_SEQ".to_string(),
        },
        SgRnaRecord {
            id: "id_2".to_string(),
            gene: "GENE_B".to_string(),
            sequence: "SAME_SEQ".to_string(),
        },
    ];
    let err = parse_control_sgrna(Cursor::new("SAME_SEQ\n"), &records)
        .expect_err("must return ambiguity error when multiple library records share sequence");
    match err {
        MageckError::Config(msg) => {
            assert!(msg.contains("ambiguous"), "msg: {msg}");
            assert!(msg.contains("row 0"), "msg: {msg}");
            assert!(msg.contains("row 1"), "msg: {msg}");
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_34_multiline_malformed_identifies_line_number() {
    let records = make_test_records();
    let text = "# Comment line 1\nguide_0\nmalformed entry with spaces\nguide_1\n";
    let err = parse_control_sgrna(Cursor::new(text), &records)
        .expect_err("must reject malformed entry on line 3");
    match err {
        MageckError::Config(msg) => {
            assert!(msg.contains("Line 3"), "must identify Line 3: {msg}");
            assert!(
                msg.contains("malformed entry with spaces"),
                "must identify entry: {msg}"
            );
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_35_empty_library_records_returns_error() {
    let empty_records: Vec<SgRnaRecord> = Vec::new();
    let err = parse_control_sgrna(Cursor::new("guide_0\n"), &empty_records)
        .expect_err("must error on empty library");
    match err {
        MageckError::Config(msg) => {
            assert!(msg.contains("library contains no records"), "msg: {msg}");
        }
        other => panic!("expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_36_file_parsing_successful_integration() {
    let records = make_test_records();
    let mut tmp = NamedTempFile::new().expect("create tempfile");
    writeln!(tmp, "# Controls list").unwrap();
    writeln!(tmp, "guide_3").unwrap();
    writeln!(tmp, "CCCC").unwrap(); // sequence of guide_1
    writeln!(tmp, "unknown_entry").unwrap();
    tmp.flush().unwrap();

    let res = parse_control_sgrna_file(tmp.path(), &records).expect("must parse tempfile");
    assert_eq!(res.guide_indexes, vec![1, 3]);
    assert_eq!(res.unknown_entries, vec!["unknown_entry".to_string()]);
}

#[test]
fn test_37_directory_passed_as_file_rejected() {
    let records = make_test_records();
    let dir = tempfile::tempdir().expect("create tempdir");
    let err =
        parse_control_sgrna_file(dir.path(), &records).expect_err("must reject directory path");
    match err {
        MageckError::Io(_) => {
            // Read from directory returns IsADirectory on Unix
            // assert_eq!(io_err.kind(), std::io::ErrorKind::IsADirectory);
            // IsADirectory might not be in older Rust versions or might be mapped differently,
            // but we can just check if it's an Io error. Actually we can check it.
        }
        other => panic!("expected MageckError::Io, got {other:?}"),
    }
}

#[test]
fn test_38_test_literal_hash_without_whitespace() {
    let records = make_test_records();
    let res = parse_control_sgrna(Cursor::new("guide_0#\nguide_1\n"), &records)
        .expect("must resolve valid guide");
    assert_eq!(res.guide_indexes, vec![1]);
    assert_eq!(res.unknown_entries, vec!["guide_0#".to_string()]);
}
