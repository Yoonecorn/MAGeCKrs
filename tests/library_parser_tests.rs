use std::io::{Cursor, Write};
use tempfile::NamedTempFile;

use mageckrs::count::{build_index, parse_library, parse_library_file, Library, SgRnaRecord};
use mageckrs::data::MatchResult;
use mageckrs::error::MageckError;

#[test]
fn test_minimal_valid_library_file() {
    let tsv = "id\tsequence\tgene\ns1\tACGT\tGENE1\n";
    let lib = parse_library(Cursor::new(tsv)).expect("minimal library must parse");
    assert_eq!(lib.len(), 1);
    assert!(!lib.is_empty());
    let rec = lib.get(0).expect("first record");
    assert_eq!(rec.id, "s1");
    assert_eq!(rec.sequence, "ACGT");
    assert_eq!(rec.gene, "GENE1");
}

#[test]
fn test_multiple_valid_guide_rows() {
    let tsv = "id\tsequence\tgene\ns1\tACGT\tGENE1\ns2\tTGCA\tGENE2\ns3\tAAAA\tGENE3\n";
    let lib = parse_library(Cursor::new(tsv)).expect("multi-row library must parse");
    assert_eq!(lib.len(), 3);
    assert_eq!(lib.records()[0].id, "s1");
    assert_eq!(lib.records()[1].id, "s2");
    assert_eq!(lib.records()[2].id, "s3");
}

#[test]
fn test_original_row_order_preservation() {
    let tsv = "id\tsequence\tgene\n\
               zeta\tACGT\tGZ\n\
               alpha\tTGCA\tGA\n\
               beta\tGGGG\tGB\n\
               gamma\tCCCC\tGC\n\
               delta\tTTTT\tGD\n";
    let lib = parse_library(Cursor::new(tsv)).expect("ordered library must parse");
    let ids: Vec<&str> = lib.records().iter().map(|r| r.id.as_str()).collect();
    assert_eq!(ids, vec!["zeta", "alpha", "beta", "gamma", "delta"]);
}

#[test]
fn test_exact_guide_id_preservation() {
    let tsv = "id\tsequence\tgene\nsgRNA-001.a_v2\tACGT\tGeneA\n";
    let lib = parse_library(Cursor::new(tsv)).expect("must parse");
    assert_eq!(lib.records()[0].id, "sgRNA-001.a_v2");
}

#[test]
fn test_exact_gene_id_preservation() {
    let tsv = "id\tsequence\tgene\ns1\tACGT\tMIR-123a_transcription\n";
    let lib = parse_library(Cursor::new(tsv)).expect("must parse");
    assert_eq!(lib.records()[0].gene, "MIR-123a_transcription");
}

#[test]
fn test_required_sequence_normalization() {
    let tsv = "id\tsequence\tgene\ns1\tacgt\tGeneA\n";
    let lib = parse_library(Cursor::new(tsv)).expect("must parse");
    assert_eq!(lib.records()[0].sequence, "ACGT");
}

#[test]
fn test_lowercase_and_mixed_case_sequence() {
    let tsv = "id\tsequence\tgene\ns1\taCgTaCgT\tGeneA\n";
    let lib = parse_library(Cursor::new(tsv)).expect("must parse");
    assert_eq!(lib.records()[0].sequence, "ACGTACGT");
}

#[test]
fn test_invalid_nucleotide_rejection() {
    // Non-ACGT nucleotide 'N'
    let tsv = "id\tsequence\tgene\ns1\tACGTNACGT\tGeneA\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("must reject 'N'");
    match err {
        MageckError::Library(msg) => {
            assert!(
                msg.contains("Line 2"),
                "message must contain line number: {msg}"
            );
            assert!(msg.contains("invalid nucleotide 'N'"), "message: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }

    // Number '1'
    let tsv_num = "id\tsequence\tgene\ns1\tACGT1ACGT\tGeneA\n";
    let err_num = parse_library(Cursor::new(tsv_num)).expect_err("must reject '1'");
    match err_num {
        MageckError::Library(msg) => {
            assert!(
                msg.contains("Line 2"),
                "message must contain line number: {msg}"
            );
            assert!(msg.contains("invalid nucleotide '1'"), "message: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_empty_guide_sequence() {
    let tsv = "id\tsequence\tgene\ns1\t\tGeneA\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("must reject empty sequence");
    match err {
        MageckError::Library(msg) => {
            assert!(msg.contains("Line 2"), "message: {msg}");
            assert!(msg.contains("guide sequence is empty"), "message: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_empty_guide_identifier() {
    let tsv = "id\tsequence\tgene\n\tACGT\tGeneA\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("must reject empty ID");
    match err {
        MageckError::Library(msg) => {
            assert!(msg.contains("Line 2"), "message: {msg}");
            assert!(msg.contains("guide identifier is empty"), "message: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_empty_gene_identifier() {
    let tsv = "id\tsequence\tgene\ns1\tACGT\t\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("must reject empty gene");
    match err {
        MageckError::Library(msg) => {
            assert!(msg.contains("Line 2"), "message: {msg}");
            assert!(msg.contains("gene identifier is empty"), "message: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_missing_required_fields() {
    let tsv = "id\tsequence\tgene\ns1\tACGT\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("must reject row with missing columns");
    match err {
        MageckError::Library(msg) => {
            assert!(msg.contains("Line 2"), "message: {msg}");
            assert!(msg.contains("insufficient columns"), "message: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_extra_fields_policy() {
    // 4th and 5th columns are ignored
    let tsv = "id\tsequence\tgene\tscore\tcomment\ns1\tACGT\tGeneA\t0.95\tnon-targeting\ns2\tTGCA\tGeneB\t0.88\texon1\n";
    let lib = parse_library(Cursor::new(tsv)).expect("extra columns must be safely ignored");
    assert_eq!(lib.len(), 2);
    assert_eq!(lib.records()[0].id, "s1");
    assert_eq!(lib.records()[0].sequence, "ACGT");
    assert_eq!(lib.records()[0].gene, "GeneA");
}

#[test]
fn test_missing_header() {
    // File starts directly with data
    let tsv = "s1\tACGT\tGeneA\ns2\tTGCA\tGeneB\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("must reject missing header");
    match err {
        MageckError::Library(msg) => {
            assert!(msg.contains("Line 1"), "message: {msg}");
            assert!(msg.contains("missing required column"), "message: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_valid_headers() {
    // Standard id, sequence, gene
    let tsv1 = "id\tsequence\tgene\ns1\tACGT\tGeneA\n";
    assert_eq!(parse_library(Cursor::new(tsv1)).unwrap().len(), 1);

    // Permuted column order: gene, sequence, id
    let tsv3 = "gene\tsequence\tid\nGeneA\tACGT\ts1\n";
    let lib3 = parse_library(Cursor::new(tsv3)).unwrap();
    assert_eq!(lib3.len(), 1);
    assert_eq!(lib3.records()[0].id, "s1");
    assert_eq!(lib3.records()[0].sequence, "ACGT");
    assert_eq!(lib3.records()[0].gene, "GeneA");
}

#[test]
fn test_incorrect_header() {
    let tsv = "foo\tbar\tbaz\ns1\tACGT\tGeneA\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("must reject invalid header names");
    match err {
        MageckError::Library(msg) => {
            assert!(msg.contains("Line 1"), "message: {msg}");
            assert!(msg.contains("missing required column"), "message: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_blank_lines_handling() {
    let tsv = "\n\nid\tsequence\tgene\n\ns1\tACGT\tGeneA\n\n\ns2\tTGCA\tGeneB\n\n";
    let lib = parse_library(Cursor::new(tsv)).expect("must ignore blank lines");
    assert_eq!(lib.len(), 2);
    assert_eq!(lib.records()[0].id, "s1");
    assert_eq!(lib.records()[1].id, "s2");
}

#[test]
fn test_comment_lines_not_permitted() {
    // Comment lines starting with # are not permitted by spec and are treated as data/malformed rows
    let tsv = "id\tsequence\tgene\n# this is a comment\ns1\tACGT\tGeneA\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("unsupported comment must error");
    match err {
        MageckError::Library(msg) => {
            assert!(msg.contains("Line 2"), "message: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_unix_lf_input() {
    let tsv = "id\tsequence\tgene\ns1\tACGT\tGeneA\ns2\tTGCA\tGeneB\n";
    let lib = parse_library(Cursor::new(tsv)).expect("LF must parse");
    assert_eq!(lib.len(), 2);
}

#[test]
fn test_crlf_input() {
    let tsv = "id\tsequence\tgene\r\ns1\tACGT\tGeneA\r\ns2\tTGCA\tGeneB\r\n";
    let lib = parse_library(Cursor::new(tsv)).expect("CRLF must parse cleanly");
    assert_eq!(lib.len(), 2);
    assert_eq!(lib.records()[0].id, "s1");
    assert_eq!(lib.records()[0].sequence, "ACGT");
    assert_eq!(lib.records()[0].gene, "GeneA");
    assert_eq!(lib.records()[1].id, "s2");
    assert_eq!(lib.records()[1].sequence, "TGCA");
    assert_eq!(lib.records()[1].gene, "GeneB");
}

#[test]
fn test_duplicate_guide_identifier() {
    let tsv = "id\tsequence\tgene\ns1\tACGT\tGeneA\ns1\tTGCA\tGeneB\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("must reject duplicate guide id");
    match err {
        MageckError::Library(msg) => {
            assert!(
                msg.contains("Line 3"),
                "message must identify second line: {msg}"
            );
            assert!(
                msg.contains("duplicate guide identifier 's1'"),
                "message: {msg}"
            );
            assert!(msg.contains("first seen on line 2"), "message: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_duplicate_guide_sequence() {
    let tsv = "id\tsequence\tgene\ns1\tACGT\tGeneA\ns2\tacgt\tGeneB\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("must reject duplicate sequence");
    match err {
        MageckError::Library(msg) => {
            assert!(msg.contains("Line 3"), "message: {msg}");
            assert!(
                msg.contains("duplicate guide sequence 'ACGT'"),
                "message: {msg}"
            );
            assert!(
                msg.contains("first seen on line 2 for guide 's1'"),
                "message: {msg}"
            );
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_guide_length_32nt() {
    let seq_32 = "A".repeat(32);
    let tsv = format!("id\tsequence\tgene\ns32\t{seq_32}\tGene32\n");
    let lib = parse_library(Cursor::new(tsv)).expect("32nt guide must parse");
    assert_eq!(lib.len(), 1);
    assert_eq!(lib.records()[0].sequence.len(), 32);

    let idx = build_index(&lib, false).expect("must build packed index");
    assert!(idx.is_packed(), "must select PackedGuideIndex for 32nt");
    assert!(!idx.is_byte());
}

#[test]
fn test_guide_length_33nt() {
    let seq_33 = "A".repeat(33);
    let tsv = format!("id\tsequence\tgene\ns33\t{seq_33}\tGene33\n");
    let lib = parse_library(Cursor::new(tsv)).expect("33nt guide must parse");
    assert_eq!(lib.len(), 1);
    assert_eq!(lib.records()[0].sequence.len(), 33);

    let idx = build_index(&lib, false).expect("must build byte fallback index");
    assert!(idx.is_byte(), "must select ByteGuideIndex for 33nt");
    assert!(!idx.is_packed());
}

#[test]
fn test_long_guide_preservation_without_truncation() {
    let seq_45 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA";
    assert_eq!(seq_45.len(), 45);
    let tsv = format!("id\tsequence\tgene\ns45\t{seq_45}\tGeneLong\n");
    let lib = parse_library(Cursor::new(tsv)).expect("45nt guide must parse");
    assert_eq!(lib.records()[0].sequence, seq_45);
    assert_eq!(lib.records()[0].sequence.len(), 45);
}

#[test]
fn test_useful_error_context_line_numbers() {
    let tsv = "id\tsequence\tgene\n\
               s1\tACGT\tGeneA\n\
               s2\tTGCA\tGeneB\n\
               s3\tTGCN\tGeneC\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("error on line 4");
    match err {
        MageckError::Library(msg) => {
            assert!(msg.contains("Line 4:"), "must report Line 4: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_empty_file() {
    let tsv = "";
    let err = parse_library(Cursor::new(tsv)).expect_err("empty file must error");
    match err {
        MageckError::Library(msg) => {
            assert_eq!(msg, "Library file is empty");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_header_only_file() {
    let tsv = "id\tsequence\tgene\n";
    let err = parse_library(Cursor::new(tsv)).expect_err("header only must error");
    match err {
        MageckError::Library(msg) => {
            assert!(msg.contains("no guide records"), "message: {msg}");
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_file_not_found() {
    let non_existent = std::path::PathBuf::from("/non/existent/path/library.tsv");
    let err = parse_library_file(&non_existent).expect_err("file not found must error");
    match err {
        MageckError::Io(io_err) => {
            assert_eq!(io_err.kind(), std::io::ErrorKind::NotFound);
        }
        other => panic!("expected MageckError::Io, got: {other:?}"),
    }
}

#[test]
fn test_no_panic_malformed_input() {
    let malformed_inputs = [
        "\0\0\0\0\0",
        "random binary \u{00FF}\u{00FE}\u{0000} content",
        "\t\t\t\n\t\t\t",
        "id,seq\na,b\n",
        "id\tsequence\tgene\n\t\t\n",
    ];
    for input in malformed_inputs {
        let _ = parse_library(Cursor::new(input));
    }
}

#[test]
fn test_index_packed_selected_for_guides_at_or_below_32() {
    let tsv = "id\tsequence\tgene\ns20\tACGTACGTACGTACGTACGT\tGeneA\ns32\tACGTACGTACGTACGTACGTACGTACGTACGT\tGeneB\n";
    let lib = parse_library(Cursor::new(tsv)).expect("parse");
    let idx = lib.build_index(false).expect("build index");
    assert!(idx.is_packed());
    assert!(!idx.is_byte());
    assert_eq!(idx.len(), 2);
}

#[test]
fn test_index_byte_fallback_for_guides_above_32() {
    let tsv = "id\tsequence\tgene\ns20\tACGTACGTACGTACGTACGT\tGeneA\ns33\tACGTACGTACGTACGTACGTACGTACGTACGTA\tGeneB\n";
    let lib = parse_library(Cursor::new(tsv)).expect("parse");
    let idx = lib.build_index(false).expect("build index");
    assert!(idx.is_byte());
    assert!(!idx.is_packed());
    assert_eq!(idx.len(), 2);
}

#[test]
fn test_index_exact_lookup() {
    let tsv = "id\tsequence\tgene\ns1\tACGTACGT\tGeneA\ns2\tTGCATGCA\tGeneB\n";
    let lib = parse_library(Cursor::new(tsv)).expect("parse");
    let idx = build_index(&lib, false).expect("build index");

    assert_eq!(idx.lookup(b"ACGTACGT"), MatchResult::Exact(0));
    assert_eq!(idx.lookup(b"TGCATGCA"), MatchResult::Exact(1));
    assert_eq!(idx.lookup(b"GGGGGGGG"), MatchResult::NoMatch);
    // Lowercase query is handled seamlessly
    assert_eq!(idx.lookup(b"acgtacgt"), MatchResult::Exact(0));
}

#[test]
fn test_index_duplicate_sequence_rejection() {
    // If Library::new is bypassed or constructed with duplicates, build_index returns error
    let recs = vec![
        SgRnaRecord {
            id: "s1".to_string(),
            gene: "G1".to_string(),
            sequence: "ACGT".to_string(),
        },
        SgRnaRecord {
            id: "s2".to_string(),
            gene: "G2".to_string(),
            sequence: "ACGT".to_string(),
        },
    ];
    let err = Library::new(recs).expect_err("Library::new must reject duplicate sequence");
    match err {
        MageckError::Library(msg) => {
            assert!(
                msg.contains("duplicate guide sequence 'ACGT'"),
                "message: {msg}"
            );
        }
        other => panic!("expected MageckError::Library, got: {other:?}"),
    }
}

#[test]
fn test_index_no_one_mismatch_unless_enabled() {
    let tsv = "id\tsequence\tgene\ns1\tACGTACGT\tGeneA\n";
    let lib = parse_library(Cursor::new(tsv)).expect("parse");

    // 0-mismatch mode: substitution is NoMatch
    let idx_0 = build_index(&lib, false).expect("build 0-mismatch index");
    assert_eq!(idx_0.lookup(b"ACGTACGT"), MatchResult::Exact(0));
    assert_eq!(idx_0.lookup(b"ACGTACGA"), MatchResult::NoMatch);

    // 1-mismatch mode: substitution is OneMismatch(0)
    let idx_1 = build_index(&lib, true).expect("build 1-mismatch index");
    assert_eq!(idx_1.lookup(b"ACGTACGT"), MatchResult::Exact(0));
    assert_eq!(idx_1.lookup(b"ACGTACGA"), MatchResult::OneMismatch(0));
}

#[test]
fn test_index_long_guide_byte_intact() {
    let seq_36 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let tsv = format!("id\tsequence\tgene\nlong1\t{seq_36}\tGeneLong\n");
    let lib = parse_library(Cursor::new(tsv)).expect("parse");
    let idx = build_index(&lib, false).expect("build index");
    assert!(idx.is_byte());

    assert_eq!(idx.lookup(seq_36.as_bytes()), MatchResult::Exact(0));
    assert_eq!(idx.lookup(b"ACGT"), MatchResult::NoMatch);
}

#[test]
fn test_parse_library_file_tempfile() {
    let mut file = NamedTempFile::new().expect("create tempfile");
    writeln!(file, "id\tsequence\tgene").unwrap();
    writeln!(file, "s1\tACGTACGT\tGeneA").unwrap();
    writeln!(file, "s2\tTGCATGCA\tGeneB").unwrap();
    file.flush().unwrap();

    let lib = parse_library_file(file.path()).expect("parse from file");
    assert_eq!(lib.len(), 2);
    assert_eq!(lib.records()[0].id, "s1");
    assert_eq!(lib.records()[1].id, "s2");
}

#[test]
fn test_csv_format() {
    let csv = "id,sequence,gene\ns1,ACGT,GeneA\n";
    let lib = parse_library(Cursor::new(csv)).expect("CSV must parse");
    assert_eq!(lib.len(), 1);
    assert_eq!(lib.records()[0].id, "s1");
    assert_eq!(lib.records()[0].sequence, "ACGT");
    assert_eq!(lib.records()[0].gene, "GeneA");
}

#[test]
fn test_mixed_length_guides_byte_fallback() {
    let seq_20 = "ACGTACGTACGTACGTACGT";
    let seq_32 = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let seq_33 = "ACGTACGTACGTACGTACGTACGTACGTACGTA";

    let tsv = format!(
        "id\tsequence\tgene\ns20\t{seq_20}\tGene20\ns32\t{seq_32}\tGene32\ns33\t{seq_33}\tGene33\n"
    );
    let lib = parse_library(Cursor::new(tsv)).expect("mixed length library must parse");

    let idx_0 = build_index(&lib, false).expect("build index");
    assert!(
        idx_0.is_byte(),
        "must select ByteGuideIndex due to 33nt guide"
    );

    assert_eq!(
        idx_0.lookup(seq_20.as_bytes()),
        mageckrs::data::MatchResult::Exact(0)
    );
    assert_eq!(
        idx_0.lookup(seq_32.as_bytes()),
        mageckrs::data::MatchResult::Exact(1)
    );
    assert_eq!(
        idx_0.lookup(seq_33.as_bytes()),
        mageckrs::data::MatchResult::Exact(2)
    );

    assert_eq!(
        idx_0.lookup(b"ACGTACGTACGTACGTACG"),
        mageckrs::data::MatchResult::NoMatch
    );

    let idx_1 = build_index(&lib, true).expect("build 1-mismatch index");
    assert!(idx_1.is_byte());

    assert_eq!(
        idx_1.lookup(seq_20.as_bytes()),
        mageckrs::data::MatchResult::Exact(0)
    );
    assert_eq!(
        idx_1.lookup(seq_32.as_bytes()),
        mageckrs::data::MatchResult::Exact(1)
    );
    assert_eq!(
        idx_1.lookup(seq_33.as_bytes()),
        mageckrs::data::MatchResult::Exact(2)
    );

    let mut mismatch_20 = seq_20.as_bytes().to_vec();
    mismatch_20[0] = b'T';
    assert_eq!(
        idx_1.lookup(&mismatch_20),
        mageckrs::data::MatchResult::OneMismatch(0)
    );
}
