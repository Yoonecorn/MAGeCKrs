use std::io::{Cursor, Write};
use std::path::PathBuf;
use tempfile::NamedTempFile;

use mageckrs::config::FastqSampleGroup;
use mageckrs::count::matcher::{build_index, parse_library, Library};
use mageckrs::count::reader::{
    count_all_fastq, count_grouped_fastq, process_one_fastq, GroupedCountResult,
    PhysicalFastqSummary, ReadSummary,
};

use mageckrs::error::MageckError;

/// Helper to write temporary FASTQ content.
fn make_temp_fastq(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().expect("create temp fastq");
    file.write_all(content.as_bytes())
        .expect("write temp fastq");
    file.flush().expect("flush temp fastq");
    file
}

/// Helper to write raw bytes to a temporary file (e.g. for gzip).
fn make_temp_fastq_bytes(bytes: &[u8], suffix: &str) -> NamedTempFile {
    let mut file = tempfile::Builder::new()
        .suffix(suffix)
        .tempfile()
        .expect("create temp binary fastq");
    file.write_all(bytes).expect("write temp binary fastq");
    file.flush().expect("flush temp binary fastq");
    file
}

/// Helper to parse a library from a TSV string.
fn make_library(tsv: &str) -> Library {
    parse_library(Cursor::new(tsv)).expect("parse library TSV")
}

const STD_LIB_TSV: &str = "id\tsequence\tgene\n\
g0\tACGTACGTACGTACGTACGT\tGENE_A\n\
g1\tTGCATGCATGCATGCATGCA\tGENE_B\n\
g2\tGGGGCCCCAAAATTTTGGGG\tGENE_C\n";

// Pre-compressed gzip FASTQ containing:
// @read_exact\nACGTACGTACGTACGTACGT\n+\nFFFFFFFFFFFFFFFFFFFF\n
const GZ_FASTQ_BYTES: &[u8] = &[
    31, 139, 8, 0, 98, 130, 152, 106, 2, 255, 115, 40, 74, 77, 76, 137, 79, 173, 72, 76, 46, 225,
    114, 116, 118, 15, 65, 199, 92, 218, 92, 110, 88, 0, 23, 0, 141, 149, 224, 203, 56, 0, 0, 0,
];

#[test]
fn test_one_valid_record() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq = "@read1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, false).expect("process fastq");

    assert_eq!(counts, vec![1.0, 0.0, 0.0]);
    assert_eq!(summary.total_reads, 1);
    assert_eq!(summary.mapped_reads, 1);
    assert_eq!(summary.exact_matches, 1);
    assert_eq!(summary.one_mismatch_matches, 0);
    assert_eq!(summary.ambiguous_reads, 0);
    assert_eq!(summary.unmapped_reads, 0);
    assert_eq!(summary.short_reads, 0);
    assert!((summary.mapped_percentage() - 1.0).abs() < 1e-6);
}

#[test]
fn test_multiple_valid_records() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @r2\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @r3\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @r4\nGGGGCCCCAAAATTTTGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, false).expect("process fastq");

    assert_eq!(counts, vec![2.0, 1.0, 1.0]);
    assert_eq!(summary.total_reads, 4);
    assert_eq!(summary.mapped_reads, 4);
    assert_eq!(summary.exact_matches, 4);
    assert_eq!(summary.one_mismatch_matches, 0);
    assert_eq!(summary.ambiguous_reads, 0);
    assert_eq!(summary.unmapped_reads, 0);
    assert_eq!(summary.short_reads, 0);
}

#[test]
fn test_exact_guide_match() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq = "@r1\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, false).expect("process fastq");

    assert_eq!(counts[1], 1.0);
    assert_eq!(counts[0], 0.0);
    assert_eq!(counts[2], 0.0);
    assert_eq!(summary.exact_matches, 1);
    assert_eq!(summary.mapped_reads, 1);
}

#[test]
fn test_exact_match_precedence_over_mismatch_candidate() {
    // Guide 0 and Guide 1 differ by 1 base:
    // g0: ACGTACGTACGTACGTACGT
    // g1: CCGTACGTACGTACGTACGT (pos 0 is C instead of A)
    let tsv = "id\tsequence\tgene\n\
               g0\tACGTACGTACGTACGTACGT\tGENE_A\n\
               g1\tCCGTACGTACGTACGTACGT\tGENE_B\n";
    let lib = make_library(tsv);
    // Index built with one_mismatch = true
    let index = build_index(&lib, true).expect("build index");

    // Read matches g0 exactly. It is also a 1-mismatch variant of g1.
    // Exact match must unconditionally take precedence and increment g0.
    let fq = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, true).expect("process fastq");

    assert_eq!(counts, vec![1.0, 0.0]);
    assert_eq!(summary.exact_matches, 1);
    assert_eq!(summary.one_mismatch_matches, 0);
    assert_eq!(summary.mapped_reads, 1);
}

#[test]
fn test_one_mismatch_unique_match_when_enabled() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, true).expect("build index with mismatch");

    // Read differs from g0 by 1 base at position 0: T instead of A
    let fq = "@r1\nTCGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, true).expect("process fastq");

    assert_eq!(counts, vec![1.0, 0.0, 0.0]);
    assert_eq!(summary.total_reads, 1);
    assert_eq!(summary.mapped_reads, 1);
    assert_eq!(summary.exact_matches, 0);
    assert_eq!(summary.one_mismatch_matches, 1);
    assert_eq!(summary.unmapped_reads, 0);
}

#[test]
fn test_one_mismatch_disabled() {
    let lib = make_library(STD_LIB_TSV);
    // Index built with mismatch enabled, but one_mismatch flag passed as false to reader
    let index = build_index(&lib, true).expect("build index with mismatch");

    // Read differs from g0 by 1 base at position 0: T instead of A
    let fq = "@r1\nTCGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, false).expect("process fastq");

    assert_eq!(counts, vec![0.0, 0.0, 0.0]);
    assert_eq!(summary.total_reads, 1);
    assert_eq!(summary.mapped_reads, 0);
    assert_eq!(summary.exact_matches, 0);
    assert_eq!(summary.one_mismatch_matches, 0);
    assert_eq!(summary.unmapped_reads, 1);
}

#[test]
fn test_ambiguous_one_mismatch_read_silently_skipped() {
    // Two guides differing by 2 bases:
    // g0: ACGTACGTACGTACGTACGA (ends in A)
    // g1: ACGTACGTACGTACGTACGG (ends in G)
    // A read ending in C (ACGTACGTACGTACGTACGC) is 1 mismatch away from both g0 and g1.
    let tsv = "id\tsequence\tgene\n\
               g0\tACGTACGTACGTACGTACGA\tGENE_A\n\
               g1\tACGTACGTACGTACGTACGG\tGENE_B\n";
    let lib = make_library(tsv);
    let index = build_index(&lib, true).expect("build index");

    let fq = "@ambiguous\nACGTACGTACGTACGTACGC\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, true).expect("process fastq");

    // Must be silently skipped: no count increment for either guide!
    assert_eq!(counts, vec![0.0, 0.0]);
    assert_eq!(summary.total_reads, 1);
    assert_eq!(summary.mapped_reads, 0);
    assert_eq!(summary.ambiguous_reads, 1);
    assert_eq!(summary.unmapped_reads, 0);
}

#[test]
fn test_no_match_read() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // Read contains sequence not present in library
    let fq = "@nomatch\nTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, false).expect("process fastq");

    assert_eq!(counts, vec![0.0, 0.0, 0.0]);
    assert_eq!(summary.total_reads, 1);
    assert_eq!(summary.mapped_reads, 0);
    assert_eq!(summary.unmapped_reads, 1);
}

#[test]
fn test_read_shorter_than_required_extraction_window() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // 15 nt sequence when window requires 20 nt
    let fq = "@short\nACGTACGTACGTACG\n+\nIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, false).expect("process fastq");

    assert_eq!(counts, vec![0.0, 0.0, 0.0]);
    assert_eq!(summary.total_reads, 1);
    assert_eq!(summary.mapped_reads, 0);
    assert_eq!(summary.short_reads, 1);
    assert_eq!(summary.unmapped_reads, 0);
}

#[test]
fn test_guide_window_at_beginning_of_read() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // 25 nt read with guide at beginning (0..20) followed by 5 trailing bases
    let fq = "@r1\nACGTACGTACGTACGTACGTNNNNN\n+\nIIIIIIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, false).expect("process fastq");

    assert_eq!(counts[0], 1.0);
    assert_eq!(summary.exact_matches, 1);
}

#[test]
fn test_guide_window_at_nonzero_offset() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // 5 nt prefix (e.g. barcode/adapter) + 20 nt guide = 25 nt read
    let fq = "@r1\nTATATACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 5, 20, false).expect("process fastq");

    assert_eq!(counts[0], 1.0);
    assert_eq!(summary.exact_matches, 1);
    assert_eq!(summary.mapped_reads, 1);
}

#[test]
fn test_guide_window_ending_at_final_read_base() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // Read length is exactly start (5) + len (20) = 25 nt
    let fq = "@r1\nCCCCCACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 5, 20, false).expect("process fastq");

    assert_eq!(counts[0], 1.0);
    assert_eq!(summary.exact_matches, 1);
}

#[test]
fn test_off_by_one_boundary_behavior() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // start = 5, len = 20 requires 25 nt.
    // Read 1 is 24 nt (off by 1 short -> short_reads)
    // Read 2 is 25 nt (exact length -> matches)
    let fq = "@r_short\nCCCCCACGTACGTACGTACGTACG\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n\
              @r_exact\nCCCCCACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 5, 20, false).expect("process fastq");

    assert_eq!(counts[0], 1.0);
    assert_eq!(summary.total_reads, 2);
    assert_eq!(summary.short_reads, 1);
    assert_eq!(summary.exact_matches, 1);
    assert_eq!(summary.mapped_reads, 1);
}

#[test]
fn test_invalid_offset_or_range_configuration() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    // sgrna_len = 0 must error
    let res_zero_len = process_one_fastq(tmp.path(), &index, 0, 0, false);
    assert!(matches!(res_zero_len, Err(MageckError::Fastq(_))));

    // start + len overflowing usize must error
    let res_overflow = process_one_fastq(tmp.path(), &index, usize::MAX, 20, false);
    assert!(matches!(res_overflow, Err(MageckError::Fastq(_))));
}

#[test]
fn test_malformed_fastq_record() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // Malformed FASTQ: quality line missing or header missing
    let bad_fq = "This is not a valid fastq file header\n";
    let tmp = make_temp_fastq(bad_fq);

    let res = process_one_fastq(tmp.path(), &index, 0, 20, false);
    assert!(
        matches!(res, Err(MageckError::Fastq(_))),
        "Expected MageckError::Fastq on malformed FASTQ"
    );
}

#[test]
fn test_truncated_fastq_record() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // Truncated FASTQ: header present, sequence present, but quality line truncated
    let truncated_fq = "@incomplete_record\nACGTACGTACGTACGTACGT\n+\n";
    let tmp = make_temp_fastq(truncated_fq);

    let res = process_one_fastq(tmp.path(), &index, 0, 20, false);
    assert!(
        matches!(res, Err(MageckError::Fastq(_))),
        "Expected MageckError::Fastq on truncated FASTQ"
    );
}

#[test]
fn test_lowercase_read_sequence_behavior() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // Lowercase read sequence: acgtacgtacgtacgtacgt matches g0
    let fq = "@r_lower\nacgtacgtacgtacgtacgt\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, false).expect("process fastq");

    assert_eq!(counts[0], 1.0);
    assert_eq!(summary.exact_matches, 1);
    assert_eq!(summary.mapped_reads, 1);
}

#[test]
fn test_invalid_read_nucleotide_behavior() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // Read contains 'N' in guide region: cannot match any guide
    let fq = "@r_n\nACGTACGTACGTACGTACGN\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, false).expect("process fastq");

    assert_eq!(counts, vec![0.0, 0.0, 0.0]);
    assert_eq!(summary.total_reads, 1);
    assert_eq!(summary.unmapped_reads, 1);
    assert_eq!(summary.mapped_reads, 0);
}

#[test]
fn test_empty_fastq_file() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let tmp = make_temp_fastq("");

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, false).expect("process fastq");

    assert_eq!(counts, vec![0.0, 0.0, 0.0]);
    assert_eq!(summary, ReadSummary::default());
    assert_eq!(summary.mapped_percentage(), 0.0);
}

#[test]
fn test_file_not_found_behavior() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let non_existent = PathBuf::from("/non/existent/path/to/reads.fastq");
    let res = process_one_fastq(&non_existent, &index, 0, 20, false);

    assert!(
        matches!(res, Err(MageckError::Fastq(_))),
        "Expected MageckError::Fastq on missing file"
    );
}

#[test]
fn test_correct_original_guide_order_counts() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // Reads arriving in reverse order: g2, g2, g0, g1, g2
    let fq = "@r1\nGGGGCCCCAAAATTTTGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @r2\nGGGGCCCCAAAATTTTGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @r3\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @r4\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @r5\nGGGGCCCCAAAATTTTGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, false).expect("process fastq");

    // g0: 1, g1: 1, g2: 3
    assert_eq!(counts, vec![1.0, 1.0, 3.0]);
    assert_eq!(summary.total_reads, 5);
    assert_eq!(summary.mapped_reads, 5);
}

#[test]
fn test_summary_counters_partition_invariant() {
    // Two guides differing by 2 bases at end:
    // g0: ACGTACGTACGTACGTACGA
    // g1: ACGTACGTACGTACGTACGG
    let tsv = "id\tsequence\tgene\n\
               g0\tACGTACGTACGTACGTACGA\tGENE_A\n\
               g1\tACGTACGTACGTACGTACGG\tGENE_B\n";
    let lib = make_library(tsv);
    let index = build_index(&lib, true).expect("build index");

    let fq = "@r_exact\nACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @r_mismatch\nTCGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @r_ambiguous\nACGTACGTACGTACGTACGC\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @r_unmapped\nTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @r_short\nACGTACGT\n+\nIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, true).expect("process fastq");

    assert_eq!(counts[0], 2.0); // exact + 1-mismatch
    assert_eq!(counts[1], 0.0);

    assert_eq!(summary.total_reads, 5);
    assert_eq!(summary.exact_matches, 1);
    assert_eq!(summary.one_mismatch_matches, 1);
    assert_eq!(summary.mapped_reads, 2);
    assert_eq!(summary.ambiguous_reads, 1);
    assert_eq!(summary.unmapped_reads, 1);
    assert_eq!(summary.short_reads, 1);

    // Invariant checks
    assert_eq!(
        summary.mapped_reads,
        summary.exact_matches + summary.one_mismatch_matches
    );
    assert_eq!(
        summary.total_reads,
        summary.mapped_reads
            + summary.ambiguous_reads
            + summary.unmapped_reads
            + summary.short_reads
    );
    assert!((summary.mapped_percentage() - 0.4).abs() < 1e-6);
}

#[test]
fn test_no_count_increment_for_ambiguous_and_nomatch() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, true).expect("build index");

    let fq = "@nomatch1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @nomatch2\nTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 20, true).expect("process fastq");

    assert_eq!(counts, vec![0.0, 0.0, 0.0]);
    assert_eq!(summary.mapped_reads, 0);
    assert_eq!(summary.unmapped_reads, 2);
}

#[test]
fn test_byte_guide_index_longer_than_32nt() {
    // 33 nt guide triggers ByteGuideIndex fallback
    let tsv = "id\tsequence\tgene\n\
               g_long\tACGTACGTACGTACGTACGTACGTACGTACGTA\tGENE_LONG\n";
    let lib = make_library(tsv);
    assert_eq!(lib.records()[0].sequence.len(), 33);

    let index = build_index(&lib, false).expect("build index");
    assert!(index.is_byte());

    let fq = "@r_long\nACGTACGTACGTACGTACGTACGTACGTACGTA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 33, false).expect("process fastq");

    assert_eq!(counts, vec![1.0]);
    assert_eq!(summary.exact_matches, 1);
    assert_eq!(summary.mapped_reads, 1);
}

#[test]
fn test_mixed_length_library_byte_guide_index() {
    // Library with 20 nt guide and 33 nt guide triggers ByteGuideIndex
    let tsv = "id\tsequence\tgene\n\
               g_short\tACGTACGTACGTACGTACGT\tGENE_SHORT\n\
               g_long\tACGTACGTACGTACGTACGTACGTACGTACGTA\tGENE_LONG\n";
    let lib = make_library(tsv);
    let index = build_index(&lib, false).expect("build index");
    assert!(index.is_byte());

    // FASTQ contains one read matching g_long
    let fq = "@r1\nACGTACGTACGTACGTACGTACGTACGTACGTA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let (counts, summary) =
        process_one_fastq(tmp.path(), &index, 0, 33, false).expect("process fastq");

    assert_eq!(counts, vec![0.0, 1.0]);
    assert_eq!(summary.exact_matches, 1);
}

#[test]
fn test_compressed_fastq_streaming() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // Write pre-compressed gzip FASTQ
    let gz_file = make_temp_fastq_bytes(GZ_FASTQ_BYTES, ".fastq.gz");

    let (counts, summary) =
        process_one_fastq(gz_file.path(), &index, 0, 20, false).expect("process gzipped fastq");

    assert_eq!(counts, vec![1.0, 0.0, 0.0]);
    assert_eq!(summary.total_reads, 1);
    assert_eq!(summary.mapped_reads, 1);
    assert_eq!(summary.exact_matches, 1);
}

#[test]
fn test_count_all_fastq_sequential() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq1 = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let fq2 = "@r2\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n";

    let tmp1 = make_temp_fastq(fq1);
    let tmp2 = make_temp_fastq(fq2);

    let (matrix, summaries) = count_all_fastq(&[tmp1.path(), tmp2.path()], &index, 0, 20, false)
        .expect("count all fastq");

    assert_eq!(matrix.n_samples(), 2);
    assert_eq!(matrix.n_guides(), 3);
    assert_eq!(matrix.sample_slice(0), &[1.0, 0.0, 0.0]);
    assert_eq!(matrix.sample_slice(1), &[0.0, 1.0, 0.0]);

    assert_eq!(summaries.len(), 2);
    assert_eq!(summaries[0].exact_matches, 1);
    assert_eq!(summaries[1].exact_matches, 1);
}

#[test]
fn test_read_summary_mapped_percentage_fraction_scale() {
    // 1. zero reads gives 0.0
    let s_zero = ReadSummary {
        total_reads: 0,
        mapped_reads: 0,
        ..Default::default()
    };
    assert_eq!(s_zero.mapped_percentage(), 0.0);

    // 2. zero mapped out of nonzero reads gives 0.0
    let s_zero_mapped = ReadSummary {
        total_reads: 100,
        mapped_reads: 0,
        unmapped_reads: 100,
        ..Default::default()
    };
    assert_eq!(s_zero_mapped.mapped_percentage(), 0.0);

    // 3. one-quarter mapped gives 0.25
    let s_quarter = ReadSummary {
        total_reads: 4,
        mapped_reads: 1,
        unmapped_reads: 3,
        ..Default::default()
    };
    assert_eq!(s_quarter.mapped_percentage(), 0.25);

    // 4. one-half mapped gives 0.5
    let s_half = ReadSummary {
        total_reads: 100,
        mapped_reads: 50,
        unmapped_reads: 50,
        ..Default::default()
    };
    assert_eq!(s_half.mapped_percentage(), 0.5);

    // 5. all reads mapped gives 1.0
    let s_all = ReadSummary {
        total_reads: 100,
        mapped_reads: 100,
        ..Default::default()
    };
    assert_eq!(s_all.mapped_percentage(), 1.0);

    // 6. reference values: 12,052,478 / 36,736,318
    let s_ref = ReadSummary {
        total_reads: 36_736_318,
        mapped_reads: 12_052_478,
        unmapped_reads: 36_736_318 - 12_052_478,
        ..Default::default()
    };
    let ref_frac = s_ref.mapped_percentage();
    assert!((ref_frac - (12_052_478.0 / 36_736_318.0)).abs() < 1e-12);
    assert_eq!(format!("{ref_frac:.4}"), "0.3281");

    // 7. verify mapped_percentage returns fraction in [0, 1], not percentage in [0, 100]
    assert!(s_quarter.mapped_percentage() <= 1.0);
    assert!(s_half.mapped_percentage() <= 1.0);
    assert!(s_all.mapped_percentage() <= 1.0);
    assert!(s_ref.mapped_percentage() <= 1.0);
}

// ─── Phase 3G Grouped FASTQ Counting Tests ───────────────────────────────────

#[test]
fn test_grouped_single_group_one_fastq() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq = "@read1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n\
              @read2\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let group = FastqSampleGroup {
        fastq_paths: vec![tmp.path().to_path_buf()],
        label: "control_group".to_string(),
    };

    let result: GroupedCountResult = count_grouped_fastq(&[group], &index, 0, 20, false)
        .expect("count grouped single group single file");

    // Dimensions: 1 sample group, 3 library guides
    assert_eq!(result.counts.n_samples(), 1);
    assert_eq!(result.counts.n_guides(), 3);
    assert_eq!(result.counts.sample_slice(0), &[2.0, 0.0, 0.0]);

    // Single physical summary
    assert_eq!(result.summaries.len(), 1);
    let s: &PhysicalFastqSummary = &result.summaries[0];
    assert_eq!(s.path, tmp.path());
    assert_eq!(s.label, "control_group");
    assert_eq!(s.summary.total_reads, 2);
    assert_eq!(s.summary.mapped_reads, 2);
    assert_eq!(s.summary.exact_matches, 2);
    assert_eq!(s.summary.one_mismatch_matches, 0);
    assert_eq!(s.zerocounts, 2); // guides 1 and 2 are 0
    assert!((s.summary.mapped_percentage() - 1.0).abs() < 1e-6);
}

#[test]
fn test_grouped_two_groups_one_fastq_each() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq1 = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let fq2 = "@r2\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp1 = make_temp_fastq(fq1);
    let tmp2 = make_temp_fastq(fq2);

    let groups = vec![
        FastqSampleGroup {
            fastq_paths: vec![tmp1.path().to_path_buf()],
            label: "sampleA".to_string(),
        },
        FastqSampleGroup {
            fastq_paths: vec![tmp2.path().to_path_buf()],
            label: "sampleB".to_string(),
        },
    ];

    let result = count_grouped_fastq(&groups, &index, 0, 20, false).expect("count two groups");

    assert_eq!(result.counts.n_samples(), 2);
    assert_eq!(result.counts.n_guides(), 3);

    // Counts remain in separate columns (sample rows in sample-major CountMatrix)
    assert_eq!(result.counts.sample_slice(0), &[1.0, 0.0, 0.0]);
    assert_eq!(result.counts.sample_slice(1), &[0.0, 1.0, 0.0]);

    // Summaries preserve sample-group order
    assert_eq!(result.summaries.len(), 2);
    assert_eq!(result.summaries[0].path, tmp1.path());
    assert_eq!(result.summaries[0].label, "sampleA");
    assert_eq!(result.summaries[0].summary.mapped_reads, 1);
    assert_eq!(result.summaries[0].zerocounts, 2);

    assert_eq!(result.summaries[1].path, tmp2.path());
    assert_eq!(result.summaries[1].label, "sampleB");
    assert_eq!(result.summaries[1].summary.mapped_reads, 1);
    assert_eq!(result.summaries[1].zerocounts, 2);
}

#[test]
fn test_grouped_single_group_two_technical_replicates() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // Rep 1: 2 reads for g0, 1 read for g1 (total 3, zerocounts = 1 for g2)
    let fq1 = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n\
               @r2\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n\
               @r3\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n";
    // Rep 2: 1 read for g1, 2 reads for g2 (total 3, zerocounts = 1 for g0)
    let fq2 = "@r4\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n\
               @r5\nGGGGCCCCAAAATTTTGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n\
               @r6\nGGGGCCCCAAAATTTTGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n";

    let tmp1 = make_temp_fastq(fq1);
    let tmp2 = make_temp_fastq(fq2);

    let group = FastqSampleGroup {
        fastq_paths: vec![tmp1.path().to_path_buf(), tmp2.path().to_path_buf()],
        label: "treated".to_string(),
    };

    let result = count_grouped_fastq(&[group], &index, 0, 20, false)
        .expect("count single group two technical replicates");

    // 1 group -> 1 sample row in CountMatrix
    assert_eq!(result.counts.n_samples(), 1);
    assert_eq!(result.counts.n_guides(), 3);

    // Technical-replicate counts add into one sample column
    // g0 = 2 + 0 = 2.0
    // g1 = 1 + 1 = 2.0
    // g2 = 0 + 2 = 2.0
    assert_eq!(result.counts.sample_slice(0), &[2.0, 2.0, 2.0]);

    // Combined sample row has 0 zero-count guides
    let combined_zero_count = result
        .counts
        .sample_slice(0)
        .iter()
        .filter(|&&c| c == 0.0)
        .count();
    assert_eq!(combined_zero_count, 0);

    // Per-file summaries remain unaggregated and retain supplied paths and shared label
    assert_eq!(result.summaries.len(), 2);

    // Replicate 1 summary
    assert_eq!(result.summaries[0].path, tmp1.path());
    assert_eq!(result.summaries[0].label, "treated");
    assert_eq!(result.summaries[0].summary.total_reads, 3);
    assert_eq!(result.summaries[0].summary.mapped_reads, 3);
    assert_eq!(result.summaries[0].summary.exact_matches, 3);
    assert_eq!(result.summaries[0].zerocounts, 1); // g2 was 0 in rep 1
    assert_ne!(result.summaries[0].zerocounts, combined_zero_count);

    // Replicate 2 summary
    assert_eq!(result.summaries[1].path, tmp2.path());
    assert_eq!(result.summaries[1].label, "treated");
    assert_eq!(result.summaries[1].summary.total_reads, 3);
    assert_eq!(result.summaries[1].summary.mapped_reads, 3);
    assert_eq!(result.summaries[1].summary.exact_matches, 3);
    assert_eq!(result.summaries[1].zerocounts, 1); // g0 was 0 in rep 2
    assert_ne!(result.summaries[1].zerocounts, combined_zero_count);

    // File order within group preserved
    assert_eq!(result.summaries[0].path, tmp1.path());
    assert_eq!(result.summaries[1].path, tmp2.path());
}

#[test]
fn test_grouped_two_groups_two_replicates_each() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq_a1 = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let fq_a2 = "@r2\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let fq_b1 = "@r3\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let fq_b2 = "@r4\nGGGGCCCCAAAATTTTGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n";

    let tmp_a1 = make_temp_fastq(fq_a1);
    let tmp_a2 = make_temp_fastq(fq_a2);
    let tmp_b1 = make_temp_fastq(fq_b1);
    let tmp_b2 = make_temp_fastq(fq_b2);

    let groups = vec![
        FastqSampleGroup {
            fastq_paths: vec![tmp_a1.path().to_path_buf(), tmp_a2.path().to_path_buf()],
            label: "ctrl".to_string(),
        },
        FastqSampleGroup {
            fastq_paths: vec![tmp_b1.path().to_path_buf(), tmp_b2.path().to_path_buf()],
            label: "trt".to_string(),
        },
    ];

    let result = count_grouped_fastq(&groups, &index, 0, 20, false).expect("count 2x2");

    assert_eq!(result.counts.n_samples(), 2);
    assert_eq!(result.counts.n_guides(), 3);

    // Group 0: rep A1 + rep A2 -> g0 = 2.0
    assert_eq!(result.counts.sample_slice(0), &[2.0, 0.0, 0.0]);
    // Group 1: rep B1 + rep B2 -> g1 = 1.0, g2 = 1.0
    assert_eq!(result.counts.sample_slice(1), &[0.0, 1.0, 1.0]);

    // Global summary order: group-then-file order
    assert_eq!(result.summaries.len(), 4);
    assert_eq!(result.summaries[0].path, tmp_a1.path());
    assert_eq!(result.summaries[0].label, "ctrl");
    assert_eq!(result.summaries[1].path, tmp_a2.path());
    assert_eq!(result.summaries[1].label, "ctrl");
    assert_eq!(result.summaries[2].path, tmp_b1.path());
    assert_eq!(result.summaries[2].label, "trt");
    assert_eq!(result.summaries[3].path, tmp_b2.path());
    assert_eq!(result.summaries[3].label, "trt");
}

#[test]
fn test_grouped_empty_fastq_behavior() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // 1. Standalone empty file in a group
    let empty_tmp = make_temp_fastq("");
    let group_empty = FastqSampleGroup {
        fastq_paths: vec![empty_tmp.path().to_path_buf()],
        label: "empty_sample".to_string(),
    };

    let res_empty = count_grouped_fastq(&[group_empty], &index, 0, 20, false).expect("count empty");
    assert_eq!(res_empty.counts.sample_slice(0), &[0.0, 0.0, 0.0]);
    assert_eq!(res_empty.summaries.len(), 1);
    assert_eq!(res_empty.summaries[0].summary.total_reads, 0);
    assert_eq!(res_empty.summaries[0].summary.mapped_reads, 0);
    assert_eq!(res_empty.summaries[0].summary.mapped_percentage(), 0.0);
    assert_eq!(res_empty.summaries[0].zerocounts, 3); // all 3 library guides have 0 counts

    // 2. Empty replicate alongside a non-empty replicate in the same group
    let fq_nonzero = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let nonzero_tmp = make_temp_fastq(fq_nonzero);

    let group_mixed = FastqSampleGroup {
        fastq_paths: vec![
            nonzero_tmp.path().to_path_buf(),
            empty_tmp.path().to_path_buf(),
        ],
        label: "mixed_sample".to_string(),
    };

    let res_mixed = count_grouped_fastq(&[group_mixed], &index, 0, 20, false).expect("count mixed");

    // Empty replicate must NOT erase or reset counts from the non-empty replicate
    assert_eq!(res_mixed.counts.sample_slice(0), &[1.0, 0.0, 0.0]);
    assert_eq!(res_mixed.summaries.len(), 2);

    assert_eq!(res_mixed.summaries[0].path, nonzero_tmp.path());
    assert_eq!(res_mixed.summaries[0].summary.total_reads, 1);
    assert_eq!(res_mixed.summaries[0].summary.mapped_reads, 1);
    assert_eq!(res_mixed.summaries[0].zerocounts, 2);

    assert_eq!(res_mixed.summaries[1].path, empty_tmp.path());
    assert_eq!(res_mixed.summaries[1].summary.total_reads, 0);
    assert_eq!(res_mixed.summaries[1].summary.mapped_reads, 0);
    assert_eq!(res_mixed.summaries[1].summary.mapped_percentage(), 0.0);
    assert_eq!(res_mixed.summaries[1].zerocounts, 3);
}

#[test]
fn test_grouped_gzip_and_uncompressed_mixed() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    // Rep 1: uncompressed file hitting g1
    let fq_plain = "@r1\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let plain_tmp = make_temp_fastq(fq_plain);

    // Rep 2: gzipped file hitting g0
    let gz_tmp = make_temp_fastq_bytes(GZ_FASTQ_BYTES, ".fastq.gz");

    let group = FastqSampleGroup {
        fastq_paths: vec![plain_tmp.path().to_path_buf(), gz_tmp.path().to_path_buf()],
        label: "mixed_compression".to_string(),
    };

    let result = count_grouped_fastq(&[group], &index, 0, 20, false)
        .expect("count mixed compression in one group");

    // Both files processed sequentially in supplied order and aggregated
    assert_eq!(result.counts.sample_slice(0), &[1.0, 1.0, 0.0]);
    assert_eq!(result.summaries.len(), 2);

    assert_eq!(result.summaries[0].path, plain_tmp.path());
    assert_eq!(result.summaries[0].summary.exact_matches, 1);
    assert_eq!(result.summaries[0].zerocounts, 2);

    assert_eq!(result.summaries[1].path, gz_tmp.path());
    assert_eq!(result.summaries[1].summary.exact_matches, 1);
    assert_eq!(result.summaries[1].zerocounts, 2);
}

#[test]
fn test_grouped_repeated_physical_path() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    // Repeated identical physical path in the same sample group
    let group = FastqSampleGroup {
        fastq_paths: vec![tmp.path().to_path_buf(), tmp.path().to_path_buf()],
        label: "repeat_group".to_string(),
    };

    let result = count_grouped_fastq(&[group], &index, 0, 20, false).expect("count repeated path");

    // Counts are added twice
    assert_eq!(result.counts.sample_slice(0), &[2.0, 0.0, 0.0]);

    // Two distinct summary records are generated
    assert_eq!(result.summaries.len(), 2);
    assert_eq!(result.summaries[0].path, tmp.path());
    assert_eq!(result.summaries[0].label, "repeat_group");
    assert_eq!(result.summaries[0].summary.total_reads, 1);
    assert_eq!(result.summaries[0].zerocounts, 2);

    assert_eq!(result.summaries[1].path, tmp.path());
    assert_eq!(result.summaries[1].label, "repeat_group");
    assert_eq!(result.summaries[1].summary.total_reads, 1);
    assert_eq!(result.summaries[1].zerocounts, 2);
}

#[test]
fn test_grouped_read_classifications_and_invariants() {
    // 3 guides:
    // g0: ACGTACGTACGTACGTACGA (ends in A)
    // g1: ACGTACGTACGTACGTACGG (ends in G)
    // g2: TGCATGCATGCATGCATGCA
    let tsv = "id\tsequence\tgene\n\
               g0\tACGTACGTACGTACGTACGA\tGENE_A\n\
               g1\tACGTACGTACGTACGTACGG\tGENE_B\n\
               g2\tTGCATGCATGCATGCATGCA\tGENE_C\n";
    let lib = make_library(tsv);
    let index = build_index(&lib, true).expect("build index with mismatch");

    // File 1 has:
    // - exact match to g0: ACGTACGTACGTACGTACGA
    // - 1-mismatch match to g2: TGCATGCATGCATGCATGCT (ends in T instead of A)
    // - unmapped read: GGGGGGGGGGGGGGGGGGGG
    let fq1 = "@r1_exact\nACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIII\n\
               @r2_mismatch\nTGCATGCATGCATGCATGCT\n+\nIIIIIIIIIIIIIIIIIIII\n\
               @r3_unmapped\nGGGGGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n";

    // File 2 has:
    // - exact match to g0: ACGTACGTACGTACGTACGA
    // - ambiguous read: ACGTACGTACGTACGTACGC (1 mismatch from g0 and g1)
    // - short read: ACGT
    let fq2 = "@r4_exact\nACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIII\n\
               @r5_ambiguous\nACGTACGTACGTACGTACGC\n+\nIIIIIIIIIIIIIIIIIIII\n\
               @r6_short\nACGT\n+\nIIII\n";

    let tmp1 = make_temp_fastq(fq1);
    let tmp2 = make_temp_fastq(fq2);

    // Case 1: one_mismatch = true
    let group = FastqSampleGroup {
        fastq_paths: vec![tmp1.path().to_path_buf(), tmp2.path().to_path_buf()],
        label: "classification_test".to_string(),
    };

    let result = count_grouped_fastq(std::slice::from_ref(&group), &index, 0, 20, true)
        .expect("count with mismatch enabled");

    // g0 = 1 (file 1 exact) + 1 (file 2 exact) = 2.0
    // g1 = 0.0 (ambiguous read does not contribute!)
    // g2 = 1.0 (file 1 1-mismatch)
    assert_eq!(result.counts.sample_slice(0), &[2.0, 0.0, 1.0]);

    // Summary 1
    let s1 = &result.summaries[0];
    assert_eq!(s1.summary.total_reads, 3);
    assert_eq!(s1.summary.mapped_reads, 2);
    assert_eq!(s1.summary.exact_matches, 1);
    assert_eq!(s1.summary.one_mismatch_matches, 1);
    assert_eq!(s1.summary.ambiguous_reads, 0);
    assert_eq!(s1.summary.unmapped_reads, 1);
    assert_eq!(s1.summary.short_reads, 0);
    assert_eq!(s1.zerocounts, 1); // only g1 is 0

    // Summary 2
    let s2 = &result.summaries[1];
    assert_eq!(s2.summary.total_reads, 3);
    assert_eq!(s2.summary.mapped_reads, 1);
    assert_eq!(s2.summary.exact_matches, 1);
    assert_eq!(s2.summary.one_mismatch_matches, 0);
    assert_eq!(s2.summary.ambiguous_reads, 1);
    assert_eq!(s2.summary.unmapped_reads, 0);
    assert_eq!(s2.summary.short_reads, 1);
    assert_eq!(s2.zerocounts, 2); // g1 and g2 are 0

    // Case 2: one_mismatch = false
    let result_no_mm = count_grouped_fastq(std::slice::from_ref(&group), &index, 0, 20, false)
        .expect("count with mismatch disabled");

    // 1-mismatch match is treated as unmapped when disabled
    assert_eq!(result_no_mm.counts.sample_slice(0), &[2.0, 0.0, 0.0]);
    assert_eq!(result_no_mm.summaries[0].summary.one_mismatch_matches, 0);
    assert_eq!(result_no_mm.summaries[0].summary.unmapped_reads, 2);
    assert_eq!(result_no_mm.summaries[0].summary.mapped_reads, 1);
    assert_eq!(result_no_mm.summaries[0].zerocounts, 2);
}

#[test]
fn test_grouped_failure_on_later_physical_file() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq_good = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp_good = make_temp_fastq(fq_good);
    let non_existent = PathBuf::from("/nonexistent/path/to/replicate2.fastq");

    let groups = vec![
        FastqSampleGroup {
            fastq_paths: vec![tmp_good.path().to_path_buf()],
            label: "group1".to_string(),
        },
        FastqSampleGroup {
            fastq_paths: vec![tmp_good.path().to_path_buf(), non_existent],
            label: "group2".to_string(),
        },
    ];

    let result = count_grouped_fastq(&groups, &index, 0, 20, false);
    assert!(result.is_err(), "later failing file must return error");

    match result.unwrap_err() {
        MageckError::Fastq(msg) => {
            assert!(
                msg.contains("Failed to open FASTQ file"),
                "error message should cite open failure: {msg}"
            );
        }
        other => panic!("expected MageckError::Fastq, got {other:?}"),
    }
}

#[test]
fn test_grouped_malformed_record_fails_immediately() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq_good = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let fq_bad = "@truncated_header";
    let tmp_good = make_temp_fastq(fq_good);
    let tmp_bad = make_temp_fastq(fq_bad);

    let group = FastqSampleGroup {
        fastq_paths: vec![tmp_good.path().to_path_buf(), tmp_bad.path().to_path_buf()],
        label: "malformed_test".to_string(),
    };

    let result = count_grouped_fastq(&[group], &index, 0, 20, false);
    assert!(result.is_err(), "malformed record must fail");
    assert!(matches!(result.unwrap_err(), MageckError::Fastq(_)));
}

#[test]
fn test_grouped_invalid_window_configuration() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let fq = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let tmp = make_temp_fastq(fq);

    let group = FastqSampleGroup {
        fastq_paths: vec![tmp.path().to_path_buf()],
        label: "sample1".to_string(),
    };

    // sgrna_len == 0
    let res_zero = count_grouped_fastq(std::slice::from_ref(&group), &index, 0, 0, false);
    assert!(matches!(res_zero.unwrap_err(), MageckError::Fastq(_)));

    // Window overflow
    let res_overflow =
        count_grouped_fastq(std::slice::from_ref(&group), &index, usize::MAX, 20, false);
    assert!(matches!(res_overflow.unwrap_err(), MageckError::Fastq(_)));
}

#[test]
fn test_grouped_empty_groups_list() {
    let lib = make_library(STD_LIB_TSV);
    let index = build_index(&lib, false).expect("build index");

    let result = count_grouped_fastq(&[], &index, 0, 20, false).expect("empty groups list");
    assert_eq!(result.counts.n_samples(), 0);
    assert_eq!(result.counts.n_guides(), 3);
    assert!(result.summaries.is_empty());
}
