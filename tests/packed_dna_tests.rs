//! Exhaustive correctness and property tests for PackedDna, PackedGuideIndex,
//! and ByteGuideIndex — Phase 2B.
//!
//! Tests follow the requirements from spec §16 (packed_dna_tests.rs).

use mageckrs::data::packed_dna::{
    ByteGuideIndex, EncodingError, MatchResult, PackedDna, PackedGuideIndex,
};

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// Encode a string literal; panics if encoding fails (test helper).
fn enc(s: &str) -> PackedDna {
    PackedDna::encode(s.as_bytes()).expect("valid sequence")
}

// ─── Empty sequence ───────────────────────────────────────────────────────────

#[test]
fn empty_sequence_is_valid() {
    let p = PackedDna::encode(b"").expect("empty sequence should be valid");
    assert_eq!(p.len(), 0);
    assert!(p.is_empty());
    assert_eq!(p.decode(), b"");
}

#[test]
fn base_at_accessor() {
    let p = enc("ACGT");
    assert_eq!(p.base_at(0), b'A');
    assert_eq!(p.base_at(1), b'C');
    assert_eq!(p.base_at(2), b'G');
    assert_eq!(p.base_at(3), b'T');
}

// ─── Multi-base encode/decode round-trip ─────────────────────────────────────

#[test]
fn round_trip_short() {
    let seq = b"GATTACA";
    let p = PackedDna::encode(seq).unwrap();
    assert_eq!(p.decode(), seq.to_vec());
}

#[test]
fn round_trip_len_20() {
    let seq = b"ACGTACGTACGTACGTACGT";
    assert_eq!(seq.len(), 20);
    let p = PackedDna::encode(seq).unwrap();
    assert_eq!(p.decode(), seq.to_vec());
}

#[test]
fn round_trip_len_23() {
    let seq = b"ACGTACGTACGTACGTACGTACG";
    assert_eq!(seq.len(), 23);
    let p = PackedDna::encode(seq).unwrap();
    assert_eq!(p.decode(), seq.to_vec());
}

#[test]
fn round_trip_len_24() {
    let seq = b"ACGTACGTACGTACGTACGTACGT";
    assert_eq!(seq.len(), 24);
    let p = PackedDna::encode(seq).unwrap();
    assert_eq!(p.decode(), seq.to_vec());
}

// ─── Maximum length of 32 nt ─────────────────────────────────────────────────

#[test]
fn round_trip_len_32() {
    // Use all four bases cycling to test every position up to 32.
    let seq: Vec<u8> = (0..32).map(|i| b"ACGT"[i % 4]).collect();
    assert_eq!(seq.len(), 32);
    let p = PackedDna::encode(&seq).unwrap();
    assert_eq!(p.len(), 32);
    assert_eq!(p.decode(), seq);
}

// ─── Rejection of length 33 ──────────────────────────────────────────────────

#[test]
fn reject_length_33() {
    let seq: Vec<u8> = vec![b'A'; 33];
    let result = PackedDna::encode(&seq);
    assert_eq!(result, Err(EncodingError::TooLong(33)));
}

#[test]
fn reject_length_100() {
    let seq: Vec<u8> = vec![b'G'; 100];
    let result = PackedDna::encode(&seq);
    assert_eq!(result, Err(EncodingError::TooLong(100)));
}

// ─── Invalid nucleotide rejection ─────────────────────────────────────────────

#[test]
fn reject_invalid_base_n() {
    let result = PackedDna::encode(b"ACGN");
    assert!(matches!(result, Err(EncodingError::InvalidBase(b'N'))));
}

#[test]
fn reject_invalid_base_u() {
    let result = PackedDna::encode(b"U");
    assert!(matches!(result, Err(EncodingError::InvalidBase(b'U'))));
}

#[test]
fn reject_invalid_base_digit() {
    let result = PackedDna::encode(b"1");
    assert!(matches!(result, Err(EncodingError::InvalidBase(b'1'))));
}

#[test]
fn reject_invalid_base_space() {
    let result = PackedDna::encode(b"A C");
    assert!(matches!(result, Err(EncodingError::InvalidBase(b' '))));
}

// ─── Case behaviour ──────────────────────────────────────────────────────────

#[test]
fn lowercase_normalised_to_uppercase() {
    let lower = PackedDna::encode(b"acgt").unwrap();
    let upper = PackedDna::encode(b"ACGT").unwrap();
    assert_eq!(
        lower, upper,
        "lowercase and uppercase must produce identical PackedDna"
    );
    assert_eq!(lower.decode(), b"ACGT");
}

#[test]
fn mixed_case_normalised() {
    let mixed = PackedDna::encode(b"AcGt").unwrap();
    let upper = PackedDna::encode(b"ACGT").unwrap();
    assert_eq!(mixed, upper);
}

// ─── Length participates in equality and hashing ──────────────────────────────

#[test]
fn a_neq_aa() {
    let a = enc("A");
    let aa = enc("AA");
    assert_ne!(a, aa, "\"A\" must not equal \"AA\"");
}

#[test]
fn length_in_hash_set() {
    use std::collections::HashSet;
    let mut set = HashSet::new();
    set.insert(enc("A"));
    set.insert(enc("AA"));
    assert_eq!(set.len(), 2, "A and AA must hash to different buckets");
}

// ─── replace_base ─────────────────────────────────────────────────────────────

#[test]
fn replace_base_first_position() {
    let p = enc("ACGT");
    let q = p.replace_base(0, b'T');
    assert_eq!(q.decode(), b"TCGT");
    assert_eq!(q.len(), p.len());
}

#[test]
fn replace_base_last_position() {
    let p = enc("ACGT");
    let q = p.replace_base(3, b'A');
    assert_eq!(q.decode(), b"ACGA");
    assert_eq!(q.len(), p.len());
}

#[test]
fn replace_base_internal_position_1() {
    let p = enc("ACGT");
    let q = p.replace_base(1, b'G');
    assert_eq!(q.decode(), b"AGGT");
}

#[test]
fn replace_base_internal_position_2() {
    let p = enc("ACGT");
    let q = p.replace_base(2, b'C');
    assert_eq!(q.decode(), b"ACCT");
}

#[test]
fn replace_base_preserves_all_non_target_bases() {
    let p = enc("GATTACA");
    // Replace pos 3 (T→C).
    let q = p.replace_base(3, b'C');
    assert_eq!(q.decode(), b"GATCACA");
    // Every other position must be unchanged.
    for i in 0..7 {
        if i != 3 {
            assert_eq!(
                q.base_at(i),
                p.base_at(i),
                "position {i} should be unchanged"
            );
        }
    }
}

#[test]
fn replace_base_preserves_length() {
    let p = enc("GATTACA");
    let q = p.replace_base(0, b'C');
    assert_eq!(q.len(), p.len());
}

#[test]
fn replace_base_lowercase_accepted() {
    let p = enc("ACGT");
    let q = p.replace_base(0, b't');
    assert_eq!(q.decode(), b"TCGT");
}

#[test]
#[should_panic(expected = "out of range")]
fn replace_base_out_of_range_panics() {
    let p = enc("ACGT");
    p.replace_base(4, b'A'); // pos 4 is out of range for len=4
}

#[test]
#[should_panic(expected = "invalid base")]
fn replace_base_invalid_base_panics() {
    let p = enc("ACGT");
    p.replace_base(0, b'N');
}

#[test]
fn replace_base_all_positions_acgt() {
    // Verify every position × every base combination round-trips correctly.
    let orig = enc("TTTT");
    for pos in 0..4usize {
        for &base in b"ACGT" {
            let q = orig.replace_base(pos, base);
            assert_eq!(q.base_at(pos), base.to_ascii_uppercase());
            assert_eq!(q.len(), 4);
            // All other positions remain T.
            for other in 0..4usize {
                if other != pos {
                    assert_eq!(q.base_at(other), b'T');
                }
            }
        }
    }
}

// ─── PackedGuideIndex — exact-match lookup ────────────────────────────────────

fn build_packed_exact(seqs: &[&str]) -> PackedGuideIndex {
    let guides: Vec<PackedDna> = seqs.iter().map(|s| enc(s)).collect();
    PackedGuideIndex::build(&guides, false).expect("no duplicates")
}

#[test]
fn packed_exact_hit() {
    let idx = build_packed_exact(&["ACGT", "TTTT", "GGGG"]);
    assert_eq!(idx.lookup(&enc("ACGT")), MatchResult::Exact(0));
    assert_eq!(idx.lookup(&enc("TTTT")), MatchResult::Exact(1));
    assert_eq!(idx.lookup(&enc("GGGG")), MatchResult::Exact(2));
}

#[test]
fn packed_exact_miss() {
    let idx = build_packed_exact(&["ACGT"]);
    assert_eq!(idx.lookup(&enc("AAGT")), MatchResult::NoMatch);
}

#[test]
fn packed_exact_different_lengths_distinct() {
    let idx = build_packed_exact(&["A", "AA", "AAA"]);
    assert_eq!(idx.lookup(&enc("A")), MatchResult::Exact(0));
    assert_eq!(idx.lookup(&enc("AA")), MatchResult::Exact(1));
    assert_eq!(idx.lookup(&enc("AAA")), MatchResult::Exact(2));
}

// ─── PackedGuideIndex — duplicate detection ───────────────────────────────────

#[test]
fn packed_duplicate_guide_returns_err() {
    let guides = vec![enc("ACGT"), enc("TTTT"), enc("ACGT")];
    let result = PackedGuideIndex::build(&guides, false);
    assert!(result.is_err(), "duplicate guide must cause an error");
}

// ─── PackedGuideIndex — one-mismatch lookup ───────────────────────────────────

fn build_packed_1mm(seqs: &[&str]) -> PackedGuideIndex {
    let guides: Vec<PackedDna> = seqs.iter().map(|s| enc(s)).collect();
    PackedGuideIndex::build(&guides, true).expect("no duplicates")
}

#[test]
fn packed_one_mismatch_unique_match() {
    // "ACGT" is in the index; query "ACGG" differs by 1 at pos 3.
    let idx = build_packed_1mm(&["ACGT", "CCCC"]);
    let q = enc("ACGG"); // 1 mismatch from ACGT, no match to CCCC
    assert_eq!(idx.lookup(&q), MatchResult::OneMismatch(0));
}

#[test]
fn packed_exact_takes_precedence_over_one_mismatch() {
    // "ACGG" is also in the index directly.
    let idx = build_packed_1mm(&["ACGT", "ACGG"]);
    let q = enc("ACGG");
    assert_eq!(idx.lookup(&q), MatchResult::Exact(1));
}

#[test]
fn packed_one_mismatch_ambiguous() {
    // "ACGT" and "ACGC" both match query "ACGA" at 1 mismatch (pos 3).
    let idx = build_packed_1mm(&["ACGT", "ACGC"]);
    let q = enc("ACGA"); // 1 mismatch from both guides
    assert_eq!(idx.lookup(&q), MatchResult::Ambiguous);
}

#[test]
fn packed_no_match_in_one_mismatch_mode() {
    let idx = build_packed_1mm(&["ACGT"]);
    // "TTTT" requires 4 mismatches.
    assert_eq!(idx.lookup(&enc("TTTT")), MatchResult::NoMatch);
}

// ─── ByteGuideIndex — long-guide fallback ────────────────────────────────────

fn bytes(s: &str) -> Vec<u8> {
    s.bytes().map(|b| b.to_ascii_uppercase()).collect()
}

fn build_byte_exact(seqs: &[&str]) -> ByteGuideIndex {
    let guides: Vec<Vec<u8>> = seqs.iter().map(|s| bytes(s)).collect();
    ByteGuideIndex::build(&guides, false).expect("no duplicates")
}

fn build_byte_1mm(seqs: &[&str]) -> ByteGuideIndex {
    let guides: Vec<Vec<u8>> = seqs.iter().map(|s| bytes(s)).collect();
    ByteGuideIndex::build(&guides, true).expect("no duplicates")
}

/// Guides longer than 32 nt must work correctly in the byte index.
#[test]
fn byte_index_long_guide_accepted() {
    let long: String = "ACGT".repeat(9); // 36 nt — exceeds PackedDna::MAX_LEN
    assert!(long.len() > 32);
    let idx = build_byte_exact(&[&long]);
    assert_eq!(idx.lookup(long.as_bytes()), MatchResult::Exact(0));
}

#[test]
fn byte_index_no_truncation_of_long_guide() {
    // A 33-nt guide and a 32-nt prefix must be distinct.
    let long33: String = "A".repeat(33);
    let long32: String = "A".repeat(32);
    let idx = build_byte_exact(&[&long33]);
    assert_eq!(idx.lookup(long33.as_bytes()), MatchResult::Exact(0));
    assert_eq!(
        idx.lookup(long32.as_bytes()),
        MatchResult::NoMatch,
        "truncated query must not match"
    );
}

#[test]
fn byte_index_exact_miss() {
    let idx = build_byte_exact(&["ACGT"]);
    assert_eq!(idx.lookup(b"AAGT"), MatchResult::NoMatch);
}

#[test]
fn byte_index_duplicate_returns_err() {
    let guides = vec![bytes("ACGT"), bytes("TTTT"), bytes("ACGT")];
    let result = ByteGuideIndex::build(&guides, false);
    assert!(result.is_err());
}

#[test]
fn byte_index_one_mismatch_unique() {
    let idx = build_byte_1mm(&["ACGT", "CCCC"]);
    let q = bytes("ACGG"); // 1 mismatch from ACGT
    assert_eq!(idx.lookup(&q), MatchResult::OneMismatch(0));
}

#[test]
fn byte_index_one_mismatch_ambiguous() {
    let idx = build_byte_1mm(&["ACGT", "ACGC"]);
    let q = bytes("ACGA");
    assert_eq!(idx.lookup(&q), MatchResult::Ambiguous);
}

#[test]
fn byte_index_exact_takes_precedence() {
    let idx = build_byte_1mm(&["ACGT", "ACGG"]);
    let q = bytes("ACGG");
    assert_eq!(idx.lookup(&q), MatchResult::Exact(1));
}

// ─── Differential: PackedGuideIndex vs ByteGuideIndex ────────────────────────

/// For short guides both indexes must return identical results.
///
/// This drives the §16 requirement: "PackedGuideIndex and ByteGuideIndex
/// produce identical match results for 10,000 randomized libraries."
///
/// We use a deterministic LCG (no external deps) to generate 10,000 queries
/// against a fixed 16-guide library in 1-mismatch mode.
#[test]
fn differential_packed_vs_byte_10000_queries() {
    // Fixed 16-guide library — each 8 nt, all distinct.
    let library_seqs: Vec<&str> = vec![
        "AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT", "ACGTACGT", "TGCATGCA", "CATGCATG",
        "GTACGTAC", "AAAACCCC", "TTTTGGGG", "ACACACAC", "TGTGTGTG", "GCGCGCGC", "AGTCAGTC",
        "TCGATCGA", "CAGTCAGT",
    ];
    let packed_guides: Vec<PackedDna> = library_seqs.iter().map(|s| enc(s)).collect();
    let byte_guides: Vec<Vec<u8>> = library_seqs.iter().map(|s| bytes(s)).collect();

    let packed_idx = PackedGuideIndex::build(&packed_guides, true).unwrap();
    let byte_idx = ByteGuideIndex::build(&byte_guides, true).unwrap();

    // Deterministic LCG PRNG (same as bench_utils).
    let mut state: u64 = 0xdeadbeef_cafebabe;
    let alphabet = b"ACGT";

    for _ in 0..10_000 {
        // Generate a random 8-nt query.
        let query_bytes: Vec<u8> = (0..8)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                alphabet[((state >> 33) as usize) & 3]
            })
            .collect();

        // Packed lookup.
        let packed_result = match PackedDna::encode(&query_bytes) {
            Ok(p) => packed_idx.lookup(&p),
            Err(_) => panic!("generated invalid sequence"),
        };
        // Byte lookup.
        let byte_result = byte_idx.lookup(&query_bytes);

        // The two indexes must agree.
        assert_eq!(
            packed_result,
            byte_result,
            "packed and byte indexes disagree for query {:?}",
            std::str::from_utf8(&query_bytes).unwrap_or("<invalid>")
        );
    }
}
