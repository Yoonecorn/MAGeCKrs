//! Tests for UmiDeduper — Phase 2B.
//!
//! Covers packed mode (≤ 32 nt) and byte mode (> 32 nt), boundary lengths,
//! duplicate detection, group independence, validation, and error paths.

use mageckrs::data::umi::{UmiDeduper, UmiError};

// ─── Storage-mode selection ───────────────────────────────────────────────────

#[test]
fn packed_mode_for_umi_len_32() {
    let d = UmiDeduper::new(1, 32);
    assert!(d.is_packed(), "umi_len=32 must use packed mode");
}

#[test]
fn packed_mode_for_umi_len_below_32() {
    let d = UmiDeduper::new(1, 12);
    assert!(d.is_packed(), "umi_len=12 must use packed mode");
}

#[test]
fn byte_mode_for_umi_len_33() {
    let d = UmiDeduper::new(1, 33);
    assert!(!d.is_packed(), "umi_len=33 must use byte mode");
}

#[test]
fn byte_mode_for_umi_len_above_32() {
    let d = UmiDeduper::new(1, 50);
    assert!(!d.is_packed(), "umi_len=50 must use byte mode");
}

/// Boundary: umi_len=32 is still packed; umi_len=33 is byte.
#[test]
fn boundary_32_is_packed_33_is_byte() {
    assert!(UmiDeduper::new(1, 32).is_packed());
    assert!(!UmiDeduper::new(1, 33).is_packed());
}

// ─── First observation is accepted ───────────────────────────────────────────

#[test]
fn packed_first_observation_returns_true() {
    let mut d = UmiDeduper::new(1, 8);
    assert_eq!(d.observe(0, b"ACGTACGT"), Ok(true));
}

#[test]
fn byte_first_observation_returns_true() {
    let mut d = UmiDeduper::new(1, 33);
    let umi: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGTA".to_vec(); // 33 nt
    assert_eq!(umi.len(), 33);
    assert_eq!(d.observe(0, &umi), Ok(true));
}

// ─── Second identical observation is a duplicate ──────────────────────────────

#[test]
fn packed_second_observation_returns_false() {
    let mut d = UmiDeduper::new(1, 8);
    assert_eq!(d.observe(0, b"ACGTACGT"), Ok(true));
    assert_eq!(
        d.observe(0, b"ACGTACGT"),
        Ok(false),
        "duplicate must return false"
    );
}

#[test]
fn byte_second_observation_returns_false() {
    let mut d = UmiDeduper::new(1, 33);
    let umi: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGTA".to_vec();
    d.observe(0, &umi).unwrap();
    assert_eq!(d.observe(0, &umi), Ok(false));
}

// ─── Same UMI in different groups remains independent ─────────────────────────

#[test]
fn same_umi_different_groups_independent_packed() {
    let mut d = UmiDeduper::new(3, 8);
    let umi = b"ACGTACGT";
    assert_eq!(d.observe(0, umi), Ok(true));
    assert_eq!(
        d.observe(1, umi),
        Ok(true),
        "same UMI in group 1 is a new observation"
    );
    assert_eq!(
        d.observe(2, umi),
        Ok(true),
        "same UMI in group 2 is a new observation"
    );
    // Second time in each group → duplicate.
    assert_eq!(d.observe(0, umi), Ok(false));
    assert_eq!(d.observe(1, umi), Ok(false));
    assert_eq!(d.observe(2, umi), Ok(false));
}

#[test]
fn same_umi_different_groups_independent_byte() {
    let mut d = UmiDeduper::new(2, 33);
    let umi: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGTA".to_vec();
    assert_eq!(d.observe(0, &umi), Ok(true));
    assert_eq!(d.observe(1, &umi), Ok(true));
    assert_eq!(d.observe(0, &umi), Ok(false));
}

// ─── Invalid UMI bases ────────────────────────────────────────────────────────

#[test]
fn packed_invalid_base_returns_error() {
    let mut d = UmiDeduper::new(1, 4);
    // umi_len=4, but "ACGN" contains 'N'
    let result = d.observe(0, b"ACGN");
    assert!(
        matches!(result, Err(UmiError::InvalidBase(_))),
        "expected InvalidBase, got {result:?}"
    );
}

#[test]
fn byte_invalid_base_returns_error() {
    let mut d = UmiDeduper::new(1, 33);
    let mut umi: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGTA".to_vec();
    umi[0] = b'N'; // corrupt first base
    let result = d.observe(0, &umi);
    assert!(matches!(result, Err(UmiError::InvalidBase(_))));
}

// ─── Configured-length mismatch ───────────────────────────────────────────────

#[test]
fn packed_length_mismatch_returns_error() {
    let mut d = UmiDeduper::new(1, 8);
    let result = d.observe(0, b"ACGT"); // 4 nt, configured 8
    assert_eq!(
        result,
        Err(UmiError::LengthMismatch {
            expected: 8,
            got: 4
        })
    );
}

#[test]
fn byte_length_mismatch_returns_error() {
    let mut d = UmiDeduper::new(1, 33);
    let short: Vec<u8> = b"ACGTACGT".to_vec(); // 8 nt, configured 33
    let result = d.observe(0, &short);
    assert_eq!(
        result,
        Err(UmiError::LengthMismatch {
            expected: 33,
            got: 8
        })
    );
}

// ─── Empty UMI behaviour ──────────────────────────────────────────────────────

/// Empty UMIs are accepted when umi_len == 0.
#[test]
fn empty_umi_accepted_when_len_zero() {
    let mut d = UmiDeduper::new(1, 0);
    assert_eq!(d.observe(0, b""), Ok(true));
    assert_eq!(d.observe(0, b""), Ok(false)); // second -> duplicate
}

#[test]
fn non_empty_umi_rejected_when_len_zero() {
    let mut d = UmiDeduper::new(1, 0);
    let result = d.observe(0, b"ACGT");
    assert_eq!(
        result,
        Err(UmiError::LengthMismatch {
            expected: 0,
            got: 4
        })
    );
}

/// Empty UMIs return UmiError::Empty when a fixed length is configured.
#[test]
fn empty_umi_rejected_when_fixed_length_configured() {
    let mut d = UmiDeduper::new(1, 8);
    let result = d.observe(0, b"");
    assert_eq!(result, Err(UmiError::Empty));
}

// ─── Long UMI is not truncated ────────────────────────────────────────────────

/// A 33-nt UMI and the same 32-nt prefix must be independently deduplicated.
/// Uses separate dedupers because the storage mode is fixed at construction.
#[test]
fn long_umi_not_truncated() {
    let full: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGTA".to_vec(); // 33 nt
    let prefix: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(); // 32 nt
    assert_eq!(full.len(), 33);
    assert_eq!(prefix.len(), 32);

    // Byte-mode deduper for 33-nt UMIs — must reject 32-nt (not silently truncate).
    let mut d33 = UmiDeduper::new(1, 33);
    assert!(!d33.is_packed());
    assert_eq!(d33.observe(0, &full), Ok(true));
    let result = d33.observe(0, &prefix);
    assert_eq!(
        result,
        Err(UmiError::LengthMismatch {
            expected: 33,
            got: 32
        }),
        "32-nt query must be rejected, not silently truncated"
    );

    // Packed-mode deduper for 32-nt UMIs — must reject 33-nt (not silently truncate).
    let mut d32 = UmiDeduper::new(1, 32);
    assert!(d32.is_packed());
    assert_eq!(d32.observe(0, &prefix), Ok(true));
    let result2 = d32.observe(0, &full);
    assert_eq!(
        result2,
        Err(UmiError::LengthMismatch {
            expected: 32,
            got: 33
        }),
        "33-nt query must be rejected in 32-nt deduper"
    );
}

// ─── Group out of range ───────────────────────────────────────────────────────

#[test]
fn group_out_of_range_packed() {
    let mut d = UmiDeduper::new(2, 4);
    let result = d.observe(5, b"ACGT");
    assert_eq!(
        result,
        Err(UmiError::GroupOutOfRange {
            group: 5,
            n_groups: 2
        })
    );
}

#[test]
fn group_out_of_range_byte() {
    let mut d = UmiDeduper::new(1, 33);
    let umi: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGTA".to_vec();
    let result = d.observe(1, &umi); // only group 0 exists
    assert_eq!(
        result,
        Err(UmiError::GroupOutOfRange {
            group: 1,
            n_groups: 1
        })
    );
}

// ─── Multiple groups remain independent ───────────────────────────────────────

#[test]
fn multiple_groups_independent_packed() {
    let mut d = UmiDeduper::new(5, 4);
    // Each group gets its own UMI.
    for g in 0..5usize {
        let umi: Vec<u8> = vec![b"ACGT"[g % 4]; 4];
        assert_eq!(d.observe(g, &umi), Ok(true));
    }
    // Re-observe in each group → all duplicates.
    for g in 0..5usize {
        let umi: Vec<u8> = vec![b"ACGT"[g % 4]; 4];
        assert_eq!(d.observe(g, &umi), Ok(false));
    }
    // Observe in a different group → still new.
    assert_eq!(d.observe(0, b"TTTT"), Ok(true));
    assert_eq!(d.observe(1, b"AAAA"), Ok(true));
}

#[test]
fn multiple_groups_independent_byte() {
    // Use umi_len=33 to force byte mode (> 32 nt threshold).
    let mut d = UmiDeduper::new(3, 33);
    assert!(!d.is_packed());
    let umis: Vec<Vec<u8>> = vec![
        b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(), // 33 A
        b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_vec(), // 33 C
        b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".to_vec(), // 33 G
    ];
    // First observations.
    for (g, umi) in umis.iter().enumerate() {
        assert_eq!(d.observe(g, umi), Ok(true));
    }
    // Second observations — duplicates.
    for (g, umi) in umis.iter().enumerate() {
        assert_eq!(d.observe(g, umi), Ok(false));
    }
}

// ─── Case normalisation ───────────────────────────────────────────────────────

#[test]
fn lowercase_umi_treated_same_as_uppercase_packed() {
    let mut d = UmiDeduper::new(1, 4);
    assert_eq!(d.observe(0, b"acgt"), Ok(true));
    assert_eq!(
        d.observe(0, b"ACGT"),
        Ok(false),
        "lowercase and uppercase are the same UMI"
    );
}

#[test]
fn lowercase_umi_treated_same_as_uppercase_byte() {
    let mut d = UmiDeduper::new(1, 33);
    let lower: Vec<u8> = b"acgtacgtacgtacgtacgtacgtacgtacgta".to_vec();
    let upper: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGTA".to_vec();
    assert_eq!(d.observe(0, &lower), Ok(true));
    assert_eq!(d.observe(0, &upper), Ok(false));
}

// ─── Accessor sanity ─────────────────────────────────────────────────────────

#[test]
fn n_groups_accessor() {
    let d = UmiDeduper::new(7, 8);
    assert_eq!(d.n_groups(), 7);
}

#[test]
fn umi_len_accessor() {
    let d = UmiDeduper::new(1, 20);
    assert_eq!(d.umi_len(), 20);
}

#[test]
fn campaign_key_cap_validates_then_checks_duplicate_without_insertion_on_failure() {
    use mageckrs::data::umi::{UmiDeduper, UmiError};
    for length in [2, 33] {
        let a = vec![b'a'; length];
        let c = vec![b'C'; length];
        let mut d = UmiDeduper::new(2, length);
        assert_eq!(d.observe_capped(0, &a, 1), Ok(true));
        assert_eq!(d.observe_capped(0, &vec![b'A'; length], 1), Ok(false));
        assert_eq!(
            d.observe_capped(0, &c, 1),
            Err(UmiError::KeyCapExceeded { limit: 1 })
        );
        assert_eq!(
            d.observe_capped(1, &a, 1),
            Err(UmiError::KeyCapExceeded { limit: 1 })
        );
        assert!(matches!(
            d.observe_capped(0, &vec![b'N'; length], 1),
            Err(UmiError::InvalidBase(_))
        ));
        assert_eq!(d.observe_capped(0, &c, 2), Ok(true));
        assert_eq!(d.observe_capped(0, &c, 2), Ok(false));
    }
}
