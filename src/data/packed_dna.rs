//! Packed-DNA guide key and guide-index types.
//!
//! [`PackedDna`] stores a nucleotide sequence of up to 32 bases in a `u64`,
//! using the encoding A=00 C=01 G=10 T=11 with the first base in the two
//! least-significant bits.  The stored length is always part of the value so
//! that "A" (len=1) and "AA" (len=2) are never equal.
//!
//! [`PackedGuideIndex`] and [`ByteGuideIndex`] provide exact-match and
//! one-mismatch lookup over a library of sgRNA guides.  The packed variant is
//! used when every guide in the library is ≤ 32 nt; the byte variant is the
//! fallback for libraries that contain any guide longer than 32 nt.

use std::collections::HashMap;

// ─── Encoding error ───────────────────────────────────────────────────────────

/// Errors that can occur when encoding a nucleotide sequence into [`PackedDna`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EncodingError {
    /// The input sequence exceeds the 32-nucleotide maximum.
    TooLong(usize),
    /// A byte in the input is not A, C, G, T (after case normalisation).
    InvalidBase(u8),
}

impl std::fmt::Display for EncodingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EncodingError::TooLong(n) => {
                write!(f, "sequence length {n} exceeds the 32-nt PackedDna maximum")
            }
            EncodingError::InvalidBase(b) => {
                write!(f, "invalid nucleotide byte 0x{b:02x} ({b})")
            }
        }
    }
}

impl std::error::Error for EncodingError {}

// ─── PackedDna ────────────────────────────────────────────────────────────────

/// A nucleotide sequence of 0–32 bases packed into a `u64`.
///
/// Encoding: A=00 C=01 G=10 T=11.
/// The first base occupies bits 1:0; the second base occupies bits 3:2; and so on.
/// The `len` field is stored separately so that distinct-length sequences
/// (e.g. "A" and "AA") always compare unequal.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct PackedDna {
    bits: u64,
    len: u8,
}

/// Map an uppercase ACGT byte to its two-bit code, or return an error.
#[inline]
fn encode_base(b: u8) -> Result<u64, EncodingError> {
    match b {
        b'A' => Ok(0b00),
        b'C' => Ok(0b01),
        b'G' => Ok(0b10),
        b'T' => Ok(0b11),
        _ => Err(EncodingError::InvalidBase(b)),
    }
}

/// Map a two-bit code back to an uppercase ACGT byte.
#[inline]
fn decode_base(code: u64) -> u8 {
    match code & 0b11 {
        0b00 => b'A',
        0b01 => b'C',
        0b10 => b'G',
        0b11 => b'T',
        _ => unreachable!(),
    }
}

impl PackedDna {
    /// The maximum sequence length that can be stored in a `PackedDna`.
    pub const MAX_LEN: usize = 32;

    /// Encode a nucleotide byte slice into a [`PackedDna`].
    ///
    /// - Both upper- and lower-case ACGT are accepted; all letters are
    ///   normalised to upper-case before encoding.
    /// - Returns [`EncodingError::TooLong`] for sequences longer than 32 nt.
    /// - Returns [`EncodingError::InvalidBase`] for any byte that is not
    ///   A, C, G, T (after case normalisation).
    /// - Empty sequences (len = 0) are valid.
    pub fn encode(seq: &[u8]) -> Result<Self, EncodingError> {
        let n = seq.len();
        if n > Self::MAX_LEN {
            return Err(EncodingError::TooLong(n));
        }
        let mut bits: u64 = 0;
        for (i, &raw) in seq.iter().enumerate() {
            let b = raw.to_ascii_uppercase();
            let code = encode_base(b)?;
            bits |= code << (2 * i);
        }
        let len = u8::try_from(n).expect("length ≤ 32 fits in u8");
        Ok(PackedDna { bits, len })
    }

    /// Decode back to a `Vec<u8>` of uppercase nucleotide bytes.
    ///
    /// Intended for tests and diagnostics only.
    pub fn decode(self) -> Vec<u8> {
        (0..self.len as usize)
            .map(|i| decode_base((self.bits >> (2 * i)) & 0b11))
            .collect()
    }

    /// The number of nucleotides stored.
    #[inline]
    pub fn len(&self) -> usize {
        self.len as usize
    }

    /// Returns `true` when the sequence has zero nucleotides.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Return the base at position `pos` (0-indexed) as an uppercase byte.
    ///
    /// # Panics
    /// Panics if `pos >= self.len()`.
    pub fn base_at(&self, pos: usize) -> u8 {
        assert!(
            pos < self.len as usize,
            "pos {pos} out of range for len {}",
            self.len
        );
        decode_base((self.bits >> (2 * pos)) & 0b11)
    }

    /// Return a new `PackedDna` with the base at `pos` replaced by `new_base`.
    ///
    /// - `new_base` must be an uppercase or lowercase A, C, G, or T.
    /// - The replacement is performed entirely on the packed bits — no
    ///   decode/re-encode of the full sequence occurs and no heap allocation
    ///   is made.
    /// - All other bases and the stored `len` are preserved exactly.
    ///
    /// # Panics
    /// Panics if `pos >= self.len()` or if `new_base` is not a valid ACGT
    /// byte.
    pub fn replace_base(&self, pos: usize, new_base: u8) -> PackedDna {
        assert!(
            pos < self.len as usize,
            "replace_base: pos {pos} out of range for len {}",
            self.len
        );
        let normalised = new_base.to_ascii_uppercase();
        let code = encode_base(normalised)
            .unwrap_or_else(|_| panic!("replace_base: invalid base byte 0x{new_base:02x}"));
        let shift = 2 * pos;
        let new_bits = (self.bits & !(0b11u64 << shift)) | (code << shift);
        PackedDna {
            bits: new_bits,
            len: self.len,
        }
    }

    /// The raw packed `u64` value.
    ///
    /// The `u64` alone does not carry length; the caller is responsible for
    /// tracking length separately (as `UmiDeduper` does in its config).
    /// Exposed as `pub(crate)` so that `UmiDeduper` can access the encoding.
    #[inline]
    pub(crate) fn bits(self) -> u64 {
        self.bits
    }
}

// ─── MatchResult ──────────────────────────────────────────────────────────────

/// The result of looking up a sequence window in a guide index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchResult {
    /// Exact match; the payload is the guide's 0-based library index.
    Exact(usize),
    /// Unique one-mismatch match; the payload is the guide's 0-based library index.
    OneMismatch(usize),
    /// Two or more guides match at one mismatch — the read is ambiguous and
    /// must be silently skipped (project invariant §5).
    Ambiguous,
    /// No match at any mismatch level.
    NoMatch,
}

// ─── PackedGuideIndex ─────────────────────────────────────────────────────────

/// Guide index for libraries where every guide is ≤ 32 nt.
///
/// At construction the caller chooses between exact-only (`mismatch = 0`) and
/// one-mismatch (`mismatch = 1`) modes.  In one-mismatch mode the index also
/// stores all 3 × L substitution variants; if a variant maps to two or more
/// distinct guides it is marked [`MatchResult::Ambiguous`].
///
/// Exact matches always take precedence over one-mismatch matches, as required
/// by the specification.
pub struct PackedGuideIndex {
    /// Primary map: exact `PackedDna` → guide index.
    exact: HashMap<PackedDna, usize>,
    /// Mismatch map: variant `PackedDna` → result (Exact hit wins; stored
    /// separately so that exact map lookup always wins during query).
    mismatch: HashMap<PackedDna, MatchResult>,
    /// Whether one-mismatch mode is enabled.
    one_mismatch: bool,
}

impl PackedGuideIndex {
    /// Construct a packed guide index.
    ///
    /// `guides` is a slice of encoded guide sequences (from the library TSV,
    /// already validated to be ≤ 32 nt).
    /// `one_mismatch` enables one-mismatch candidate generation.
    ///
    /// Returns `Err(usize)` containing the duplicate guide index if a duplicate
    /// sequence is detected.
    pub fn build(guides: &[PackedDna], one_mismatch: bool) -> Result<Self, usize> {
        let mut exact: HashMap<PackedDna, usize> = HashMap::with_capacity(guides.len());
        for (idx, &key) in guides.iter().enumerate() {
            if exact.insert(key, idx).is_some() {
                return Err(idx);
            }
        }

        let mut mismatch: HashMap<PackedDna, MatchResult> = HashMap::new();
        if one_mismatch {
            // Build substitution variants for every guide.
            for (idx, &key) in guides.iter().enumerate() {
                let l = key.len();
                for pos in 0..l {
                    let orig_code = (key.bits() >> (2 * pos)) & 0b11;
                    // Try all 3 alternative bases at this position.
                    for alt_code in 0u64..4 {
                        if alt_code == orig_code {
                            continue;
                        }
                        let shift = 2 * pos;
                        let variant_bits = (key.bits() & !(0b11u64 << shift)) | (alt_code << shift);
                        let variant = PackedDna {
                            bits: variant_bits,
                            len: key.len,
                        };
                        // Skip if this variant is also an exact guide — exact
                        // lookup wins unconditionally, so we don't need (or want)
                        // a mismatch entry for it.
                        if exact.contains_key(&variant) {
                            continue;
                        }
                        let entry = mismatch
                            .entry(variant)
                            .or_insert(MatchResult::OneMismatch(idx));
                        // If already occupied by a different guide → ambiguous.
                        if let MatchResult::OneMismatch(prev_idx) = *entry {
                            if prev_idx != idx {
                                *entry = MatchResult::Ambiguous;
                            }
                        }
                    }
                }
            }
        }

        Ok(PackedGuideIndex {
            exact,
            mismatch,
            one_mismatch,
        })
    }

    /// Look up a query sequence.
    ///
    /// 1. Exact match → `MatchResult::Exact(idx)` (always takes precedence).
    /// 2. If one-mismatch mode is enabled: unique one-mismatch → `OneMismatch(idx)`;
    ///    ambiguous → `Ambiguous`.
    /// 3. Otherwise → `NoMatch`.
    pub fn lookup(&self, query: &PackedDna) -> MatchResult {
        if let Some(&idx) = self.exact.get(query) {
            return MatchResult::Exact(idx);
        }
        if self.one_mismatch {
            if let Some(&result) = self.mismatch.get(query) {
                return result;
            }
        }
        MatchResult::NoMatch
    }

    /// The number of guides in the index.
    pub fn len(&self) -> usize {
        self.exact.len()
    }

    /// Returns `true` when the index contains no guides.
    pub fn is_empty(&self) -> bool {
        self.exact.is_empty()
    }
}

// ─── ByteGuideIndex ───────────────────────────────────────────────────────────

/// Guide index fallback for libraries that contain any guide longer than 32 nt.
///
/// Uses `Vec<u8>` (uppercase-normalised bytes) as the key.  Exact-match and
/// one-mismatch semantics are identical to [`PackedGuideIndex`].
pub struct ByteGuideIndex {
    /// Primary map: exact byte sequence → guide index.
    exact: HashMap<Vec<u8>, usize>,
    /// Mismatch map (populated only in one-mismatch mode).
    mismatch: HashMap<Vec<u8>, MatchResult>,
    /// Whether one-mismatch mode is enabled.
    one_mismatch: bool,
}

impl ByteGuideIndex {
    /// Construct a byte guide index.
    ///
    /// `guides` is a slice of raw guide byte sequences (from the library TSV,
    /// already uppercase-normalised and validated).
    /// `one_mismatch` enables one-mismatch candidate generation.
    ///
    /// Returns `Err(usize)` containing the duplicate guide index if a duplicate
    /// sequence is detected.
    pub fn build(guides: &[Vec<u8>], one_mismatch: bool) -> Result<Self, usize> {
        let mut exact: HashMap<Vec<u8>, usize> = HashMap::with_capacity(guides.len());
        for (idx, key) in guides.iter().enumerate() {
            // Use borrowed lookup to avoid cloning during the duplicate check.
            if exact.contains_key(key.as_slice()) {
                return Err(idx);
            }
            exact.insert(key.clone(), idx);
        }

        let mut mismatch: HashMap<Vec<u8>, MatchResult> = HashMap::new();
        if one_mismatch {
            for (idx, key) in guides.iter().enumerate() {
                let l = key.len();
                for pos in 0..l {
                    let orig = key[pos];
                    for &alt in b"ACGT" {
                        if alt == orig {
                            continue;
                        }
                        // Build variant — one allocation per candidate, avoided
                        // during lookup via borrowed HashMap access.
                        let mut variant = key.clone();
                        variant[pos] = alt;
                        // Skip if the variant is an exact guide key.
                        if exact.contains_key(variant.as_slice()) {
                            continue;
                        }
                        let entry = mismatch
                            .entry(variant)
                            .or_insert(MatchResult::OneMismatch(idx));
                        if let MatchResult::OneMismatch(prev_idx) = *entry {
                            if prev_idx != idx {
                                *entry = MatchResult::Ambiguous;
                            }
                        }
                    }
                }
            }
        }

        Ok(ByteGuideIndex {
            exact,
            mismatch,
            one_mismatch,
        })
    }

    /// Look up a query byte sequence (must be already uppercase-normalised).
    ///
    /// Uses borrowed lookups to avoid heap allocation per query.
    pub fn lookup(&self, query: &[u8]) -> MatchResult {
        if let Some(&idx) = self.exact.get(query) {
            return MatchResult::Exact(idx);
        }
        if self.one_mismatch {
            if let Some(&result) = self.mismatch.get(query) {
                return result;
            }
        }
        MatchResult::NoMatch
    }

    /// The number of guides in the index.
    pub fn len(&self) -> usize {
        self.exact.len()
    }

    /// Returns `true` when the index contains no guides.
    pub fn is_empty(&self) -> bool {
        self.exact.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enc(s: &str) -> PackedDna {
        PackedDna::encode(s.as_bytes()).expect("valid sequence")
    }

    #[test]
    fn single_base_a() {
        let p = enc("A");
        assert_eq!(p.len(), 1);
        assert_eq!(p.decode(), b"A");
        assert_eq!(p.bits() & 0b11, 0b00);
    }

    #[test]
    fn single_base_c() {
        let p = enc("C");
        assert_eq!(p.len(), 1);
        assert_eq!(p.decode(), b"C");
        assert_eq!(p.bits() & 0b11, 0b01);
    }

    #[test]
    fn single_base_g() {
        let p = enc("G");
        assert_eq!(p.len(), 1);
        assert_eq!(p.decode(), b"G");
        assert_eq!(p.bits() & 0b11, 0b10);
    }

    #[test]
    fn single_base_t() {
        let p = enc("T");
        assert_eq!(p.len(), 1);
        assert_eq!(p.decode(), b"T");
        assert_eq!(p.bits() & 0b11, 0b11);
    }

    #[test]
    fn binary_mapping_a_eq_00() {
        assert_eq!(enc("A").bits() & 0b11, 0b00);
    }

    #[test]
    fn binary_mapping_c_eq_01() {
        assert_eq!(enc("C").bits() & 0b11, 0b01);
    }

    #[test]
    fn binary_mapping_g_eq_10() {
        assert_eq!(enc("G").bits() & 0b11, 0b10);
    }

    #[test]
    fn binary_mapping_t_eq_11() {
        assert_eq!(enc("T").bits() & 0b11, 0b11);
    }

    #[test]
    fn first_base_in_lsb_pair() {
        let p = enc("CA");
        assert_eq!(p.bits() & 0b11, 0b01, "first base C must be in bits[1:0]");
        assert_eq!(
            (p.bits() >> 2) & 0b11,
            0b00,
            "second base A must be in bits[3:2]"
        );
    }

    #[test]
    fn base_ordering_acgt() {
        let p = enc("ACGT");
        assert_eq!(p.bits() & 0b11, 0b00, "A at pos 0");
        assert_eq!((p.bits() >> 2) & 0b11, 0b01, "C at pos 1");
        assert_eq!((p.bits() >> 4) & 0b11, 0b10, "G at pos 2");
        assert_eq!((p.bits() >> 6) & 0b11, 0b11, "T at pos 3");
    }

    #[test]
    fn bits_beyond_length_are_zero() {
        let p = enc("A");
        assert_eq!(p.bits() >> 2, 0);
    }

    #[test]
    fn length_distinguishes_otherwise_identical_bits() {
        let a = enc("A");
        let aa = enc("AA");
        assert_eq!(a.bits(), 0, "A should have bits=0");
        assert_eq!(aa.bits(), 0, "AA should have bits=0");
        assert_ne!(a, aa);
    }
}
