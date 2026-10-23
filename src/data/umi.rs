//! UMI (Unique Molecular Identifier) deduplication.
//!
//! [`UmiDeduper`] maintains per-group sets of already-seen UMIs and answers
//! whether a newly observed UMI has been seen before in the same group.
//!
//! The storage representation is chosen once at construction based on the
//! configured UMI length:
//!
//! - **Packed mode** (`umi_len ≤ 32`): each group's set is a
//!   `HashSet<u64>` holding the packed two-bit encoding produced by
//!   [`PackedDna`].
//! - **Byte mode** (`umi_len > 32`): each group's set is a
//!   `HashSet<Vec<u8>>` holding the uppercase-normalised raw bytes.
//!
//! The representation cannot change after construction.

use std::collections::HashSet;

use crate::data::packed_dna::{EncodingError, PackedDna};

// ─── Configuration error ──────────────────────────────────────────────────────

/// Errors returned by [`UmiDeduper`] operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UmiError {
    /// A UMI byte is not A, C, G, or T (after case normalisation).
    InvalidBase(u8),
    /// Inserting a new key would exceed the active scope cap.
    KeyCapExceeded { limit: usize },
    /// The supplied UMI length does not match the configured fixed length.
    LengthMismatch { expected: usize, got: usize },
    /// The group index is out of range.
    GroupOutOfRange { group: usize, n_groups: usize },
    /// The UMI is empty but a non-zero length was configured.
    Empty,
}

impl std::fmt::Display for UmiError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            UmiError::KeyCapExceeded { limit } => write!(f, "UMI unique-key cap {limit} exceeded"),
            UmiError::InvalidBase(b) => write!(f, "invalid UMI base byte 0x{b:02x}"),
            UmiError::LengthMismatch { expected, got } => {
                write!(f, "UMI length mismatch: expected {expected}, got {got}")
            }
            UmiError::GroupOutOfRange { group, n_groups } => {
                write!(
                    f,
                    "group index {group} out of range (n_groups = {n_groups})"
                )
            }
            UmiError::Empty => write!(f, "UMI is empty but a non-zero length is configured"),
        }
    }
}

impl std::error::Error for UmiError {}

impl From<EncodingError> for UmiError {
    fn from(e: EncodingError) -> Self {
        match e {
            EncodingError::InvalidBase(b) => UmiError::InvalidBase(b),
            // TooLong should never reach here — packed mode is only selected
            // for umi_len ≤ 32, so encode() cannot return TooLong.
            EncodingError::TooLong(_) => {
                unreachable!("PackedDna::encode called with len > 32 in packed mode")
            }
        }
    }
}

// ─── Internal storage ─────────────────────────────────────────────────────────

/// Type-erased per-group storage.
enum Storage {
    /// Packed mode: one `HashSet<u64>` per group.
    Packed(Vec<HashSet<u64>>),
    /// Byte mode: one `HashSet<Vec<u8>>` per group.
    Byte(Vec<HashSet<Vec<u8>>>),
}

// ─── UmiDeduper ───────────────────────────────────────────────────────────────

/// Deduplicates UMI observations within independent groups (e.g. per sgRNA).
///
/// # Construction
///
/// ```
/// use mageckrs::data::umi::UmiDeduper;
///
/// // 5 groups (e.g. 5 guides), 12-nt UMIs (packed mode)
/// let deduper = UmiDeduper::new(5, 12);
/// ```
///
/// # Group semantics
///
/// Each group is independent: the same UMI sequence observed in different
/// groups does not count as a duplicate.
///
/// # Fixed-length enforcement
///
/// When `umi_len > 0` every submitted UMI must have exactly `umi_len` bytes
/// (after case normalisation).  A length mismatch returns
/// [`UmiError::LengthMismatch`] without inserting.
///
/// When `umi_len == 0` empty UMIs are accepted (all-zero bits packed, or
/// empty `Vec<u8>`).
pub struct UmiDeduper {
    /// Configured UMI length (0 = unset / accept-empty).
    umi_len: usize,
    unique_keys: usize,
    /// Type-dispatched storage.
    storage: Storage,
}

impl UmiDeduper {
    /// Construct a new `UmiDeduper` with `n_groups` independent deduplication
    /// sets, each expecting UMIs of exactly `umi_len` nucleotides.
    ///
    /// The storage mode is selected at this point and does not change.
    pub fn new(n_groups: usize, umi_len: usize) -> Self {
        let storage = if umi_len <= PackedDna::MAX_LEN {
            Storage::Packed((0..n_groups).map(|_| HashSet::new()).collect())
        } else {
            Storage::Byte((0..n_groups).map(|_| HashSet::new()).collect())
        };
        UmiDeduper {
            umi_len,
            storage,
            unique_keys: 0,
        }
    }

    /// The number of independent deduplication groups.
    pub fn n_groups(&self) -> usize {
        match &self.storage {
            Storage::Packed(v) => v.len(),
            Storage::Byte(v) => v.len(),
        }
    }

    /// The configured UMI length (0 means unset).
    pub fn umi_len(&self) -> usize {
        self.umi_len
    }

    /// Returns `true` when UMIs are stored in packed (u64) mode.
    pub fn is_packed(&self) -> bool {
        matches!(&self.storage, Storage::Packed(_))
    }

    /// Observe `umi` in `group`.
    ///
    /// Returns `Ok(true)` if this is the first observation of `umi` in
    /// `group` (i.e. it was **newly inserted**), or `Ok(false)` if the UMI
    /// was already present (i.e. it is a **duplicate**).
    ///
    /// # Errors
    ///
    /// - [`UmiError::GroupOutOfRange`] — `group` ≥ `n_groups`.
    /// - [`UmiError::LengthMismatch`] — the UMI length does not match the
    ///   configured `umi_len` (only checked when `umi_len > 0`).
    /// - [`UmiError::Empty`] — the UMI is empty and `umi_len > 0`.
    /// - [`UmiError::InvalidBase`] — a byte is not ACGT.
    pub fn observe(&mut self, group: usize, umi: &[u8]) -> Result<bool, UmiError> {
        self.observe_capped(group, umi, usize::MAX)
    }

    /// Observe exactly; validate and check duplicates before enforcing the total key cap.
    /// A rejected key is never inserted. The cap does not bound process memory.
    pub fn observe_capped(
        &mut self,
        group: usize,
        umi: &[u8],
        limit: usize,
    ) -> Result<bool, UmiError> {
        let n_groups = self.n_groups();
        if group >= n_groups {
            return Err(UmiError::GroupOutOfRange { group, n_groups });
        }

        // Length validation.
        if self.umi_len == 0 {
            if !umi.is_empty() {
                return Err(UmiError::LengthMismatch {
                    expected: 0,
                    got: umi.len(),
                });
            }
        } else {
            if umi.is_empty() {
                return Err(UmiError::Empty);
            }
            if umi.len() != self.umi_len {
                return Err(UmiError::LengthMismatch {
                    expected: self.umi_len,
                    got: umi.len(),
                });
            }
        }

        match &mut self.storage {
            Storage::Packed(sets) => {
                // Validate and encode.
                let packed = PackedDna::encode(
                    // Uppercase-normalise before encoding.
                    &umi.iter()
                        .map(|b| b.to_ascii_uppercase())
                        .collect::<Vec<_>>(),
                )?;
                let key = packed.bits();
                if sets[group].contains(&key) {
                    return Ok(false);
                }
                if self.unique_keys >= limit {
                    return Err(UmiError::KeyCapExceeded { limit });
                }
                sets[group].insert(key);
                self.unique_keys += 1;
                Ok(true)
            }
            Storage::Byte(sets) => {
                // Validate each byte and build an uppercase copy.
                let mut normalised = Vec::with_capacity(umi.len());
                for &raw in umi {
                    let b = raw.to_ascii_uppercase();
                    match b {
                        b'A' | b'C' | b'G' | b'T' => normalised.push(b),
                        _ => return Err(UmiError::InvalidBase(raw)),
                    }
                }
                if sets[group].contains(&normalised) {
                    return Ok(false);
                }
                if self.unique_keys >= limit {
                    return Err(UmiError::KeyCapExceeded { limit });
                }
                sets[group].insert(normalised);
                self.unique_keys += 1;
                Ok(true)
            }
        }
    }
}
