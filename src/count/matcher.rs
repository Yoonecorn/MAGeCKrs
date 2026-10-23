use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::data::{ByteGuideIndex, MatchResult, PackedDna, PackedGuideIndex};
use crate::error::MageckError;

/// A single sgRNA record from a library definition file.
///
/// Corresponds to specification §8.2: retains the original identifier,
/// target gene symbol/name, and nucleotide sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SgRnaRecord {
    /// Unique sgRNA identifier (e.g., "s_1001", "A1CF_guide1").
    pub id: String,
    /// Target gene symbol or name (e.g., "A1CF", "Control").
    pub gene: String,
    /// Uppercase nucleotide sequence (e.g., "ACGTACGTACGTACGTACGT").
    pub sequence: String,
}

/// In-memory representation of a parsed sgRNA reference library.
///
/// Preserves the original file order of the guides and guarantees
/// uniqueness of guide identifiers and sequences.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Library {
    records: Vec<SgRnaRecord>,
}

impl Library {
    /// Creates a new `Library` from a list of records.
    ///
    /// Validates uniqueness of guide IDs and sequences, validates ACGT bases,
    /// and ensures non-empty records.
    pub fn new(records: Vec<SgRnaRecord>) -> Result<Self, MageckError> {
        if records.is_empty() {
            return Err(MageckError::Library(
                "Library must contain at least one sgRNA record".to_string(),
            ));
        }

        let mut seen_ids: HashMap<String, usize> = HashMap::with_capacity(records.len());
        let mut seen_seqs: HashMap<String, (String, usize)> = HashMap::with_capacity(records.len());

        for (idx, record) in records.iter().enumerate() {
            let row_num = idx + 1;
            if record.id.trim().is_empty() {
                return Err(MageckError::Library(format!(
                    "Record {row_num}: guide identifier is empty"
                )));
            }
            if record.gene.trim().is_empty() {
                return Err(MageckError::Library(format!(
                    "Record {row_num}: gene identifier is empty for guide '{}'",
                    record.id
                )));
            }
            if record.sequence.trim().is_empty() {
                return Err(MageckError::Library(format!(
                    "Record {row_num}: guide sequence is empty for guide '{}'",
                    record.id
                )));
            }

            for b in record.sequence.bytes() {
                match b {
                    b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' => {}
                    _ => {
                        return Err(MageckError::Library(format!(
                            "Record {row_num}: invalid nucleotide '{}' in guide sequence for guide '{}'",
                            b as char, record.id
                        )));
                    }
                }
            }

            let upper_seq = record.sequence.to_ascii_uppercase();

            if let Some(&first_row) = seen_ids.get(&record.id) {
                return Err(MageckError::Library(format!(
                    "Record {row_num}: duplicate guide identifier '{}' (first seen at record {first_row})",
                    record.id
                )));
            }
            seen_ids.insert(record.id.clone(), row_num);

            if let Some((first_id, first_row)) = seen_seqs.get(&upper_seq) {
                return Err(MageckError::Library(format!(
                    "Record {row_num}: duplicate guide sequence '{upper_seq}' for guide '{}' (first seen at record {first_row} for guide '{first_id}')",
                    record.id
                )));
            }
            seen_seqs.insert(upper_seq, (record.id.clone(), row_num));
        }

        Ok(Self { records })
    }

    /// Read-only slice of all parsed sgRNA records in their original file order.
    #[inline]
    pub fn records(&self) -> &[SgRnaRecord] {
        &self.records
    }

    /// Number of sgRNAs in the library.
    #[inline]
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Returns `true` if the library contains zero sgRNAs.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Get a reference to the record at 0-based index `idx`.
    #[inline]
    pub fn get(&self, idx: usize) -> Option<&SgRnaRecord> {
        self.records.get(idx)
    }

    /// Build a [`GuideIndex`] from this library.
    ///
    /// Delegates to [`build_index`].
    pub fn build_index(&self, one_mismatch: bool) -> Result<GuideIndex, MageckError> {
        build_index(self, one_mismatch)
    }
}

/// Unified guide index enum wrapping either a [`PackedGuideIndex`] (for guides <= 32 nt)
/// or a [`ByteGuideIndex`] fallback (when any guide exceeds 32 nt).
///
/// Implements spec §18 Phase 3:
/// `src/count/matcher.rs: GuideIndex enum wrapping PackedGuideIndex or ByteGuideIndex; build_index()`
pub enum GuideIndex {
    /// All guides in the library are <= 32 nt; bit-packed representation with no per-query heap allocations.
    Packed(PackedGuideIndex),
    /// At least one guide in the library exceeds 32 nt; byte slice fallback representation.
    Byte(ByteGuideIndex),
}

impl std::fmt::Debug for GuideIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GuideIndex::Packed(idx) => f
                .debug_struct("GuideIndex::Packed")
                .field("len", &idx.len())
                .finish(),
            GuideIndex::Byte(idx) => f
                .debug_struct("GuideIndex::Byte")
                .field("len", &idx.len())
                .finish(),
        }
    }
}

impl GuideIndex {
    /// Look up a query nucleotide sequence slice.
    ///
    /// Returns:
    /// - [`MatchResult::Exact`] if the query is an exact match for guide `idx`.
    /// - [`MatchResult::OneMismatch`] if the query uniquely matches guide `idx` with 1 mismatch.
    /// - [`MatchResult::Ambiguous`] if the query matches multiple guides with 1 mismatch.
    /// - [`MatchResult::NoMatch`] if the query has no match (or invalid characters).
    pub fn lookup(&self, query: &[u8]) -> MatchResult {
        match self {
            GuideIndex::Packed(idx) => {
                if query.is_empty() || query.len() > 32 {
                    MatchResult::NoMatch
                } else {
                    match PackedDna::encode(query) {
                        Ok(packed) => idx.lookup(&packed),
                        Err(_) => MatchResult::NoMatch,
                    }
                }
            }
            GuideIndex::Byte(idx) => {
                let mut needs_upper = false;
                for &b in query {
                    if b.is_ascii_lowercase() {
                        needs_upper = true;
                        break;
                    }
                }
                if needs_upper {
                    let upper = query.to_ascii_uppercase();
                    idx.lookup(&upper)
                } else {
                    idx.lookup(query)
                }
            }
        }
    }

    /// Look up a pre-encoded [`PackedDna`] query.
    pub fn lookup_packed(&self, query: &PackedDna) -> MatchResult {
        match self {
            GuideIndex::Packed(idx) => idx.lookup(query),
            GuideIndex::Byte(idx) => {
                let decoded = query.decode();
                idx.lookup(&decoded)
            }
        }
    }

    /// Number of guides indexed.
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            GuideIndex::Packed(idx) => idx.len(),
            GuideIndex::Byte(idx) => idx.len(),
        }
    }

    /// Returns `true` if the index contains zero guides.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns `true` if this index is using the 64-bit packed representation.
    #[inline]
    pub fn is_packed(&self) -> bool {
        matches!(self, GuideIndex::Packed(_))
    }

    /// Returns `true` if this index is using the byte fallback representation.
    #[inline]
    pub fn is_byte(&self) -> bool {
        matches!(self, GuideIndex::Byte(_))
    }
}

/// Internal helper storing 0-based column indices for required fields.
#[derive(Debug, Clone, Copy)]
struct ColumnMapping {
    id_col: usize,
    seq_col: usize,
    gene_col: usize,
    min_cols: usize,
}

/// Detect delimiter: prefers tab `\t`, falls back to comma `,` if tab is absent.
fn detect_delimiter(first_line: &str) -> char {
    if first_line.contains('\t') {
        '\t'
    } else if first_line.contains(',') {
        ','
    } else {
        '\t'
    }
}

/// Parse and validate the library file header row.
fn parse_header(
    header_line: &str,
    line_num: usize,
    delimiter: char,
) -> Result<ColumnMapping, MageckError> {
    let cols: Vec<&str> = header_line.split(delimiter).map(|s| s.trim()).collect();
    let mut id_col = None;
    let mut seq_col = None;
    let mut gene_col = None;

    for (idx, &col) in cols.iter().enumerate() {
        let lower = col.to_ascii_lowercase();
        if id_col.is_none() && lower == "id" {
            id_col = Some(idx);
        } else if seq_col.is_none() && lower == "sequence" {
            seq_col = Some(idx);
        } else if gene_col.is_none() && lower == "gene" {
            gene_col = Some(idx);
        }
    }

    let mut missing = Vec::new();
    if id_col.is_none() {
        missing.push("id");
    }
    if seq_col.is_none() {
        missing.push("sequence");
    }
    if gene_col.is_none() {
        missing.push("gene");
    }

    if !missing.is_empty() {
        return Err(MageckError::Library(format!(
            "Line {line_num}: library header is missing required column(s): {}",
            missing.join(", ")
        )));
    }

    let id = id_col.unwrap();
    let seq = seq_col.unwrap();
    let gene = gene_col.unwrap();
    let min_cols = id.max(seq).max(gene) + 1;

    Ok(ColumnMapping {
        id_col: id,
        seq_col: seq,
        gene_col: gene,
        min_cols,
    })
}

/// Parse an sgRNA reference library from a reader implementing [`BufRead`].
///
/// Accepts tab-separated or comma-separated formats, detects headers case-insensitively,
/// handles Unix (`\n`) and CRLF (`\r\n`) line endings, converts sequence case to uppercase,
/// validates nucleotide characters (A, C, G, T), enforces uniqueness of guide IDs and sequences,
/// and reports 1-based line numbers in error messages.
pub fn parse_library<R: BufRead>(reader: R) -> Result<Library, MageckError> {
    let mut lines = reader.lines();
    let mut line_num = 0usize;

    // 1. Locate the first non-blank line (must be the header).
    let mut header_mapping = None;
    let mut delimiter = '\t';

    for line_res in lines.by_ref() {
        line_num += 1;
        let line = line_res.map_err(MageckError::Io)?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        delimiter = detect_delimiter(trimmed);
        let mapping = parse_header(trimmed, line_num, delimiter)?;
        header_mapping = Some(mapping);
        break;
    }

    let mapping = match header_mapping {
        Some(m) => m,
        None => {
            return Err(MageckError::Library("Library file is empty".to_string()));
        }
    };

    // 2. Parse data rows.
    let mut records = Vec::new();
    let mut seen_ids: HashMap<String, usize> = HashMap::new();
    let mut seen_seqs: HashMap<String, (String, usize)> = HashMap::new();

    for line_res in lines {
        line_num += 1;
        let line = line_res.map_err(MageckError::Io)?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        let cols: Vec<&str> = line.split(delimiter).map(|s| s.trim()).collect();
        if cols.len() < mapping.min_cols {
            return Err(MageckError::Library(format!(
                "Line {line_num}: insufficient columns (expected at least {}, got {})",
                mapping.min_cols,
                cols.len()
            )));
        }

        let id = cols[mapping.id_col];
        if id.is_empty() {
            return Err(MageckError::Library(format!(
                "Line {line_num}: guide identifier is empty"
            )));
        }

        let gene = cols[mapping.gene_col];
        if gene.is_empty() {
            return Err(MageckError::Library(format!(
                "Line {line_num}: gene identifier is empty for guide '{id}'"
            )));
        }

        let raw_seq = cols[mapping.seq_col];
        if raw_seq.is_empty() {
            return Err(MageckError::Library(format!(
                "Line {line_num}: guide sequence is empty for guide '{id}'"
            )));
        }

        // Validate nucleotide symbols.
        for b in raw_seq.bytes() {
            match b {
                b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' => {}
                _ => {
                    return Err(MageckError::Library(format!(
                        "Line {line_num}: invalid nucleotide '{}' in guide sequence for guide '{id}'",
                        b as char
                    )));
                }
            }
        }

        let seq = raw_seq.to_ascii_uppercase();

        // Check duplicate identifier.
        if let Some(&first_line) = seen_ids.get(id) {
            return Err(MageckError::Library(format!(
                "Line {line_num}: duplicate guide identifier '{id}' (first seen on line {first_line})"
            )));
        }
        seen_ids.insert(id.to_string(), line_num);

        // Check duplicate sequence.
        if let Some((first_id, first_line)) = seen_seqs.get(&seq) {
            return Err(MageckError::Library(format!(
                "Line {line_num}: duplicate guide sequence '{seq}' for guide '{id}' (first seen on line {first_line} for guide '{first_id}')"
            )));
        }
        seen_seqs.insert(seq.clone(), (id.to_string(), line_num));

        records.push(SgRnaRecord {
            id: id.to_string(),
            gene: gene.to_string(),
            sequence: seq,
        });
    }

    if records.is_empty() {
        return Err(MageckError::Library(
            "Library file contains no guide records (header only)".to_string(),
        ));
    }

    Ok(Library { records })
}

/// Parse an sgRNA reference library directly from a file path.
pub fn parse_library_file<P: AsRef<Path>>(path: P) -> Result<Library, MageckError> {
    let file = File::open(path.as_ref()).map_err(MageckError::Io)?;
    let reader = BufReader::new(file);
    parse_library(reader)
}

/// Build a unified [`GuideIndex`] from a parsed [`Library`].
///
/// If every guide sequence is <= 32 nt, builds a [`PackedGuideIndex`] ([`GuideIndex::Packed`]).
/// If any guide sequence exceeds 32 nt, builds a [`ByteGuideIndex`] fallback ([`GuideIndex::Byte`]).
/// Invariant: Sequences longer than 32 nt are never truncated.
pub fn build_index(library: &Library, one_mismatch: bool) -> Result<GuideIndex, MageckError> {
    if library.is_empty() {
        return Err(MageckError::Library(
            "Cannot build index for empty library".to_string(),
        ));
    }

    let any_long = library.records().iter().any(|r| r.sequence.len() > 32);

    if any_long {
        let byte_guides: Vec<Vec<u8>> = library
            .records()
            .iter()
            .map(|r| r.sequence.as_bytes().to_vec())
            .collect();
        let idx = ByteGuideIndex::build(&byte_guides, one_mismatch).map_err(|dup_idx| {
            MageckError::Library(format!(
                "Duplicate guide sequence in library at record index {dup_idx} ('{}')",
                library.records()[dup_idx].id
            ))
        })?;
        Ok(GuideIndex::Byte(idx))
    } else {
        let mut packed_guides = Vec::with_capacity(library.len());
        for r in library.records() {
            let packed = PackedDna::encode(r.sequence.as_bytes()).map_err(|e| {
                MageckError::Library(format!(
                    "Failed to encode guide '{}' (sequence '{}'): {e}",
                    r.id, r.sequence
                ))
            })?;
            packed_guides.push(packed);
        }
        let idx = PackedGuideIndex::build(&packed_guides, one_mismatch).map_err(|dup_idx| {
            MageckError::Library(format!(
                "Duplicate guide sequence in library at record index {dup_idx} ('{}')",
                library.records()[dup_idx].id
            ))
        })?;
        Ok(GuideIndex::Packed(idx))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_minimal_valid() {
        let tsv = "id\tsequence\tgene\ns1\tACGT\tGENE1\n";
        let lib = parse_library(Cursor::new(tsv)).expect("parse library");
        assert_eq!(lib.len(), 1);
        assert_eq!(lib.records()[0].id, "s1");
        assert_eq!(lib.records()[0].sequence, "ACGT");
        assert_eq!(lib.records()[0].gene, "GENE1");
    }

    #[test]
    fn test_parse_case_normalization() {
        let tsv = "id\tsequence\tgene\ns1\tacgt\tGeneA\ns2\tTgCa\tGeneB\n";
        let lib = parse_library(Cursor::new(tsv)).expect("parse library");
        assert_eq!(lib.len(), 2);
        assert_eq!(lib.records()[0].sequence, "ACGT");
        assert_eq!(lib.records()[1].sequence, "TGCA");
    }
}
