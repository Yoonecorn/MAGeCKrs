//! Control-sgRNA list parsing and deterministic library resolution.
//!
//! Implements Phase 3H of the count pipeline:
//! - UTF-8 text file reading and line-by-line parsing.
//! - Blank lines and full-line comments (`#`) skipped.
//! - Rejection of internal whitespace (spaces, tabs, carriage returns, newlines).
//! - Exact, case-sensitive independent resolution against sgRNA IDs and sequences.
//! - Ambiguity detection (token matching sequence and different ID, or multiple records).
//! - Deduplicated resolved guide indexes returned in ascending library-row order.
//! - Unmatched entries preserved in original control-file order for warning diagnostics.
//! - Typed errors for empty resolution, malformed entries, ambiguity, or file read/encoding failures.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use crate::count::matcher::SgRnaRecord;
use crate::error::MageckError;

/// Resolved control sgRNA guide indexes and unmatched entries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedControlGuides {
    /// 0-based library row indices of matched control sgRNAs, deduplicated and in ascending library-row order.
    pub guide_indexes: Vec<usize>,
    /// Unmatched prospective control entries in original control-file order, preserving duplicates.
    pub unknown_entries: Vec<String>,
}

/// Parses control sgRNA identifiers from a reader implementing [`BufRead`] and resolves them against library records.
///
/// Leading and trailing whitespace is trimmed. Blank lines and full-line comments (`#`) are skipped.
/// Non-comment entries containing internal whitespace are rejected with [`MageckError::Config`].
/// Each entry is resolved independently against exact sgRNA ID and exact sequence.
/// Ambiguous entries matching multiple records return [`MageckError::Config`].
/// Resolved indexes are deduplicated and returned in ascending library-row order.
/// If no guide indices resolve, returns [`MageckError::Config`].
pub fn parse_control_sgrna<R: BufRead>(
    reader: R,
    records: &[SgRnaRecord],
) -> Result<ResolvedControlGuides, MageckError> {
    parse_control_sgrna_inner(reader, records, None)
}

/// Parses a control sgRNA file from disk and resolves entries against library records.
///
/// Returns [`MageckError::Io`] on file missing, file access errors, or invalid UTF-8 encoding.
/// Rejects non-comment entries containing internal whitespace with [`MageckError::Config`].
/// Ambiguous entries and empty resolutions return [`MageckError::Config`].
pub fn parse_control_sgrna_file<P: AsRef<Path>>(
    path: P,
    records: &[SgRnaRecord],
) -> Result<ResolvedControlGuides, MageckError> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref).map_err(|err| {
        MageckError::Io(io::Error::new(
            err.kind(),
            format!(
                "Failed to open control sgRNA file '{}': {err}",
                path_ref.display()
            ),
        ))
    })?;
    let reader = BufReader::new(file);
    parse_control_sgrna_inner(reader, records, Some(path_ref))
}

fn parse_control_sgrna_inner<R: BufRead>(
    reader: R,
    records: &[SgRnaRecord],
    file_path: Option<&Path>,
) -> Result<ResolvedControlGuides, MageckError> {
    let mut line_num = 0usize;
    let mut prospective_tokens = Vec::new();

    for line_res in reader.lines() {
        line_num += 1;
        let line = line_res.map_err(|err| {
            MageckError::Io(io::Error::new(
                err.kind(),
                match file_path {
                    Some(path) => format!(
                        "Failed to read line {line_num} from control sgRNA file '{}': {err}",
                        path.display()
                    ),
                    None => {
                        format!("Failed to read line {line_num} from control sgRNA input: {err}")
                    }
                },
            ))
        })?;

        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed.chars().any(char::is_whitespace) {
            return Err(MageckError::Config(format!(
                "Line {line_num}: malformed control sgRNA entry '{trimmed}' contains internal whitespace"
            )));
        }

        prospective_tokens.push(trimmed.to_string());
    }

    if records.is_empty() {
        return Err(MageckError::Config(match file_path {
            Some(path) => format!(
                "No control sgRNAs in '{}' matched any library record (library contains no records)",
                path.display()
            ),
            None => "No control sgRNAs matched any library record (library contains no records)".to_string(),
        }));
    }

    // Precompute exact lookup maps for ID and sequence
    let mut id_map: HashMap<&str, Vec<usize>> = HashMap::with_capacity(records.len());
    let mut seq_map: HashMap<&str, Vec<usize>> = HashMap::with_capacity(records.len());

    for (idx, record) in records.iter().enumerate() {
        id_map.entry(record.id.as_str()).or_default().push(idx);
        seq_map
            .entry(record.sequence.as_str())
            .or_default()
            .push(idx);
    }

    let mut resolved_indices = Vec::new();
    let mut unknown_entries = Vec::new();

    for token in prospective_tokens {
        let id_matches = id_map.get(token.as_str());
        let seq_matches = seq_map.get(token.as_str());

        let mut matched_indices = Vec::new();
        if let Some(matches) = id_matches {
            matched_indices.extend_from_slice(matches);
        }
        if let Some(matches) = seq_matches {
            matched_indices.extend_from_slice(matches);
        }

        matched_indices.sort_unstable();
        matched_indices.dedup();

        match matched_indices.len() {
            0 => {
                unknown_entries.push(token);
            }
            1 => {
                resolved_indices.push(matched_indices[0]);
            }
            _ => {
                let matched_descs: Vec<String> = matched_indices
                    .iter()
                    .map(|&idx| format!("row {idx} ('{}')", records[idx].id))
                    .collect();
                return Err(MageckError::Config(format!(
                    "Control sgRNA entry '{token}' is ambiguous: matches multiple library records ({})",
                    matched_descs.join(", ")
                )));
            }
        }
    }

    resolved_indices.sort_unstable();
    resolved_indices.dedup();

    if resolved_indices.is_empty() {
        return Err(MageckError::Config(match file_path {
            Some(path) => format!(
                "No control sgRNAs from control file '{}' matched the library",
                path.display()
            ),
            None => "No control sgRNAs matched the library".to_string(),
        }));
    }

    Ok(ResolvedControlGuides {
        guide_indexes: resolved_indices,
        unknown_entries,
    })
}
