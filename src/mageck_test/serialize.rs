//! RRA summary formatting and failure-safe publication for one comparison.
#![allow(dead_code)] // Crate-private entry points are consumed by M4 orchestration.

use std::collections::HashSet;
use std::ffi::CString;
use std::fmt::Write as _;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use super::aggregation::AggregatedGeneResult;
use crate::data::ParsedCountTable;
use crate::error::MageckError;

const GENE_HEADER: &str = "id\tnum\tneg|score\tneg|p-value\tneg|fdr\tneg|rank\tneg|lfc\tneg|goodsgrna\tpos|score\tpos|p-value\tpos|fdr\tpos|rank\tpos|lfc\tpos|goodsgrna\n";
const SGRNA_HEADER: &str = "sgRNA\tGene\tcontrol_count\ttreatment_count\tcontrol_mean\ttreat_mean\tLFC\tcontrol_var\tadj_var\tscore\tp.low\tp.high\tp.twosided\tFDR\thigh_in_treatment\n";
const MAX_TEMP_CREATE_ATTEMPTS: usize = 32;
static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

fn finite(value: f64, context: &str) -> Result<f64, MageckError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MageckError::Math(format!(
            "RRA serialization {context}: non-finite value"
        )))
    }
}

fn identifier(value: &str) -> Result<(), MageckError> {
    if value.is_empty() || value.contains(['\t', '\r', '\n']) {
        return Err(MageckError::Config(format!(
            "Invalid TSV identifier: {value:?}"
        )));
    }
    Ok(())
}

/// Retains aggregation row order; each direction has independent 1-based ranks.
pub(crate) fn serialize_gene_summary(
    genes: &[AggregatedGeneResult],
) -> Result<String, MageckError> {
    let mut seen = HashSet::new();
    for gene in genes {
        identifier(&gene.gene_id)?;
        if !seen.insert(&gene.gene_id) {
            return Err(MageckError::Math(format!(
                "Duplicate summary gene {:?}",
                gene.gene_id
            )));
        }
        for value in [
            gene.neg_score,
            gene.neg_pvalue,
            gene.neg_fdr,
            gene.neg_lfc,
            gene.pos_score,
            gene.pos_pvalue,
            gene.pos_fdr,
            gene.pos_lfc,
        ] {
            finite(value, &gene.gene_id)?;
        }
    }
    let ranks = |negative: bool| {
        let mut order: Vec<usize> = (0..genes.len()).collect();
        order.sort_by(|&a, &b| {
            let score = |i: usize| {
                if negative {
                    genes[i].neg_score
                } else {
                    genes[i].pos_score
                }
            };
            // Numeric equality includes signed zero; opaque IDs break score ties.
            score(a)
                .partial_cmp(&score(b))
                .expect("validated finite scores")
                .then_with(|| genes[a].gene_id.cmp(&genes[b].gene_id))
        });
        let mut ranks = vec![0; genes.len()];
        for (rank, index) in order.into_iter().enumerate() {
            ranks[index] = rank + 1;
        }
        ranks
    };
    let neg = ranks(true);
    let pos = ranks(false);
    let mut output = GENE_HEADER.to_owned();
    for (i, g) in genes.iter().enumerate() {
        writeln!(
            output,
            "{}\t{}\t{:.8e}\t{:.8e}\t{:.8e}\t{}\t{:.4}\t{}\t{:.8e}\t{:.8e}\t{:.8e}\t{}\t{:.4}\t{}",
            g.gene_id,
            g.num_sgrnas,
            g.neg_score,
            g.neg_pvalue,
            g.neg_fdr,
            neg[i],
            g.neg_lfc,
            g.neg_goodsgrna,
            g.pos_score,
            g.pos_pvalue,
            g.pos_fdr,
            pos[i],
            g.pos_lfc,
            g.pos_goodsgrna
        )
        .expect("writing to String");
    }
    Ok(output)
}

fn validate_columns(columns: &[usize], samples: usize) -> Result<(), MageckError> {
    let mut seen = HashSet::new();
    if columns.is_empty() || columns.iter().any(|&c| c >= samples || !seen.insert(c)) {
        return Err(MageckError::Math(
            "RRA serialization: empty, duplicate or out-of-range sample selection".into(),
        ));
    }
    Ok(())
}

fn condition(
    table: &ParsedCountTable,
    columns: &[usize],
    guide: usize,
) -> Result<(f64, f64), MageckError> {
    let mut sum = 0.0;
    for &column in columns {
        let count = table
            .counts()
            .get(column, guide)
            .expect("validated dimensions");
        sum = finite(sum + count, "replicate sum")?;
    }
    Ok((sum, finite(sum / columns.len() as f64, "replicate mean")?))
}

/// `scores` are P4F4 masked rank fractions, aligned to physical table rows.
/// Guide empirical p-values and FDR are unavailable in P4F10 and always NA.
pub(crate) fn serialize_sgrna_summary(
    table: &ParsedCountTable,
    control_columns: &[usize],
    treatment_columns: &[usize],
    scores: &[Option<f64>],
) -> Result<String, MageckError> {
    validate_columns(control_columns, table.n_samples())?;
    validate_columns(treatment_columns, table.n_samples())?;
    if control_columns
        .iter()
        .any(|c| treatment_columns.contains(c))
        || scores.len() != table.n_guides()
    {
        return Err(MageckError::Math(
            "RRA serialization: overlapping conditions or mismatched guide scores".into(),
        ));
    }
    let mut output = SGRNA_HEADER.to_owned();
    for (guide, score) in scores.iter().enumerate() {
        let id = &table.sgrna_ids()[guide];
        let gene = &table.genes()[guide];
        identifier(id)?;
        identifier(gene)?;
        let score = match score {
            Some(value) if value.is_finite() && *value > 0.0 && *value <= 1.0 => {
                format!("{value:.8e}")
            }
            Some(_) => return Err(MageckError::Math(format!("Invalid guide rank for {id:?}"))),
            None => "NA".into(),
        };
        let (control_sum, control_mean) = condition(table, control_columns, guide)?;
        let (treat_sum, treat_mean) = condition(table, treatment_columns, guide)?;
        let lfc = finite(
            ((treat_mean + 1.0) / (control_mean + 1.0)).log2(),
            "guide LFC",
        )?;
        let mut squares = 0.0;
        if control_columns.len() > 1 {
            for &column in control_columns {
                let delta = table
                    .counts()
                    .get(column, guide)
                    .expect("validated dimensions")
                    - control_mean;
                squares = finite(squares + delta * delta, "control variance")?;
            }
        }
        let variance = if control_columns.len() > 1 {
            finite(
                squares / (control_columns.len() - 1) as f64,
                "control variance",
            )?
        } else {
            0.0
        };
        writeln!(output, "{id}\t{gene}\t{control_sum:.8e}\t{treat_sum:.8e}\t{control_mean:.8e}\t{treat_mean:.8e}\t{lfc:.4}\t{variance:.8e}\t{:.8e}\t{score}\tNA\tNA\tNA\tNA\t{}",
            0.0, if lfc > 0.0 { "True" } else { "False" }).expect("writing to String");
    }
    Ok(output)
}

trait Filesystem {
    fn create_new(&self, path: &Path) -> io::Result<File>;
    fn write(&self, file: File, contents: &str) -> io::Result<()>;
    fn rename(&self, from: &Path, to: &Path) -> io::Result<()>;
    fn remove(&self, path: &Path) -> io::Result<()>;
}

struct RealFilesystem;
impl Filesystem for RealFilesystem {
    fn create_new(&self, path: &Path) -> io::Result<File> {
        OpenOptions::new().write(true).create_new(true).open(path)
    }
    fn write(&self, file: File, contents: &str) -> io::Result<()> {
        let mut writer = BufWriter::new(file);
        writer.write_all(contents.as_bytes())?;
        writer.flush()
    }
    fn rename(&self, from: &Path, to: &Path) -> io::Result<()> {
        // macOS SDK sys/stdio.h: RENAME_EXCL = 0x00000004.
        // Exclusive rename atomically refuses a destination created concurrently.
        unsafe extern "C" {
            fn renamex_np(
                from: *const std::ffi::c_char,
                to: *const std::ffi::c_char,
                flags: std::ffi::c_uint,
            ) -> std::ffi::c_int;
        }
        let path = |p: &Path| {
            CString::new(p.as_os_str().as_bytes())
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "NUL in output path"))
        };
        let from = path(from)?;
        let to = path(to)?;
        // SAFETY: both C strings live through the call; pointers are NUL-terminated.
        if unsafe { renamex_np(from.as_ptr(), to.as_ptr(), 0x00000004) } == 0 {
            Ok(())
        } else {
            Err(io::Error::last_os_error())
        }
    }
    fn remove(&self, path: &Path) -> io::Result<()> {
        fs::remove_file(path)
    }
}

fn reject_existing(path: &Path) -> io::Result<()> {
    match fs::symlink_metadata(path) {
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!("Output already exists: {}", path.display()),
        )),
        Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(e),
    }
}

fn create_temp<F: Filesystem>(filesystem: &F, target: &Path) -> io::Result<(PathBuf, File)> {
    let name = target
        .file_name()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "Output has no filename"))?;
    for _ in 0..MAX_TEMP_CREATE_ATTEMPTS {
        let mut candidate = name.to_os_string();
        candidate.push(format!(
            ".tmp.{}.{}",
            std::process::id(),
            TEMP_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        let candidate = target.with_file_name(candidate);
        match filesystem.create_new(&candidate) {
            Ok(file) => return Ok((candidate, file)),
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e),
        }
    }
    Err(io::Error::new(
        io::ErrorKind::AlreadyExists,
        "Temporary output creation exhausted 32 attempts",
    ))
}

fn publish_with<F: Filesystem>(
    filesystem: &F,
    prefix: &Path,
    gene: &str,
    sgrna: &str,
) -> Result<(), MageckError> {
    let target = |suffix: &str| {
        let mut path = prefix.as_os_str().to_os_string();
        path.push(suffix);
        PathBuf::from(path)
    };
    let gene_target = target(".gene_summary.txt");
    let sgrna_target = target(".sgrna_summary.txt");
    reject_existing(&gene_target)?;
    reject_existing(&sgrna_target)?;
    let mut owned: Vec<PathBuf> = Vec::new();
    let operation = (|| -> io::Result<()> {
        let (gene_temp, gene_file) = create_temp(filesystem, &gene_target)?;
        owned.push(gene_temp.clone());
        let (sgrna_temp, sgrna_file) = create_temp(filesystem, &sgrna_target)?;
        owned.push(sgrna_temp.clone());
        filesystem.write(gene_file, gene)?;
        filesystem.write(sgrna_file, sgrna)?;
        filesystem.rename(&gene_temp, &gene_target)?;
        owned[0] = gene_target;
        filesystem.rename(&sgrna_temp, &sgrna_target)?;
        Ok(())
    })();
    if let Err(primary) = operation {
        let mut message = primary.to_string();
        for path in owned {
            if let Err(cleanup) = filesystem.remove(&path) {
                message.push_str(&format!(
                    "; cleanup failed for {}: {cleanup}",
                    path.display()
                ));
            }
        }
        return Err(MageckError::Io(io::Error::new(primary.kind(), message)));
    }
    Ok(())
}

/// Formats and validates both payloads before creating any output file.
/// Output parent directory must exist. Publication follows the RUN1 policy.
pub(crate) fn write_rra_summaries(
    prefix: &Path,
    genes: &[AggregatedGeneResult],
    table: &ParsedCountTable,
    control_columns: &[usize],
    treatment_columns: &[usize],
    scores: &[Option<f64>],
) -> Result<(), MageckError> {
    let gene = serialize_gene_summary(genes)?;
    let sgrna = serialize_sgrna_summary(table, control_columns, treatment_columns, scores)?;
    publish_with(&RealFilesystem, prefix, &gene, &sgrna)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;
    use std::io::Cursor;

    fn gene(id: &str, neg: f64, pos: f64) -> AggregatedGeneResult {
        AggregatedGeneResult {
            gene_id: id.into(),
            num_sgrnas: 2,
            num_valid_sgrnas: 2,
            neg_score: neg,
            neg_pvalue: 0.0123456789,
            neg_fdr: 0.2,
            neg_lfc: -1.234567,
            neg_goodsgrna: 1,
            pos_score: pos,
            pos_pvalue: 0.9,
            pos_fdr: 1.0,
            pos_lfc: -1.234567,
            pos_goodsgrna: 0,
        }
    }

    fn table() -> ParsedCountTable {
        ParsedCountTable::from_reader(Cursor::new("sgRNA\tGene\tc1\tc2\tt1\tt2\ng2\tB\t2\t4\t6\t8\ng1\tA\t6\t8\t2\t4\ng3\tB\t0\t0\t0\t0\n")).unwrap()
    }

    fn rows(output: &str) -> Vec<Vec<&str>> {
        output
            .lines()
            .map(|line| line.split('\t').collect())
            .collect()
    }

    #[test]
    fn gene_schema_precision_and_source_order() {
        let output = serialize_gene_summary(&[gene("B", 0.1, 0.9), gene("A", 0.2, 0.1)]).unwrap();
        let rows = rows(&output);
        assert_eq!(
            rows[0],
            [
                "id",
                "num",
                "neg|score",
                "neg|p-value",
                "neg|fdr",
                "neg|rank",
                "neg|lfc",
                "neg|goodsgrna",
                "pos|score",
                "pos|p-value",
                "pos|fdr",
                "pos|rank",
                "pos|lfc",
                "pos|goodsgrna"
            ]
        );
        assert_eq!(
            rows[1],
            [
                "B",
                "2",
                "1.00000000e-1",
                "1.23456789e-2",
                "2.00000000e-1",
                "1",
                "-1.2346",
                "1",
                "9.00000000e-1",
                "9.00000000e-1",
                "1.00000000e0",
                "2",
                "-1.2346",
                "0"
            ]
        );
        assert_eq!(rows[2][0], "A");
        assert!(output.ends_with('\n'));
        assert!(!output.contains('\r'));
    }

    #[test]
    fn independent_ranks_lexical_ties_and_signed_zero() {
        let output = serialize_gene_summary(&[
            gene("Z", 0.0, 0.1),
            gene("A", -0.0, 0.2),
            gene("B", 0.5, 0.1),
        ])
        .unwrap();
        let rows = rows(&output);
        assert_eq!(
            rows.iter()
                .skip(1)
                .map(|r| (r[5], r[11]))
                .collect::<Vec<_>>(),
            [("2", "2"), ("1", "3"), ("3", "1")]
        );
    }

    #[test]
    fn zero_valid_sentinel_passthrough() {
        let mut sentinel = gene("empty", 1.0, 1.0);
        sentinel.num_valid_sgrnas = 0;
        sentinel.neg_pvalue = 1.0;
        sentinel.pos_pvalue = 1.0;
        sentinel.neg_fdr = 1.0;
        sentinel.pos_fdr = 1.0;
        sentinel.neg_goodsgrna = 0;
        let output = serialize_gene_summary(&[sentinel]).unwrap();
        let rows = rows(&output);
        for field in [2, 3, 4, 8, 9, 10] {
            assert_eq!(rows[1][field], "1.00000000e0");
        }
        assert_eq!(
            (rows[1][1], rows[1][6], rows[1][7], rows[1][12], rows[1][13]),
            ("2", "-1.2346", "0", "-1.2346", "0")
        );
    }

    #[test]
    fn empty_gene_payload_has_header_only() {
        assert_eq!(serialize_gene_summary(&[]).unwrap(), GENE_HEADER);
    }

    #[test]
    fn gene_invalid_values_duplicates_and_tsv_identifiers_rejected() {
        let mut g = gene("A", f64::NAN, 1.0);
        assert!(serialize_gene_summary(&[g.clone()]).is_err());
        g.neg_score = 0.1;
        g.pos_lfc = f64::INFINITY;
        assert!(serialize_gene_summary(&[g]).is_err());
        assert!(serialize_gene_summary(&[gene("A", 0.1, 1.0), gene("A", 0.2, 1.0)]).is_err());
        for id in ["", "A\tB", "A\nB", "A\rB"] {
            assert!(serialize_gene_summary(&[gene(id, 0.1, 1.0)]).is_err());
        }
        assert!(serialize_gene_summary(&[gene(" α|A,B ", 0.1, 1.0)])
            .unwrap()
            .contains(" α|A,B \t"));
    }

    #[test]
    fn sgrna_schema_raw_counts_variance_lfc_and_na() {
        let output =
            serialize_sgrna_summary(&table(), &[0, 1], &[2, 3], &[Some(0.25), None, Some(1.0)])
                .unwrap();
        let rows = rows(&output);
        assert_eq!(
            rows[0],
            [
                "sgRNA",
                "Gene",
                "control_count",
                "treatment_count",
                "control_mean",
                "treat_mean",
                "LFC",
                "control_var",
                "adj_var",
                "score",
                "p.low",
                "p.high",
                "p.twosided",
                "FDR",
                "high_in_treatment"
            ]
        );
        assert_eq!(
            rows[1],
            [
                "g2",
                "B",
                "6.00000000e0",
                "1.40000000e1",
                "3.00000000e0",
                "7.00000000e0",
                "1.0000",
                "2.00000000e0",
                "0.00000000e0",
                "2.50000000e-1",
                "NA",
                "NA",
                "NA",
                "NA",
                "True"
            ]
        );
        assert_eq!(rows[2][6], "-1.0000");
        assert_eq!(rows[2][9], "NA");
        assert_eq!(rows[2][14], "False");
        assert_eq!(rows[3][6], "0.0000");
        assert_eq!(rows[3][14], "False");
        assert_eq!(
            rows.iter().skip(1).map(|r| r[0]).collect::<Vec<_>>(),
            ["g2", "g1", "g3"]
        );
        for row in rows.iter().skip(1) {
            assert_eq!(&row[10..14], &["NA"; 4]);
        }
    }

    #[test]
    fn single_control_variance_is_zero() {
        let output = serialize_sgrna_summary(&table(), &[0], &[2], &[None; 3]).unwrap();
        assert_eq!(rows(&output)[1][7], "0.00000000e0");
    }

    #[test]
    fn high_in_treatment_uses_unrounded_lfc() {
        let table =
            ParsedCountTable::from_reader(Cursor::new("sgRNA\tGene\tc\tt\ng\tA\t100000\t100001\n"))
                .unwrap();
        let output = serialize_sgrna_summary(&table, &[0], &[1], &[None]).unwrap();
        assert_eq!(
            (rows(&output)[1][6], rows(&output)[1][14]),
            ("0.0000", "True")
        );
    }

    #[test]
    fn malformed_guide_inputs_rejected() {
        for scores in [
            vec![],
            vec![Some(0.0); 3],
            vec![Some(1.1); 3],
            vec![Some(f64::NAN); 3],
            vec![Some(f64::INFINITY); 3],
        ] {
            assert!(serialize_sgrna_summary(&table(), &[0], &[2], &scores).is_err());
        }
        for (c, t) in [
            (vec![], vec![2]),
            (vec![0, 0], vec![2]),
            (vec![4], vec![2]),
            (vec![0], vec![0]),
            (vec![0], vec![]),
        ] {
            assert!(serialize_sgrna_summary(&table(), &c, &t, &[None; 3]).is_err());
        }
    }

    #[test]
    fn nonfinite_raw_intermediates_rejected() {
        let table = ParsedCountTable::from_reader(Cursor::new(
            "sgRNA\tGene\tc1\tc2\tt\ng\tA\t1e308\t1e308\t1\n",
        ))
        .unwrap();
        assert!(serialize_sgrna_summary(&table, &[0, 1], &[2], &[None]).is_err());
        let table = ParsedCountTable::from_reader(Cursor::new(
            "sgRNA\tGene\tc1\tc2\tt\ng\tA\t0\t1e200\t1\n",
        ))
        .unwrap();
        assert!(serialize_sgrna_summary(&table, &[0, 1], &[2], &[None]).is_err());
    }

    struct Directory(PathBuf);
    impl Directory {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "p4f10-{}-{}",
                std::process::id(),
                TEMP_COUNTER.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }
        fn prefix(&self) -> PathBuf {
            self.0.join("result")
        }
        fn entries(&self) -> usize {
            fs::read_dir(&self.0).unwrap().count()
        }
    }
    impl Drop for Directory {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    struct FaultFilesystem {
        fault: &'static str,
        nth: usize,
        calls: Cell<usize>,
        cleanup_fails: bool,
    }
    impl FaultFilesystem {
        fn fail(&self, stage: &str) -> io::Result<()> {
            if stage == self.fault {
                let call = self.calls.get() + 1;
                self.calls.set(call);
                if call == self.nth {
                    return Err(io::Error::new(
                        io::ErrorKind::PermissionDenied,
                        format!("injected {stage}"),
                    ));
                }
            }
            Ok(())
        }
    }
    impl Filesystem for FaultFilesystem {
        fn create_new(&self, path: &Path) -> io::Result<File> {
            self.fail("create")?;
            RealFilesystem.create_new(path)
        }
        fn write(&self, file: File, contents: &str) -> io::Result<()> {
            self.fail("write")?;
            RealFilesystem.write(file, contents)
        }
        fn rename(&self, from: &Path, to: &Path) -> io::Result<()> {
            self.fail("rename")?;
            // Both payloads must be completely written before the first publication.
            assert!(!fs::read_to_string(from)?.is_empty());
            if self.calls.get() == 1 && self.fault == "rename" {
                assert!(to.ends_with("result.gene_summary.txt"));
            }
            RealFilesystem.rename(from, to)
        }
        fn remove(&self, path: &Path) -> io::Result<()> {
            if self.cleanup_fails {
                Err(io::Error::new(
                    io::ErrorKind::PermissionDenied,
                    "injected cleanup",
                ))
            } else {
                RealFilesystem.remove(path)
            }
        }
    }

    #[test]
    fn writes_real_summaries_and_leaves_only_two_final_files() {
        let dir = Directory::new();
        write_rra_summaries(
            &dir.prefix(),
            &[gene("B", 0.1, 0.9), gene("A", 0.2, 0.1)],
            &table(),
            &[0, 1],
            &[2, 3],
            &[None; 3],
        )
        .unwrap();
        assert_eq!(dir.entries(), 2);
        assert!(fs::read_to_string(dir.0.join("result.gene_summary.txt"))
            .unwrap()
            .starts_with(GENE_HEADER));
        assert!(fs::read_to_string(dir.0.join("result.sgrna_summary.txt"))
            .unwrap()
            .starts_with(SGRNA_HEADER));
    }

    #[test]
    fn rejects_either_existing_output_and_dangling_symlink() {
        for suffix in ["gene_summary.txt", "sgrna_summary.txt"] {
            let dir = Directory::new();
            let output = dir.0.join(format!("result.{suffix}"));
            fs::write(&output, "keep").unwrap();
            assert!(publish_with(&RealFilesystem, &dir.prefix(), "gene", "guide").is_err());
            assert_eq!(fs::read_to_string(output).unwrap(), "keep");
            assert_eq!(dir.entries(), 1);
        }
        let dir = Directory::new();
        std::os::unix::fs::symlink("missing", dir.0.join("result.gene_summary.txt")).unwrap();
        assert!(publish_with(&RealFilesystem, &dir.prefix(), "gene", "guide").is_err());
        assert_eq!(dir.entries(), 1);
    }

    #[test]
    fn create_write_and_both_rename_failures_roll_back() {
        for fault in ["create", "write", "rename"] {
            for nth in 1..=2 {
                let dir = Directory::new();
                let filesystem = FaultFilesystem {
                    fault,
                    nth,
                    calls: Cell::new(0),
                    cleanup_fails: false,
                };
                let error = publish_with(&filesystem, &dir.prefix(), "gene", "guide").unwrap_err();
                assert!(error.to_string().contains(&format!("injected {fault}")));
                assert_eq!(dir.entries(), 0, "{fault} {nth}");
            }
        }
    }

    #[test]
    fn rollback_failures_preserve_primary_error_and_report_cleanup() {
        let dir = Directory::new();
        let filesystem = FaultFilesystem {
            fault: "rename",
            nth: 2,
            calls: Cell::new(0),
            cleanup_fails: true,
        };
        let error = publish_with(&filesystem, &dir.prefix(), "gene", "guide")
            .unwrap_err()
            .to_string();
        assert!(error.contains("injected rename"));
        assert!(error.contains("cleanup failed"));
        assert_eq!(dir.entries(), 2);
    }

    #[test]
    fn temporary_collision_limit_and_sibling_naming() {
        struct Collisions(Cell<usize>);
        impl Filesystem for Collisions {
            fn create_new(&self, path: &Path) -> io::Result<File> {
                self.0.set(self.0.get() + 1);
                assert_eq!(path.parent().unwrap(), Path::new("/tmp"));
                assert!(path
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .contains(&format!(".tmp.{}.", std::process::id())));
                Err(io::Error::from(io::ErrorKind::AlreadyExists))
            }
            fn write(&self, _: File, _: &str) -> io::Result<()> {
                unreachable!()
            }
            fn rename(&self, _: &Path, _: &Path) -> io::Result<()> {
                unreachable!()
            }
            fn remove(&self, _: &Path) -> io::Result<()> {
                unreachable!()
            }
        }
        let filesystem = Collisions(Cell::new(0));
        assert!(create_temp(&filesystem, Path::new("/tmp/result")).is_err());
        assert_eq!(filesystem.0.get(), 32);
        let dir = Directory::new();
        let (a, fa) = create_temp(&RealFilesystem, &dir.prefix()).unwrap();
        let (b, fb) = create_temp(&RealFilesystem, &dir.prefix()).unwrap();
        assert_ne!(a, b);
        drop((fa, fb));
    }

    #[test]
    fn exclusive_rename_preserves_concurrently_created_destination() {
        let dir = Directory::new();
        let source = dir.0.join("temporary");
        let target = dir.0.join("final");
        fs::write(&source, "new").unwrap();
        fs::write(&target, "existing").unwrap();
        assert!(RealFilesystem.rename(&source, &target).is_err());
        assert_eq!(fs::read_to_string(&target).unwrap(), "existing");
        assert_eq!(fs::read_to_string(&source).unwrap(), "new");
    }

    #[test]
    fn formatting_failure_creates_no_outputs() {
        let dir = Directory::new();
        assert!(write_rra_summaries(
            &dir.prefix(),
            &[gene("B", 0.1, 0.9)],
            &table(),
            &[0],
            &[2],
            &[]
        )
        .is_err());
        assert_eq!(dir.entries(), 0);
    }
}
