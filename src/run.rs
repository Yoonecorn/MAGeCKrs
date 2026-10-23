//! Shared staged publication for opt-in COUNT and combined COUNT-to-test.
use crate::error::MageckError;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static STAGE_COUNTER: AtomicU64 = AtomicU64::new(0);

pub(crate) struct OutputStage {
    directory: PathBuf,
    pub prefix: String,
}

impl OutputStage {
    pub fn new(prefix: &str, suffixes: &[&str]) -> Result<Self, MageckError> {
        if prefix.trim().is_empty() || prefix.contains(['\t', '\r', '\n', '\0']) {
            return Err(MageckError::Config(
                "Output prefix is empty or contains control characters".into(),
            ));
        }
        for suffix in suffixes {
            reject_existing(&PathBuf::from(format!("{prefix}.{suffix}")))?;
        }
        let parent = Path::new(prefix)
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        for _ in 0..32 {
            let directory = parent.join(format!(
                ".mageckrs-stage-{}-{}",
                std::process::id(),
                STAGE_COUNTER.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::create_dir(&directory) {
                Ok(()) => {
                    let staged = directory.join("out");
                    return Ok(Self {
                        prefix: staged
                            .to_str()
                            .ok_or_else(|| MageckError::Config("Output path must be UTF-8".into()))?
                            .into(),
                        directory,
                    });
                }
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e.into()),
            }
        }
        Err(MageckError::Config(
            "Output staging directory allocation exhausted".into(),
        ))
    }

    pub fn publish(&self, destination: &str, suffixes: &[&str]) -> Result<(), MageckError> {
        self.publish_with(destination, suffixes, |_, _| Ok(()))
    }

    fn publish_with<F>(
        &self,
        destination: &str,
        suffixes: &[&str],
        mut before: F,
    ) -> Result<(), MageckError>
    where
        F: FnMut(usize, &Path) -> Result<(), MageckError>,
    {
        let mut published = Vec::new();
        for (i, suffix) in suffixes.iter().enumerate() {
            let from = PathBuf::from(format!("{}.{}", self.prefix, suffix));
            let to = PathBuf::from(format!("{destination}.{suffix}"));
            let result = before(i, &to).and_then(|()| {
                crate::mageck_test::run::exclusive_rename(&from, &to).map_err(MageckError::Io)
            });
            if let Err(error) = result {
                let mut details = format!("Publishing '{}': {error}", to.display());
                for owned in published {
                    if let Err(cleanup) = fs::remove_file(&owned) {
                        details.push_str(&format!(
                            "; rollback failed for {}: {cleanup}",
                            Path::new(&owned).display()
                        ));
                    }
                }
                return Err(MageckError::Io(std::io::Error::other(details)));
            }
            published.push(to);
        }
        Ok(())
    }
}

impl Drop for OutputStage {
    fn drop(&mut self) {
        if let Err(e) = fs::remove_dir_all(&self.directory) {
            log::warn!(
                "Staging cleanup failed for {}: {e}",
                self.directory.display()
            );
        }
    }
}

fn reject_existing(path: &Path) -> Result<(), MageckError> {
    match fs::symlink_metadata(path) {
        Ok(_) => Err(MageckError::Config(format!(
            "Output already exists: {}",
            path.display()
        ))),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(e.into()),
    }
}

/// Combined run derives each stage's normalization from the same raw neutral counts.
pub(crate) fn run(config: &crate::config::RunConfig) -> Result<(), MageckError> {
    use std::io::Write;
    crate::config::validate_test_threads(config.threads)?;
    let mut suffixes = crate::count::run::output_suffixes(&config.count);
    suffixes.retain(|s| *s != "run_info.txt");
    suffixes.extend([
        "count_run_info.txt",
        "gene_summary.txt",
        "sgrna_summary.txt",
        "run_info.txt",
    ]);
    let stage = OutputStage::new(&config.output_prefix, &suffixes)?;
    let mut count = config.count.clone();
    count.output_prefix = stage.prefix.clone();
    let mut unmapped = if count.unmapped_to_file {
        Some(std::io::BufWriter::new(fs::File::create(format!(
            "{}.unmapped.fastq",
            stage.prefix
        ))?))
    } else {
        None
    };
    let prepared =
        crate::count::run::prepare(&count, unmapped.as_mut().map(|w| w as &mut dyn Write))?;
    if let Some(writer) = unmapped.as_mut() {
        writer.flush()?;
    }
    let raw = crate::data::ParsedCountTable::from_parts(
        prepared
            .library
            .records()
            .iter()
            .map(|r| r.id.clone())
            .collect(),
        prepared
            .library
            .records()
            .iter()
            .map(|r| r.gene.clone())
            .collect(),
        count
            .sample_groups
            .iter()
            .map(|g| g.label.clone())
            .collect(),
        prepared.grouped.counts.clone(),
    )?;
    crate::count::run::write_prepared(&count, &prepared)?;
    fs::rename(
        format!("{}.run_info.txt", stage.prefix),
        format!("{}.count_run_info.txt", stage.prefix),
    )?;
    let test = crate::config::TestConfig {
        count_table: format!("{}.count.txt", config.output_prefix).into(),
        treatment_id: config.treatment_id.clone(),
        day0_label: config.day0_label.clone(),
        output_prefix: stage.prefix.clone(),
        norm_method: count.norm_method,
        control_sgrna: count.control_sgrna.clone(),
        gene_test_fdr_threshold: 0.25,
        adjust_method: "fdr".into(),
        permutation: config.permutation,
        remove_zero: crate::cli::RemoveZeroMode::Both,
        remove_zero_threshold: crate::cli::DEFAULT_REMOVE_ZERO_THRESHOLD,
        seed: config.seed,
        rra_seed: config
            .seed
            .unwrap_or(crate::mageck_test::rra::DEFAULT_RRA_PERMUTATION_SEED),
        threads: config.threads,
    };
    crate::mageck_test::run::run_from_table(&test, raw, &config.output_prefix)?;
    stage.publish(&config.output_prefix, &suffixes)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn late_conflict_never_overwrites_and_rolls_back_owned_files() {
        let d = tempfile::tempdir().unwrap();
        let prefix = d.path().join("result").display().to_string();
        let s = OutputStage::new(&prefix, &["a", "b"]).unwrap();
        fs::write(format!("{}.a", s.prefix), "our a").unwrap();
        fs::write(format!("{}.b", s.prefix), "our b").unwrap();
        let error = s
            .publish_with(&prefix, &["a", "b"], |i, path| {
                if i == 1 {
                    fs::write(path, "concurrent sentinel")?;
                }
                Ok(())
            })
            .unwrap_err();
        assert!(error.to_string().contains("Publishing"));
        assert!(!Path::new(&format!("{prefix}.a")).exists());
        assert_eq!(
            fs::read_to_string(format!("{prefix}.b")).unwrap(),
            "concurrent sentinel"
        );
    }
    #[test]
    fn missing_staged_file_rolls_back_and_cleanup_failure_is_reported() {
        let d = tempfile::tempdir().unwrap();
        let prefix = d.path().join("result").display().to_string();
        let s = OutputStage::new(&prefix, &["a", "b"]).unwrap();
        fs::write(format!("{}.a", s.prefix), "our a").unwrap();
        let error = s.publish(&prefix, &["a", "b"]).unwrap_err();
        assert!(error.to_string().contains("Publishing"));
        assert!(!Path::new(&format!("{prefix}.a")).exists());
        fs::write(format!("{}.a", s.prefix), "our a").unwrap();
        let error = s
            .publish_with(&prefix, &["a", "b"], |i, _| {
                if i == 1 {
                    fs::remove_file(format!("{prefix}.a"))?;
                    fs::create_dir(format!("{prefix}.a"))?;
                }
                Ok(())
            })
            .unwrap_err();
        assert!(error.to_string().contains("rollback failed"));
        fs::remove_dir(format!("{prefix}.a")).unwrap();
    }

    #[test]
    fn first_late_conflict_preserves_symlinks_and_cleans_staging() {
        use std::os::unix::fs::symlink;
        let d = tempfile::tempdir().unwrap();
        let prefix = d.path().join("result").display().to_string();
        let s = OutputStage::new(&prefix, &["a"]).unwrap();
        let directory = s.directory.clone();
        fs::write(format!("{}.a", s.prefix), "ours").unwrap();
        assert!(s
            .publish_with(&prefix, &["a"], |_, path| {
                symlink("missing", path)?;
                Ok(())
            })
            .is_err());
        assert_eq!(
            fs::read_link(format!("{prefix}.a")).unwrap(),
            PathBuf::from("missing")
        );
        drop(s);
        assert!(!directory.exists());
    }
}
