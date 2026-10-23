use std::collections::HashSet;
use std::path::{Path, PathBuf};

use crate::cli::{
    CountArgs, MleArgs, NormMethod, PathwayArgs, PathwayMethod, RemoveZeroMode, RunArgs, TestArgs,
};
use crate::error::MageckError;

fn check_file_exists(path: &Path, field_name: &str) -> Result<(), MageckError> {
    if !path.exists() {
        return Err(MageckError::Config(format!(
            "{field_name} file does not exist: {}",
            path.display()
        )));
    }
    if !path.is_file() {
        return Err(MageckError::Config(format!(
            "{field_name} path is not a file: {}",
            path.display()
        )));
    }
    Ok(())
}

fn check_non_empty_str(value: &str, field_name: &str) -> Result<(), MageckError> {
    if value.trim().is_empty() {
        return Err(MageckError::Config(format!(
            "{field_name} must not be empty"
        )));
    }
    Ok(())
}

/// Reject reserved options before input I/O or output publication.
pub(crate) fn validate_count_options(
    count_n: bool,
    _unmapped_to_file: bool,
    threads: Option<crate::cli::ThreadMode>,
) -> Result<(), MageckError> {
    for (supplied, option) in [(count_n, "--count-n")] {
        if supplied {
            return Err(MageckError::Config(format!(
                "Unsupported count option {option}: retired or not implemented; omit this option"
            )));
        }
    }
    validate_test_threads(threads)
}

pub(crate) fn validate_test_threads(
    threads: Option<crate::cli::ThreadMode>,
) -> Result<(), MageckError> {
    if matches!(
        threads,
        Some(crate::cli::ThreadMode::Fixed(0) | crate::cli::ThreadMode::Auto { maximum: Some(0) })
    ) {
        return Err(MageckError::Config(
            "--threads must be greater than zero".into(),
        ));
    }
    Ok(())
}

/// A validated sample group consisting of physical FASTQ paths and an assigned sample label.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqSampleGroup {
    /// Ordered physical FASTQ paths for technical replicates within this sample group.
    pub fastq_paths: Vec<PathBuf>,
    /// Validated sample label for this sample group.
    pub label: String,
}

/// Fixed raw-sequence coordinates and explicit exact-deduplication boundary.
#[derive(Debug, Clone)]
pub struct UmiConfig {
    pub start: usize,
    pub length: usize,
    pub scope: crate::cli::UmiDedupScope,
    pub max_keys: usize,
}

impl UmiConfig {
    pub(crate) fn validate(&self) -> Result<(), MageckError> {
        if !(1..=64).contains(&self.length) || self.max_keys == 0 {
            return Err(MageckError::Config(
                "UMI length must be 1..=64 and --umi-max-keys must be positive".into(),
            ));
        }
        self.start
            .checked_add(self.length)
            .ok_or_else(|| MageckError::Config("UMI coordinate addition overflow".into()))?;
        Ok(())
    }
}

fn umi_config(args: &CountArgs) -> Result<Option<UmiConfig>, MageckError> {
    match (args.umi_start, args.umi_length, args.umi_dedup_scope, args.umi_max_keys) {
        (None, None, None, None) => Ok(None),
        (Some(start), Some(length), Some(scope), Some(max_keys)) => {
            let config = UmiConfig { start, length, scope, max_keys };
            config.validate()?;
            Ok(Some(config))
        }
        _ => Err(MageckError::Config("UMI requires explicit --umi-start, --umi-length, --umi-dedup-scope file|sample and --umi-max-keys together".into())),
    }
}

/// Validated configuration for `mageckrs count`.
#[derive(Debug, Clone)]
pub struct CountConfig {
    pub list_seq: PathBuf,
    pub sample_groups: Vec<FastqSampleGroup>,
    pub output_prefix: String,
    pub norm_method: NormMethod,
    pub control_sgrna: Option<PathBuf>,
    pub sgrna_len: usize,
    pub trim_5: crate::cli::Trim5Mode,
    pub test_mismatch: usize,
    pub count_n: bool,
    pub unmapped_to_file: bool,
    pub umi: Option<UmiConfig>,
    pub threads: Option<crate::cli::ThreadMode>,
}

impl TryFrom<CountArgs> for CountConfig {
    type Error = MageckError;

    fn try_from(args: CountArgs) -> Result<Self, Self::Error> {
        validate_count_options(args.count_n, args.unmapped_to_file, args.threads)?;
        let umi = umi_config(&args)?;
        check_non_empty_str(&args.output_prefix, "Output prefix (-n/--output-prefix)")?;
        check_file_exists(&args.list_seq, "Library (-l/--list-seq)")?;

        if args.fastq.is_empty() {
            return Err(MageckError::Config(
                "At least one FASTQ file (-f/--fastq) is required".to_string(),
            ));
        }

        let n_groups = args.fastq.len();
        let labels = match args.sample_label {
            Some(supplied_labels) => {
                for (idx, label) in supplied_labels.iter().enumerate() {
                    if label.is_empty() || label.trim().is_empty() {
                        return Err(MageckError::Config(format!(
                            "Sample label at index {} must not be empty",
                            idx + 1
                        )));
                    }
                    if label.contains('\t') || label.contains('\r') || label.contains('\n') {
                        return Err(MageckError::Config(format!(
                            "Sample label '{label}' contains invalid control characters (tab, newline, or carriage return)"
                        )));
                    }
                }

                if supplied_labels.len() != n_groups {
                    return Err(MageckError::Config(format!(
                        "Number of sample labels ({}) does not match number of sample groups ({n_groups})",
                        supplied_labels.len()
                    )));
                }

                let mut seen = HashSet::with_capacity(supplied_labels.len());
                for label in &supplied_labels {
                    if !seen.insert(label.clone()) {
                        return Err(MageckError::Config(format!(
                            "Duplicate sample label '{label}' is not allowed"
                        )));
                    }
                }
                supplied_labels
            }
            None => (1..=n_groups).map(|idx| format!("sample{idx}")).collect(),
        };

        let mut raw_groups = Vec::with_capacity(args.fastq.len());
        for (group_idx, group_str) in args.fastq.iter().enumerate() {
            let mut paths = Vec::new();
            for comp in group_str.split(',') {
                if comp.is_empty() {
                    return Err(MageckError::Config(format!(
                        "Empty FASTQ path in sample group {}: '{group_str}'",
                        group_idx + 1
                    )));
                }
                let path = PathBuf::from(comp);
                check_file_exists(&path, "FASTQ (-f/--fastq)")?;
                paths.push(path);
            }
            raw_groups.push(paths);
        }

        let sample_groups = raw_groups
            .into_iter()
            .zip(labels)
            .map(|(fastq_paths, label)| FastqSampleGroup { fastq_paths, label })
            .collect();

        if args.norm_method == NormMethod::Control {
            match &args.control_sgrna {
                Some(ctrl_path) => {
                    check_file_exists(ctrl_path, "Control sgRNA (--control-sgrna)")?;
                }
                None => {
                    return Err(MageckError::Config(
                        "Control normalization (--norm-method control) requires --control-sgrna"
                            .to_string(),
                    ));
                }
            }
        } else if let Some(ref ctrl_path) = args.control_sgrna {
            check_file_exists(ctrl_path, "Control sgRNA (--control-sgrna)")?;
        }

        if args.sgrna_len == 0 || args.sgrna_len > 32 {
            return Err(MageckError::Config(format!(
                "sgRNA length (--sgrna-len) must be between 1 and 32 nt, got {}",
                args.sgrna_len
            )));
        }

        if args.test_mismatch > 1 {
            return Err(MageckError::Config(format!(
                "Mismatch count (--test-mismatch) must be 0 or 1, got {}",
                args.test_mismatch
            )));
        }

        Ok(CountConfig {
            list_seq: args.list_seq,
            sample_groups,
            output_prefix: args.output_prefix,
            norm_method: args.norm_method,
            control_sgrna: args.control_sgrna,
            sgrna_len: args.sgrna_len,
            trim_5: args.trim_5,
            test_mismatch: args.test_mismatch,
            count_n: args.count_n,
            unmapped_to_file: args.unmapped_to_file,
            umi,
            threads: args.threads,
        })
    }
}

/// Validated configuration for `mageckrs test`.
#[derive(Debug, Clone)]
pub struct TestConfig {
    pub count_table: PathBuf,
    pub treatment_id: Vec<String>,
    pub day0_label: Vec<String>,
    pub output_prefix: String,
    pub norm_method: NormMethod,
    pub control_sgrna: Option<PathBuf>,
    pub gene_test_fdr_threshold: f64,
    pub adjust_method: String,
    pub permutation: usize,
    pub remove_zero: RemoveZeroMode,
    pub remove_zero_threshold: f64,
    pub seed: Option<u64>,
    pub rra_seed: u64,
    pub threads: Option<crate::cli::ThreadMode>,
}

impl TryFrom<TestArgs> for TestConfig {
    type Error = MageckError;

    fn try_from(args: TestArgs) -> Result<Self, Self::Error> {
        validate_test_threads(args.threads)?;
        check_non_empty_str(&args.output_prefix, "Output prefix (-n/--output-prefix)")?;
        check_file_exists(&args.count_table, "Count table (-k/--count-table)")?;

        if args.treatment_id.is_empty() {
            return Err(MageckError::Config(
                "Treatment sample ID(s) (-t/--treatment-id) must not be empty".to_string(),
            ));
        }

        if args.day0_label.is_empty() {
            return Err(MageckError::Config(
                "Control sample ID(s) (-c/--day0-label) must not be empty".to_string(),
            ));
        }

        let permutation = args.num_permutations.unwrap_or(args.permutation);
        let rra_seed = args
            .rra_seed
            .or(args.seed)
            .unwrap_or(crate::mageck_test::rra::DEFAULT_RRA_PERMUTATION_SEED);
        if permutation == 0 {
            return Err(MageckError::Config(
                "Number of permutations (--permutation/--num-permutations) must be greater than 0"
                    .to_string(),
            ));
        }

        if !args.remove_zero_threshold.is_finite() || args.remove_zero_threshold < 0.0 {
            return Err(MageckError::Config(format!(
                "Remove-zero threshold (--remove-zero-threshold) must be finite and non-negative, got {}",
                args.remove_zero_threshold
            )));
        }

        if !(0.0..=1.0).contains(&args.gene_test_fdr_threshold) {
            return Err(MageckError::Config(format!(
                "FDR threshold (--gene-test-fdr-threshold) must be between 0.0 and 1.0, got {}",
                args.gene_test_fdr_threshold
            )));
        }

        if args.norm_method == NormMethod::Control {
            match &args.control_sgrna {
                Some(ctrl_path) => {
                    check_file_exists(ctrl_path, "Control sgRNA (--control-sgrna)")?;
                }
                None => {
                    return Err(MageckError::Config(
                        "Control normalization (--norm-method control) requires --control-sgrna"
                            .to_string(),
                    ));
                }
            }
        } else if let Some(ref ctrl_path) = args.control_sgrna {
            check_file_exists(ctrl_path, "Control sgRNA (--control-sgrna)")?;
        }

        Ok(TestConfig {
            count_table: args.count_table,
            treatment_id: args.treatment_id,
            day0_label: args.day0_label,
            output_prefix: args.output_prefix,
            norm_method: args.norm_method,
            control_sgrna: args.control_sgrna,
            gene_test_fdr_threshold: args.gene_test_fdr_threshold,
            adjust_method: args.adjust_method,
            permutation,
            remove_zero: args.remove_zero,
            remove_zero_threshold: args.remove_zero_threshold,
            seed: args.seed,
            rra_seed,
            threads: args.threads,
        })
    }
}

/// Validated configuration for `mageckrs mle`.
#[derive(Debug, Clone)]
pub struct MleConfig {
    pub count_table: PathBuf,
    pub design_matrix: PathBuf,
    pub output_prefix: String,
    pub norm_method: NormMethod,
    pub max_sgrnas_per_gene: usize,
    pub genes_varmodeling: usize,
    pub control_sgrna: Option<PathBuf>,
    pub cnv_norm: Option<PathBuf>,
    pub threads: Option<crate::cli::ThreadMode>,
    pub seed: Option<u64>,
}

impl TryFrom<MleArgs> for MleConfig {
    type Error = MageckError;

    fn try_from(args: MleArgs) -> Result<Self, Self::Error> {
        validate_test_threads(args.threads)?;
        check_non_empty_str(&args.output_prefix, "Output prefix (-n/--output-prefix)")?;
        check_file_exists(&args.count_table, "Count table (-k/--count-table)")?;
        check_file_exists(&args.design_matrix, "Design matrix (-d/--design-matrix)")?;

        if args.norm_method == NormMethod::Control {
            match &args.control_sgrna {
                Some(ctrl_path) => {
                    check_file_exists(ctrl_path, "Control sgRNA (--control-sgrna)")?;
                }
                None => {
                    return Err(MageckError::Config(
                        "Control normalization (--norm-method control) requires --control-sgrna"
                            .to_string(),
                    ));
                }
            }
        } else if let Some(ref ctrl_path) = args.control_sgrna {
            check_file_exists(ctrl_path, "Control sgRNA (--control-sgrna)")?;
        }

        if let Some(ref cnv_path) = args.cnv_norm {
            check_file_exists(cnv_path, "CNV normalization profile (--cnv-norm)")?;
        }

        if args.max_sgrnas_per_gene == 0 {
            return Err(MageckError::Config(
                "--max-sgrnas-per-gene must be >= 1".to_owned(),
            ));
        }

        Ok(MleConfig {
            count_table: args.count_table,
            design_matrix: args.design_matrix,
            output_prefix: args.output_prefix,
            norm_method: args.norm_method,
            max_sgrnas_per_gene: args.max_sgrnas_per_gene,
            genes_varmodeling: args.genes_varmodeling,
            control_sgrna: args.control_sgrna,
            cnv_norm: args.cnv_norm,
            threads: args.threads,
            seed: args.seed,
        })
    }
}

/// Validated configuration for `mageckrs pathway`.
#[derive(Debug, Clone)]
pub struct PathwayConfig {
    pub gene_summary: PathBuf,
    pub gmt_file: PathBuf,
    pub output_prefix: String,
    pub method: PathwayMethod,
    pub permutation: usize,
    pub seed: Option<u64>,
    pub threads: Option<usize>,
}

impl TryFrom<PathwayArgs> for PathwayConfig {
    type Error = MageckError;

    fn try_from(args: PathwayArgs) -> Result<Self, Self::Error> {
        check_non_empty_str(&args.output_prefix, "Output prefix (-n/--output-prefix)")?;
        check_file_exists(&args.gene_summary, "Gene summary (-g/--gene-summary)")?;
        check_file_exists(&args.gmt_file, "GMT file (--gmt-file)")?;

        if args.permutation == 0 {
            return Err(MageckError::Config(
                "Number of permutations (--permutation) must be greater than 0".to_string(),
            ));
        }

        Ok(PathwayConfig {
            gene_summary: args.gene_summary,
            gmt_file: args.gmt_file,
            output_prefix: args.output_prefix,
            method: args.method,
            permutation: args.permutation,
            seed: args.seed,
            threads: args.threads,
        })
    }
}

/// Validated configuration for `mageckrs run` (convenience chaining count → test).
#[derive(Debug, Clone)]
pub struct RunConfig {
    pub count: CountConfig,
    pub output_prefix: String,
    pub treatment_id: Vec<String>,
    pub day0_label: Vec<String>,
    pub permutation: usize,
    pub seed: Option<u64>,
    pub threads: Option<crate::cli::ThreadMode>,
}

impl TryFrom<RunArgs> for RunConfig {
    type Error = MageckError;
    fn try_from(args: RunArgs) -> Result<Self, Self::Error> {
        validate_test_threads(args.threads)?;
        let count = CountConfig::try_from(CountArgs {
            list_seq: args.list_seq,
            fastq: args.fastq,
            output_prefix: args.output_prefix.clone(),
            sample_label: args.sample_label,
            norm_method: args.norm_method,
            control_sgrna: args.control_sgrna,
            sgrna_len: args.sgrna_len,
            trim_5: args.trim_5,
            test_mismatch: args.test_mismatch,
            count_n: args.count_n,
            unmapped_to_file: args.unmapped_to_file,
            threads: args.threads.filter(|mode| mode.is_auto()),
            umi_start: args.umi_start,
            umi_length: args.umi_length,
            umi_dedup_scope: args.umi_dedup_scope,
            umi_max_keys: args.umi_max_keys,
        })?;
        if args.treatment_id.is_empty() || args.day0_label.is_empty() {
            return Err(MageckError::Config(
                "Treatment and control sample IDs must not be empty".into(),
            ));
        }
        if args.permutation == 0 {
            return Err(MageckError::Config(
                "Number of permutations must be greater than 0".into(),
            ));
        }
        Ok(Self {
            count,
            output_prefix: args.output_prefix,
            treatment_id: args.treatment_id,
            day0_label: args.day0_label,
            permutation: args.permutation,
            seed: args.seed,
            threads: args.threads,
        })
    }
}
