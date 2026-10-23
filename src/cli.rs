use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

/// Command-line parser for MAGeCKrs.
#[derive(Parser, Debug, Clone)]
#[command(
    name = "mageckrs",
    author,
    version,
    about = "CRISPR screen analysis for macOS Apple Silicon: count, test, combined run, and classical mle",
    long_about = "MAGeCKrs: CRISPR screen analysis for macOS Apple Silicon. COUNT, explicitly scoped exact UMI counting, COUNT-to-test run, RRA/test and classical MLE are supported. MLE, test and run accept deterministic CPU worker limits; COUNT remains serial under every worker policy. Auto is opt-in; use fixed workers for publication. RRA p-values and FDRs are nominal rank-model diagnostics: interpret with rank, effect direction, guide consistency, reproducibility and biological validation. Universal biological calibration is not claimed. MLE permutation-style columns are placeholders; CNV corrects matched betas without recomputing uncertainty or significance. Pathway is reserved and returns an error. Metal execution and static distribution are deferred.",
    disable_help_subcommand = true
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

/// Available subcommands in MAGeCKrs.
#[derive(Subcommand, Debug, Clone)]
pub enum Commands {
    /// Count sgRNAs from FASTQ files using a reference library
    Count(CountArgs),
    /// Test for gene essentiality/enrichment using Robust Rank Aggregation (RRA)
    Test(TestArgs),
    /// Perform Maximum Likelihood Estimation (MLE) using a design matrix
    Mle(MleArgs),
    /// Unsupported: reserved for future pathway enrichment (returns an error)
    #[command(
        after_help = "All pathway options are reserved. This command returns a configuration error before reading inputs or writing outputs."
    )]
    Pathway(PathwayArgs),
    /// Count FASTQ samples then test their raw counts with RRA
    Run(RunArgs),
}

/// Exact UMI deduplication boundary, explicitly selected when enabled.
#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiDedupScope {
    File,
    Sample,
}

/// Normalization methods supported across counting and testing.
#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NormMethod {
    #[default]
    Median,
    Total,
    Control,
    None,
}

impl std::fmt::Display for NormMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NormMethod::Median => write!(f, "median"),
            NormMethod::Total => write!(f, "total"),
            NormMethod::Control => write!(f, "control"),
            NormMethod::None => write!(f, "none"),
        }
    }
}

/// Guide-removal policy based on normalized comparison means.
#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RemoveZeroMode {
    None,
    Control,
    Treatment,
    Any,
    #[default]
    Both,
}

impl std::fmt::Display for RemoveZeroMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RemoveZeroMode::None => write!(f, "none"),
            RemoveZeroMode::Control => write!(f, "control"),
            RemoveZeroMode::Treatment => write!(f, "treatment"),
            RemoveZeroMode::Any => write!(f, "any"),
            RemoveZeroMode::Both => write!(f, "both"),
        }
    }
}

/// Default normalized-mean threshold for guide removal.
pub const DEFAULT_REMOVE_ZERO_THRESHOLD: f64 = 0.0;

/// Guide start selection for count. Auto scans offsets 0..=40 per read for the first exact match.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trim5Mode {
    None,
    Auto,
    Fixed(usize),
}

impl std::str::FromStr for Trim5Mode {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "none" => Ok(Self::None),
            "auto" => Ok(Self::Auto),
            _ => value
                .parse::<usize>()
                .map(Self::Fixed)
                .map_err(|_| "--trim-5 must be auto, none, or a non-negative integer".to_string()),
        }
    }
}

pub use crate::execution::ThreadMode;

/// Command-line arguments for the `count` subcommand.
#[derive(Args, Debug, Clone)]
pub struct CountArgs {
    /// Library file (TSV/CSV) containing sgRNA ID, sequence, and gene name
    #[arg(short = 'l', long = "list-seq", visible_alias = "library")]
    pub list_seq: PathBuf,

    /// One or more FASTQ sample groups (comma-separated technical replicates within each argument)
    #[arg(short = 'f', long = "fastq", num_args = 1.., required = true)]
    pub fastq: Vec<String>,

    /// Prefix for output files
    #[arg(short = 'n', long = "output-prefix")]
    pub output_prefix: String,

    /// Sample labels corresponding to the sample groups (comma-separated)
    #[arg(long = "sample-label", value_delimiter = ',')]
    pub sample_label: Option<Vec<String>>,

    /// Normalization method (median, total, control, none)
    #[arg(long = "norm-method", value_enum, default_value_t = NormMethod::Median)]
    pub norm_method: NormMethod,

    /// Optional control sgRNA list file for control normalization
    #[arg(long = "control-sgrna")]
    pub control_sgrna: Option<PathBuf>,

    /// Length of sgRNA sequence (default: 20)
    #[arg(long = "sgrna-len", default_value_t = 20)]
    pub sgrna_len: usize,

    /// Guide start: auto (first exact match at offsets 0..=40), none (offset 0), or integer
    #[arg(long = "trim-5", default_value = "none", value_name = "MODE")]
    pub trim_5: Trim5Mode,

    /// Mismatches: 0 or 1 at fixed/none offsets; auto always uses exact matches
    #[arg(long = "test-mismatch", default_value_t = 0)]
    pub test_mismatch: usize,

    /// Unsupported: --count-n is retired; wildcard/ambiguous-base counting is not supported
    #[arg(long = "count-n")]
    pub count_n: bool,

    /// Write no-match, ambiguous and short-guide records to `<prefix>.unmapped.fastq`
    #[arg(long = "unmapped-to-file")]
    pub unmapped_to_file: bool,

    /// Absolute zero-based UMI start in the original FASTQ sequence (requires all UMI options)
    #[arg(long = "umi-start")]
    pub umi_start: Option<usize>,
    /// Fixed UMI length, 1..=64; exact A/C/G/T equality, no clustering
    #[arg(long = "umi-length")]
    pub umi_length: Option<usize>,
    /// Required UMI scope: file resets per physical file; sample spans grouped files
    #[arg(long = "umi-dedup-scope", value_enum)]
    pub umi_dedup_scope: Option<UmiDedupScope>,
    /// Required positive unique-key cap per active scope; not a process-memory bound
    #[arg(long = "umi-max-keys")]
    pub umi_max_keys: Option<usize>,

    /// Worker policy N, auto, or auto:MAX; COUNT is qualified only for serial execution
    #[arg(long = "threads")]
    pub threads: Option<ThreadMode>,
}

/// Command-line arguments for the `test` subcommand.
#[derive(Args, Debug, Clone)]
pub struct TestArgs {
    /// sgRNA count table file
    #[arg(short = 'k', long = "count-table")]
    pub count_table: PathBuf,

    /// Treatment sample label(s) or 0-based column indices (comma-separated)
    #[arg(short = 't', long = "treatment-id", value_delimiter = ',')]
    pub treatment_id: Vec<String>,

    /// Control / Day0 sample label(s) or 0-based column indices (comma-separated)
    #[arg(
        short = 'c',
        long = "day0-label",
        visible_alias = "control-id",
        value_delimiter = ','
    )]
    pub day0_label: Vec<String>,

    /// Prefix for output files
    #[arg(short = 'n', long = "output-prefix")]
    pub output_prefix: String,

    /// Normalization method (median, total, control, none)
    #[arg(long = "norm-method", value_enum, default_value_t = NormMethod::Median)]
    pub norm_method: NormMethod,

    /// Optional control sgRNA list file
    #[arg(long = "control-sgrna")]
    pub control_sgrna: Option<PathBuf>,

    /// Guide p-value threshold for directional RRA alpha selection; does not filter output
    #[arg(long = "gene-test-fdr-threshold", default_value_t = 0.25)]
    pub gene_test_fdr_threshold: f64,

    /// P-value adjustment: fdr or bh (both apply Benjamini-Hochberg)
    #[arg(long = "adjust-method", default_value = "fdr")]
    pub adjust_method: String,

    /// Number of permutations for RRA FDR calculation
    #[arg(long = "permutation", default_value_t = 1000)]
    pub permutation: usize,

    /// Number of RRA permutations (default: 1000); alternative to --permutation
    #[arg(long = "num-permutations", conflicts_with = "permutation")]
    pub num_permutations: Option<usize>,

    /// RRA permutation seed (default: 0x0123456789abcdef = 81985529216486895; enter decimal u64)
    #[arg(long = "rra-seed", conflicts_with = "seed")]
    pub rra_seed: Option<u64>,

    /// Remove guides by normalized means: any requires both groups to pass; both requires either
    #[arg(long = "remove-zero", value_enum, default_value_t = RemoveZeroMode::Both)]
    pub remove_zero: RemoveZeroMode,

    /// Strict normalized-mean threshold used by --remove-zero
    #[arg(
        long = "remove-zero-threshold",
        default_value_t = DEFAULT_REMOVE_ZERO_THRESHOLD
    )]
    pub remove_zero_threshold: f64,

    /// Alternative to --rra-seed; fixed default 81985529216486895 (decimal u64 input)
    #[arg(long = "seed")]
    pub seed: Option<u64>,

    /// Worker limit N, auto, or auto:MAX; omitted or 1 is serial; auto is opt-in
    #[arg(long = "threads")]
    pub threads: Option<ThreadMode>,
}

/// Command-line arguments for the `mle` subcommand.
#[derive(Args, Debug, Clone)]
pub struct MleArgs {
    /// sgRNA count table file
    #[arg(short = 'k', long = "count-table")]
    pub count_table: PathBuf,

    /// Design matrix file (TSV)
    #[arg(short = 'd', long = "design-matrix")]
    pub design_matrix: PathBuf,

    /// Prefix for output files
    #[arg(short = 'n', long = "output-prefix")]
    pub output_prefix: String,

    /// Normalization method (median, total, control, none)
    #[arg(long = "norm-method", value_enum, default_value_t = NormMethod::Median)]
    pub norm_method: NormMethod,

    /// Skip fitting genes with more guides than this limit (must be >= 1)
    #[arg(long = "max-sgrnas-per-gene", default_value_t = 500)]
    pub max_sgrnas_per_gene: usize,

    /// Number of genes used to fit the mean-variance dispersion model
    #[arg(
        long = "genes-varmodeling",
        default_value_t = crate::mle::run::DEFAULT_GENES_VARMODELING,
        help = "Number of genes used to fit the mean-variance dispersion model (default: 0, scalar alpha=0.01). Set to 2000 to match the MAGeCK standalone MLE entry point behavior."
    )]
    pub genes_varmodeling: usize,

    /// Control guide IDs for control normalization; unused with other normalization methods
    #[arg(long = "control-sgrna")]
    pub control_sgrna: Option<PathBuf>,

    /// Log2 CNV TSV: SYMBOL and exact design-condition profiles; correct betas after MLE
    #[arg(long = "cnv-norm")]
    pub cnv_norm: Option<PathBuf>,

    /// Worker limit N, auto, or auto:MAX; omitted or 1 is serial; auto is opt-in
    #[arg(long = "threads")]
    pub threads: Option<ThreadMode>,

    /// Unsupported: MLE permutations are deferred (supplying this option errors)
    #[arg(long = "seed")]
    pub seed: Option<u64>,
}

/// Pathway analysis method options.
#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PathwayMethod {
    #[default]
    Both,
    Rra,
    Gsea,
}

/// Command-line arguments for the `pathway` subcommand.
#[derive(Args, Debug, Clone)]
pub struct PathwayArgs {
    /// Gene summary file (.gene_summary.txt) from mageckrs test or mle
    #[arg(short = 'g', long = "gene-summary")]
    pub gene_summary: PathBuf,

    /// GMT format pathway database file
    #[arg(long = "gmt-file")]
    pub gmt_file: PathBuf,

    /// Prefix for output files
    #[arg(short = 'n', long = "output-prefix")]
    pub output_prefix: String,

    /// Pathway enrichment method (rra, gsea, both)
    #[arg(long = "method", value_enum, default_value_t = PathwayMethod::Both)]
    pub method: PathwayMethod,

    /// Number of permutations for pathway statistical testing
    #[arg(long = "permutation", default_value_t = 1000)]
    pub permutation: usize,

    /// Reserved seed; this command is unsupported
    #[arg(long = "seed")]
    pub seed: Option<u64>,

    /// Reserved worker count; this command is unsupported
    #[arg(long = "threads")]
    pub threads: Option<usize>,
}

/// Command-line arguments for the `run` convenience subcommand (count + test).
#[derive(Args, Debug, Clone)]
pub struct RunArgs {
    /// Library file (TSV/CSV) containing sgRNA ID, sequence, and gene name
    #[arg(short = 'l', long = "list-seq", visible_alias = "library")]
    pub list_seq: PathBuf,

    /// FASTQ sample groups; comma-separated technical replicates within each argument
    #[arg(short = 'f', long = "fastq", num_args = 1.., required = true)]
    pub fastq: Vec<String>,

    /// Treatment sample label(s) or indices (comma-separated)
    #[arg(short = 't', long = "treatment-id", value_delimiter = ',')]
    pub treatment_id: Vec<String>,

    /// Control / Day0 sample label(s) or indices (comma-separated)
    #[arg(
        short = 'c',
        long = "day0-label",
        visible_alias = "control-id",
        value_delimiter = ','
    )]
    pub day0_label: Vec<String>,

    /// Prefix for output files
    #[arg(short = 'n', long = "output-prefix")]
    pub output_prefix: String,

    /// Labels corresponding to FASTQ sample groups (comma-separated)
    #[arg(long = "sample-label", value_delimiter = ',')]
    pub sample_label: Option<Vec<String>>,

    /// Normalization method (median, total, control, none)
    #[arg(long = "norm-method", value_enum, default_value_t = NormMethod::Median)]
    pub norm_method: NormMethod,

    /// Optional control sgRNA list file
    #[arg(long = "control-sgrna")]
    pub control_sgrna: Option<PathBuf>,

    /// Length of sgRNA sequence (default: 20)
    #[arg(long = "sgrna-len", default_value_t = 20)]
    pub sgrna_len: usize,

    /// Guide start: auto (first exact offset 0..=40), none (0), or integer
    #[arg(long = "trim-5", default_value = "none", value_name = "MODE")]
    pub trim_5: Trim5Mode,

    /// Mismatches: 0 or 1 at fixed/none offsets; auto remains exact-only
    #[arg(long = "test-mismatch", default_value_t = 0)]
    pub test_mismatch: usize,

    /// Unsupported: --count-n is retired; wildcard/ambiguous-base counting is not supported
    #[arg(long = "count-n")]
    pub count_n: bool,

    /// Write unmapped FASTQ records with protected publication
    #[arg(long = "unmapped-to-file")]
    pub unmapped_to_file: bool,

    /// Absolute zero-based UMI start in the original FASTQ sequence (requires all UMI options)
    #[arg(long = "umi-start")]
    pub umi_start: Option<usize>,
    /// Fixed UMI length, 1..=64; exact A/C/G/T equality, no clustering
    #[arg(long = "umi-length")]
    pub umi_length: Option<usize>,
    /// Required UMI scope: file resets per physical file; sample spans grouped files
    #[arg(long = "umi-dedup-scope", value_enum)]
    pub umi_dedup_scope: Option<UmiDedupScope>,
    /// Required positive unique-key cap per active scope; not a process-memory bound
    #[arg(long = "umi-max-keys")]
    pub umi_max_keys: Option<usize>,

    /// Number of permutations for RRA FDR calculation
    #[arg(long = "permutation", default_value_t = 1000)]
    pub permutation: usize,

    /// RRA seed (decimal u64); omitted uses the fixed existing test default
    #[arg(long = "seed")]
    pub seed: Option<u64>,

    /// RRA worker limit N, auto, or auto:MAX; counting remains serial
    #[arg(long = "threads")]
    pub threads: Option<ThreadMode>,
}
