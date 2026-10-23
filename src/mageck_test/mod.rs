//! Test subcommand module (Phase 4).
//!
//! Implements count-table parsing, sample resolution, and Negative Binomial modeling.

pub(crate) mod aggregation;
pub mod count_table;
pub mod nb_model;
pub mod rra;
pub(crate) mod run;
pub(crate) mod serialize;
pub mod sgrna_pvalue;
pub(crate) mod validity_mask;

pub use count_table::{
    parse_count_table, parse_count_table_file, resolve_sample_indices, CountTable, CountTableParts,
};
pub use nb_model::{
    fit_meanvar, NbModel, RegressionCoefficients, POISSON_DISPERSION, POISSON_SIZE,
};
pub use rra::{
    compute_gene_lo_values, compute_gene_lo_values_from_count_table, ControlResolutionDiagnostics,
    GeneEmpiricalResult, GeneLoValueResult,
};
pub use sgrna_pvalue::{calculate_sgrna_pvalues, SgRnaPvalueResult};
