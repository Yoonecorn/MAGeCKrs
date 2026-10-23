mod cnv_input;
#[allow(dead_code)]
mod correction;
#[allow(dead_code)]
mod design_matrix;
mod diagnostics;
#[allow(dead_code)]
mod fit;
#[allow(dead_code)]
mod gene_result;
#[allow(dead_code)]
mod ingest;
#[allow(
    dead_code,
    unused_assignments,
    clippy::needless_range_loop,
    clippy::too_many_arguments
)]
mod irls;
pub(crate) mod likelihood;
mod meanvar;
#[allow(dead_code)]
mod normalize;
#[allow(dead_code)]
mod output;
pub mod run;
#[allow(dead_code)]
mod serial;
pub use run::run;
#[allow(dead_code, clippy::excessive_precision)]
mod wald;
