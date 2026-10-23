pub mod cli;
pub mod config;
pub mod count;
pub mod data;
pub mod error;
pub mod mageck_test;
pub(crate) mod mle;
#[allow(dead_code)] // PRNG1 core is registered before its later in-crate consumers.
pub(crate) mod prng;
pub(crate) mod run;
pub mod stats;

use clap::Parser;

use crate::cli::{Cli, Commands};
use crate::config::{CountConfig, MleConfig, TestConfig};
use crate::error::MageckError;

/// Main library entry point: parses CLI arguments, validates configuration,
/// and dispatches to the appropriate subcommand pipeline.
pub fn run() -> Result<(), MageckError> {
    let cli = Cli::parse();
    run_with_cli(cli)
}

/// Dispatches a parsed CLI invocation. Separated from `run()` to facilitate testing.
pub fn run_with_cli(cli: Cli) -> Result<(), MageckError> {
    match cli.command {
        Commands::Count(args) => {
            let config = CountConfig::try_from(args)?;
            crate::count::run::run(&config)?;
            Ok(())
        }
        Commands::Test(args) => {
            let config = TestConfig::try_from(args)?;
            crate::mageck_test::run::run(&config)
        }
        Commands::Mle(args) => {
            let config = MleConfig::try_from(args)?;
            crate::mle::run(&config)?;
            Ok(())
        }
        Commands::Pathway(_) => Err(MageckError::Config(
            "Unsupported command 'pathway': pathway enrichment is not implemented; no analysis was run"
                .into(),
        )),
        Commands::Run(args) => {
            let config = crate::config::RunConfig::try_from(args)?;
            crate::run::run(&config)
        },
    }
}

mod execution;
