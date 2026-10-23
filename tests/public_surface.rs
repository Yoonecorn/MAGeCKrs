//! M10 process-boundary checks: reserved features fail cleanly without I/O,
//! while the retained supported defaults remain usable.
use clap::{CommandFactory, Parser};
use mageckrs::{cli::Cli, config::CountConfig, error::MageckError};
use std::{fs, process::Command};

fn invoke(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_mageckrs"))
        .args(args)
        .output()
        .unwrap()
}

fn assert_unsupported(args: &[&str], needle: &str) {
    let output = invoke(args);
    assert_eq!(output.status.code(), Some(1), "{args:?}: {output:?}");
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.starts_with("Configuration error:"), "{stderr}");
    assert!(stderr.contains(needle), "{stderr}");
    assert!(!stderr.contains("panicked"), "{stderr}");
    assert!(output.stdout.is_empty());
}

#[test]
fn reserved_commands_reject_before_input_io_or_output_creation() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("nested/out");
    let prefix = prefix.to_str().unwrap();
    for method in ["rra", "gsea", "both"] {
        assert_unsupported(
            &[
                "pathway",
                "-g",
                "missing.genes",
                "--gmt-file",
                "missing.gmt",
                "-n",
                prefix,
                "--method",
                method,
                "--threads",
                "2",
                "--seed",
                "7",
            ],
            "Unsupported command 'pathway'",
        );
    }
    assert_unsupported(
        &[
            "run",
            "-l",
            "missing.lib",
            "-f",
            "missing.fq",
            "-t",
            "T",
            "-c",
            "C",
            "-n",
            prefix,
            "--count-n",
            "--unmapped-to-file",
            "--threads",
            "2",
        ],
        "--count-n",
    );
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
}

#[test]
fn reserved_commands_with_existing_files_also_return_errors() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input");
    fs::write(&input, "not parsed").unwrap();
    let input = input.to_str().unwrap();
    let prefix = dir.path().join("out");
    let prefix = prefix.to_str().unwrap();
    assert_unsupported(
        &["pathway", "-g", input, "--gmt-file", input, "-n", prefix],
        "Unsupported command 'pathway'",
    );
    let output = invoke(&[
        "run", "-l", input, "-f", input, "-t", "T", "-c", "C", "-n", prefix,
    ]);
    assert_eq!(output.status.code(), Some(1));
    let error = String::from_utf8(output.stderr).unwrap();
    assert!(error.contains("Library"), "{error}");
    assert!(!error.contains("Unsupported command 'run'"));
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
}

#[test]
fn count_retired_option_and_zero_workers_are_errors() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    let prefix = prefix.to_str().unwrap();
    for extra in [vec!["--count-n"], vec!["--threads", "0"]] {
        let mut args = vec![
            "count",
            "-l",
            "missing.lib",
            "-f",
            "missing.fq",
            "-n",
            prefix,
        ];
        args.extend(extra.iter().copied());
        if extra[0] == "--count-n" {
            assert_unsupported(&args, extra[0]);
        } else {
            let output = invoke(&args);
            assert_eq!(output.status.code(), Some(1));
            assert!(String::from_utf8(output.stderr)
                .unwrap()
                .contains("greater than zero"));
        }
    }
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
}

#[test]
fn public_count_config_cannot_bypass_inactive_option_rejection() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input");
    fs::write(&input, "").unwrap();
    let cli = Cli::try_parse_from([
        "mageckrs",
        "count",
        "-l",
        input.to_str().unwrap(),
        "-f",
        input.to_str().unwrap(),
        "-n",
        "unused",
    ])
    .unwrap();
    let mageckrs::cli::Commands::Count(args) = cli.command else {
        unreachable!()
    };
    let config = CountConfig::try_from(args).unwrap();
    for option in ["--count-n", "--threads"] {
        let mut config = config.clone();
        match option {
            "--count-n" => config.count_n = true,
            _ => config.threads = Some(mageckrs::cli::ThreadMode::Fixed(0)),
        }
        let error = mageckrs::count::run(&config).unwrap_err();
        assert!(matches!(error, MageckError::Config(_)));
        assert!(error.to_string().contains(option));
    }
}

#[test]
fn zero_test_threads_are_rejected_before_io() {
    {
        let threads = "0";
        assert_unsupported(
            &[
                "test",
                "-k",
                "missing.count",
                "-t",
                "T",
                "-c",
                "C",
                "-n",
                "unused",
                "--threads",
                threads,
            ],
            "--threads must be greater than zero",
        );
    }
}

#[test]
fn mle_inactive_options_still_fail_without_publishing() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input");
    fs::write(&input, "not parsed").unwrap();
    let input = input.to_str().unwrap();
    let prefix = dir.path().join("out");
    let prefix = prefix.to_str().unwrap();
    for extra in [vec!["--threads", "0"], vec!["--seed", "0"]] {
        let mut args = vec!["mle", "-k", input, "-d", input, "-n", prefix];
        args.extend(extra.iter().copied());
        assert_unsupported(&args, extra[0]);
    }
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
}

#[test]
fn help_reports_reserved_commands_and_options() {
    let root = invoke(&["--help"]);
    assert!(root.status.success());
    let text = String::from_utf8(root.stdout).unwrap();
    assert!(text.contains("Unsupported: reserved"));
    let long = invoke(&["--help"]);
    assert!(!String::from_utf8(long.stdout)
        .unwrap()
        .contains("statically linked"));
    for sub in ["count", "test", "mle", "pathway", "run"] {
        let help = invoke(&[sub, "--help"]);
        assert!(help.status.success());
        let text = String::from_utf8(help.stdout).unwrap();
        if sub == "test" {
            assert!(text.contains("Worker limit N, auto, or auto:MAX"), "{text}");
        } else {
            assert!(text.contains("Unsupported:"), "{sub}: {text}");
        }
        if sub == "pathway" {
            assert!(text.contains("All"));
            assert!(text.contains("reserved"));
        }
    }
}

#[test]
fn cli_audit_covers_every_advertised_long_option() {
    let audit = format!(
        "{}\n{}",
        include_str!("../docs/public-surface/M10-CLI-AUDIT.md"),
        include_str!("../docs/public-surface/CAMPAIGN-A.md")
    );
    let root = Cli::command();
    for command in std::iter::once(&root).chain(root.get_subcommands()) {
        for arg in command.get_arguments() {
            if let Some(long) = arg.get_long() {
                assert!(
                    audit.contains(&format!("`--{long}`")),
                    "Missing {} --{} in audit",
                    command.get_name(),
                    long
                );
            }
        }
    }
}
