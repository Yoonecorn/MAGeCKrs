use clap::{CommandFactory, Parser};
use mageckrs::cli::{
    Cli, Commands, CountArgs, MleArgs, NormMethod, PathwayArgs, PathwayMethod, RemoveZeroMode,
    RunArgs, TestArgs, DEFAULT_REMOVE_ZERO_THRESHOLD,
};
use mageckrs::config::{
    CountConfig, FastqSampleGroup, MleConfig, PathwayConfig, RunConfig, TestConfig,
};
use mageckrs::error::MageckError;
use std::io::Write;
use std::path::PathBuf;
use tempfile::NamedTempFile;

#[test]
fn test_top_level_help() {
    let mut cmd = Cli::command();
    let mut help_buf = Vec::new();
    cmd.write_help(&mut help_buf).expect("Failed to write help");
    let help_text = String::from_utf8(help_buf).expect("Help is not UTF-8");

    assert!(
        help_text.contains("count"),
        "Missing count subcommand in help"
    );
    assert!(
        help_text.contains("test"),
        "Missing test subcommand in help"
    );
    assert!(help_text.contains("mle"), "Missing mle subcommand in help");
    assert!(
        help_text.contains("pathway"),
        "Missing pathway subcommand in help"
    );
    assert!(help_text.contains("run"), "Missing run subcommand in help");

    let subcommands: Vec<_> = cmd.get_subcommands().map(|s| s.get_name()).collect();
    assert_eq!(
        subcommands,
        vec!["count", "test", "mle", "pathway", "run"],
        "Subcommand list must match exactly"
    );
}

#[test]
fn test_subcommand_helps() {
    let cmd = Cli::command();
    for sub_name in ["count", "test", "mle", "pathway", "run"] {
        let sub = cmd
            .find_subcommand(sub_name)
            .unwrap_or_else(|| panic!("Subcommand {sub_name} not found"));
        let mut help_buf = Vec::new();
        let mut sub_clone = sub.clone();
        sub_clone.write_help(&mut help_buf).unwrap();
        let help_text = String::from_utf8(help_buf).unwrap();
        assert!(
            !help_text.is_empty(),
            "Help text for {sub_name} should not be empty"
        );
    }
}

#[test]
fn test_count_missing_required_args() {
    // Missing all args
    let result = Cli::try_parse_from(["mageckrs", "count"]);
    assert!(result.is_err());

    // Missing fastq
    let result = Cli::try_parse_from(["mageckrs", "count", "-l", "lib.tsv", "-n", "out"]);
    assert!(result.is_err());

    // Missing output-prefix
    let result = Cli::try_parse_from(["mageckrs", "count", "-l", "lib.tsv", "-f", "s1.fq"]);
    assert!(result.is_err());
}

#[test]
fn test_count_config_validation() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();

    let valid_args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![
            fq1.path().display().to_string(),
            fq2.path().display().to_string(),
        ],
        output_prefix: "sample_out".to_string(),
        sample_label: Some(vec!["ctrl".to_string(), "trt".to_string()]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    let config =
        CountConfig::try_from(valid_args).expect("Valid CountArgs must parse into CountConfig");
    assert_eq!(config.output_prefix, "sample_out");
    assert_eq!(config.sample_groups.len(), 2);
    assert_eq!(config.sample_groups[0].label, "ctrl");
    assert_eq!(
        config.sample_groups[0].fastq_paths,
        vec![fq1.path().to_path_buf()]
    );
    assert_eq!(config.sample_groups[1].label, "trt");
    assert_eq!(
        config.sample_groups[1].fastq_paths,
        vec![fq2.path().to_path_buf()]
    );
    assert_eq!(config.sgrna_len, 20);

    // Missing library file
    let invalid_lib_args = CountArgs {
        list_seq: PathBuf::from("/nonexistent/path/lib.tsv"),
        fastq: vec![fq1.path().display().to_string()],
        output_prefix: "out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };
    match CountConfig::try_from(invalid_lib_args) {
        Err(MageckError::Config(msg)) => assert!(msg.contains("does not exist")),
        other => panic!("Expected MageckError::Config, got {other:?}"),
    }

    // Mismatched sample labels count
    let mismatched_labels_args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![
            fq1.path().display().to_string(),
            fq2.path().display().to_string(),
        ],
        output_prefix: "out".to_string(),
        sample_label: Some(vec!["single_label".to_string()]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };
    match CountConfig::try_from(mismatched_labels_args) {
        Err(MageckError::Config(msg)) => assert!(msg.contains("does not match")),
        other => panic!("Expected MageckError::Config, got {other:?}"),
    }

    // Invalid sgrna length (> 32)
    let invalid_sgrna_len_args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![fq1.path().display().to_string()],
        output_prefix: "out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 40,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };
    match CountConfig::try_from(invalid_sgrna_len_args) {
        Err(MageckError::Config(msg)) => assert!(msg.contains("must be between 1 and 32 nt")),
        other => panic!("Expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_test_config_validation() {
    let mut count_file = NamedTempFile::new().unwrap();
    writeln!(count_file, "sgRNA\tGene\tctrl\ttrt").unwrap();

    let valid_args = TestArgs {
        count_table: count_file.path().to_path_buf(),
        treatment_id: vec!["trt".to_string()],
        day0_label: vec!["ctrl".to_string()],
        output_prefix: "test_out".to_string(),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        gene_test_fdr_threshold: 0.25,
        adjust_method: "fdr".to_string(),
        rra_seed: None,
        num_permutations: None,
        permutation: 1000,
        remove_zero: RemoveZeroMode::Both,
        remove_zero_threshold: DEFAULT_REMOVE_ZERO_THRESHOLD,
        seed: Some(42),
        threads: None,
    };

    let config =
        TestConfig::try_from(valid_args).expect("Valid TestArgs must parse into TestConfig");
    assert_eq!(config.output_prefix, "test_out");
    assert_eq!(config.permutation, 1000);
    assert_eq!(config.seed, Some(42));

    // Zero permutations error
    let zero_perm_args = TestArgs {
        count_table: count_file.path().to_path_buf(),
        treatment_id: vec!["trt".to_string()],
        day0_label: vec!["ctrl".to_string()],
        output_prefix: "test_out".to_string(),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        gene_test_fdr_threshold: 0.25,
        adjust_method: "fdr".to_string(),
        rra_seed: None,
        num_permutations: None,
        permutation: 0,
        remove_zero: RemoveZeroMode::Both,
        remove_zero_threshold: DEFAULT_REMOVE_ZERO_THRESHOLD,
        seed: None,
        threads: None,
    };
    match TestConfig::try_from(zero_perm_args) {
        Err(MageckError::Config(msg)) => assert!(msg.contains("permutations")),
        other => panic!("Expected MageckError::Config, got {other:?}"),
    }

    // Missing count table file
    let missing_table_args = TestArgs {
        count_table: PathBuf::from("/nonexistent/counts.txt"),
        treatment_id: vec!["trt".to_string()],
        day0_label: vec!["ctrl".to_string()],
        output_prefix: "test_out".to_string(),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        gene_test_fdr_threshold: 0.25,
        adjust_method: "fdr".to_string(),
        rra_seed: None,
        num_permutations: None,
        permutation: 1000,
        remove_zero: RemoveZeroMode::Both,
        remove_zero_threshold: DEFAULT_REMOVE_ZERO_THRESHOLD,
        seed: None,
        threads: None,
    };
    match TestConfig::try_from(missing_table_args) {
        Err(MageckError::Config(msg)) => assert!(msg.contains("does not exist")),
        other => panic!("Expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_remove_zero_default_is_both() {
    let cli = Cli::try_parse_from([
        "mageckrs",
        "test",
        "-k",
        "counts.txt",
        "-t",
        "trt",
        "-c",
        "ctrl",
        "-n",
        "out",
    ])
    .unwrap();
    match cli.command {
        Commands::Test(args) => {
            assert_eq!(args.remove_zero, RemoveZeroMode::Both);
            assert_eq!(args.remove_zero_threshold, 0.0);
        }
        _ => panic!("expected test arguments"),
    }
}

#[test]
fn test_remove_zero_accepts_exact_five_modes() {
    for (spelling, expected) in [
        ("none", RemoveZeroMode::None),
        ("control", RemoveZeroMode::Control),
        ("treatment", RemoveZeroMode::Treatment),
        ("any", RemoveZeroMode::Any),
        ("both", RemoveZeroMode::Both),
    ] {
        let cli = Cli::try_parse_from([
            "mageckrs",
            "test",
            "-k",
            "counts.txt",
            "-t",
            "trt",
            "-c",
            "ctrl",
            "-n",
            "out",
            "--remove-zero",
            spelling,
        ])
        .unwrap();
        match cli.command {
            Commands::Test(args) => assert_eq!(args.remove_zero, expected),
            _ => panic!("expected test arguments"),
        }
    }
}

#[test]
fn test_remove_zero_rejects_unknown_and_case_altered_modes() {
    for mode in ["unknown", "Both", "CONTROL"] {
        assert!(Cli::try_parse_from([
            "mageckrs",
            "test",
            "-k",
            "counts.txt",
            "-t",
            "trt",
            "-c",
            "ctrl",
            "-n",
            "out",
            "--remove-zero",
            mode,
        ])
        .is_err());
    }
}

#[test]
fn test_remove_zero_threshold_default_domain_and_help() {
    let count_file = NamedTempFile::new().unwrap();
    let valid = |threshold| TestArgs {
        count_table: count_file.path().to_path_buf(),
        treatment_id: vec!["trt".into()],
        day0_label: vec!["ctrl".into()],
        output_prefix: "out".into(),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        gene_test_fdr_threshold: 0.25,
        adjust_method: "fdr".into(),
        rra_seed: None,
        num_permutations: None,
        permutation: 1,
        remove_zero: RemoveZeroMode::Both,
        remove_zero_threshold: threshold,
        seed: None,
        threads: None,
    };
    assert!(TestConfig::try_from(valid(0.5)).is_ok());
    for threshold in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            TestConfig::try_from(valid(threshold)),
            Err(MageckError::Config(_))
        ));
    }

    let mut command = Cli::command();
    let test = command.find_subcommand_mut("test").unwrap();
    let mut help = Vec::new();
    test.write_help(&mut help).unwrap();
    let help = String::from_utf8(help).unwrap();
    assert!(help.contains("--remove-zero-threshold"));
    assert!(help.contains("any requires both groups to pass"));
    assert!(help.contains("both requires either"));
}

#[test]
fn mle_genes_varmodeling_default_and_opt_in_pass_through() {
    let count = NamedTempFile::new().unwrap();
    let design = NamedTempFile::new().unwrap();
    let args = [
        "mageckrs",
        "mle",
        "-k",
        count.path().to_str().unwrap(),
        "-d",
        design.path().to_str().unwrap(),
        "-n",
        "out",
    ];
    let config = |parsed: Cli| {
        let mageckrs::cli::Commands::Mle(args) = parsed.command else {
            panic!("expected MLE")
        };
        MleConfig::try_from(args).unwrap()
    };
    assert_eq!(
        config(Cli::try_parse_from(args).unwrap()).genes_varmodeling,
        0
    );
    for value in [0, 1, 2000, usize::MAX] {
        let value_text = value.to_string();
        assert_eq!(
            config(
                Cli::try_parse_from(args.into_iter().chain(["--genes-varmodeling", &value_text]))
                    .unwrap()
            )
            .genes_varmodeling,
            value
        );
    }
    assert!(Cli::try_parse_from(args.into_iter().chain(["--genes-varmodeling=-1"])).is_err());
    let mut command = Cli::command();
    let mut help = Vec::new();
    command
        .find_subcommand_mut("mle")
        .unwrap()
        .write_help(&mut help)
        .unwrap();
    let help = String::from_utf8(help).unwrap();
    assert!(help.contains("--genes-varmodeling"));
    assert!(help.contains("default: 0, scalar alpha=0.01"));
    assert!(help.contains("Set to 2000 to match the MAGeCK standalone MLE entry point behavior."));
}

#[test]
fn test_mle_config_validation() {
    let count_file = NamedTempFile::new().unwrap();
    let design_file = NamedTempFile::new().unwrap();

    let valid_args = MleArgs {
        count_table: count_file.path().to_path_buf(),
        design_matrix: design_file.path().to_path_buf(),
        output_prefix: "mle_out".to_string(),
        norm_method: NormMethod::Median,
        max_sgrnas_per_gene: 500,
        genes_varmodeling: 0,
        control_sgrna: None,
        cnv_norm: None,
        threads: None,
        seed: None,
    };

    let mut zero_args = valid_args.clone();
    zero_args.max_sgrnas_per_gene = 0;
    assert!(
        matches!(MleConfig::try_from(zero_args), Err(MageckError::Config(message)) if message == "--max-sgrnas-per-gene must be >= 1")
    );

    let config = MleConfig::try_from(valid_args).expect("Valid MleArgs must parse into MleConfig");
    assert_eq!(config.output_prefix, "mle_out");

    // Missing design matrix file
    let missing_design_args = MleArgs {
        count_table: count_file.path().to_path_buf(),
        design_matrix: PathBuf::from("/nonexistent/design.txt"),
        output_prefix: "mle_out".to_string(),
        norm_method: NormMethod::Median,
        max_sgrnas_per_gene: 500,
        genes_varmodeling: 0,
        control_sgrna: None,
        cnv_norm: None,
        threads: None,
        seed: None,
    };
    match MleConfig::try_from(missing_design_args) {
        Err(MageckError::Config(msg)) => assert!(msg.contains("does not exist")),
        other => panic!("Expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_pathway_config_validation() {
    let gene_file = NamedTempFile::new().unwrap();
    let gmt_file = NamedTempFile::new().unwrap();

    let valid_args = PathwayArgs {
        gene_summary: gene_file.path().to_path_buf(),
        gmt_file: gmt_file.path().to_path_buf(),
        output_prefix: "pathway_out".to_string(),
        method: PathwayMethod::Both,
        permutation: 1000,
        seed: None,
        threads: None,
    };

    let config = PathwayConfig::try_from(valid_args)
        .expect("Valid PathwayArgs must parse into PathwayConfig");
    assert_eq!(config.output_prefix, "pathway_out");
    assert_eq!(config.method, PathwayMethod::Both);

    // Missing GMT file
    let missing_gmt_args = PathwayArgs {
        gene_summary: gene_file.path().to_path_buf(),
        gmt_file: PathBuf::from("/nonexistent/pathways.gmt"),
        output_prefix: "pathway_out".to_string(),
        method: PathwayMethod::Both,
        permutation: 1000,
        seed: None,
        threads: None,
    };
    match PathwayConfig::try_from(missing_gmt_args) {
        Err(MageckError::Config(msg)) => assert!(msg.contains("does not exist")),
        other => panic!("Expected MageckError::Config, got {other:?}"),
    }
}

#[test]
fn test_run_config_validation() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();

    let valid_args = RunArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![fq1.path().display().to_string()],
        treatment_id: vec!["trt".to_string()],
        day0_label: vec!["ctrl".to_string()],
        output_prefix: "run_out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        permutation: 1000,
        seed: None,
        threads: None,
    };

    let config = RunConfig::try_from(valid_args).expect("Valid RunArgs must parse into RunConfig");
    assert_eq!(config.output_prefix, "run_out");
    assert_eq!(config.treatment_id, vec!["trt"]);
}

#[test]
fn test_cli_parsing_subcommands() {
    let cli = Cli::try_parse_from([
        "mageckrs",
        "count",
        "-l",
        "library.tsv",
        "-f",
        "s1.fq",
        "s2.fq",
        "-n",
        "out",
    ])
    .expect("Cli parse count");
    match cli.command {
        Commands::Count(args) => {
            assert_eq!(args.list_seq, PathBuf::from("library.tsv"));
            assert_eq!(args.fastq.len(), 2);
            assert_eq!(args.output_prefix, "out");
        }
        _ => panic!("Expected Commands::Count"),
    }

    let cli = Cli::try_parse_from([
        "mageckrs",
        "test",
        "-k",
        "counts.txt",
        "-t",
        "trt1,trt2",
        "-c",
        "ctrl1,ctrl2",
        "-n",
        "test_out",
    ])
    .expect("Cli parse test");
    match cli.command {
        Commands::Test(args) => {
            assert_eq!(args.treatment_id, vec!["trt1", "trt2"]);
            assert_eq!(args.day0_label, vec!["ctrl1", "ctrl2"]);
        }
        _ => panic!("Expected Commands::Test"),
    }

    let cli = Cli::try_parse_from([
        "mageckrs",
        "mle",
        "-k",
        "counts.txt",
        "-d",
        "design.txt",
        "-n",
        "mle_out",
    ])
    .expect("Cli parse mle");
    match cli.command {
        Commands::Mle(args) => {
            assert_eq!(args.count_table, PathBuf::from("counts.txt"));
            assert_eq!(args.design_matrix, PathBuf::from("design.txt"));
        }
        _ => panic!("Expected Commands::Mle"),
    }

    let cli = Cli::try_parse_from([
        "mageckrs",
        "pathway",
        "-g",
        "genes.txt",
        "--gmt-file",
        "hallmarks.gmt",
        "-n",
        "pathway_out",
    ])
    .expect("Cli parse pathway");
    match cli.command {
        Commands::Pathway(args) => {
            assert_eq!(args.gene_summary, PathBuf::from("genes.txt"));
            assert_eq!(args.gmt_file, PathBuf::from("hallmarks.gmt"));
        }
        _ => panic!("Expected Commands::Pathway"),
    }

    let cli = Cli::try_parse_from([
        "mageckrs",
        "run",
        "-l",
        "library.tsv",
        "-f",
        "s1.fq",
        "-t",
        "trt",
        "-c",
        "ctrl",
        "-n",
        "run_out",
    ])
    .expect("Cli parse run");
    match cli.command {
        Commands::Run(args) => {
            assert_eq!(args.list_seq, PathBuf::from("library.tsv"));
            assert_eq!(args.treatment_id, vec!["trt"]);
        }
        _ => panic!("Expected Commands::Run"),
    }
}

#[test]
fn test_req_01_single_sample_group_single_fastq() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![fq1.path().display().to_string()],
        output_prefix: "out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    let config = CountConfig::try_from(args).expect("Must succeed");
    assert_eq!(config.sample_groups.len(), 1);
    assert_eq!(
        config.sample_groups[0],
        FastqSampleGroup {
            fastq_paths: vec![fq1.path().to_path_buf()],
            label: "sample1".to_string(),
        }
    );
}

#[test]
fn test_req_02_multiple_sample_groups_single_fastq_each() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();
    let fq3 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![
            fq1.path().display().to_string(),
            fq2.path().display().to_string(),
            fq3.path().display().to_string(),
        ],
        output_prefix: "out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    let config = CountConfig::try_from(args).expect("Must succeed");
    assert_eq!(config.sample_groups.len(), 3);
    assert_eq!(
        config.sample_groups[0].fastq_paths,
        vec![fq1.path().to_path_buf()]
    );
    assert_eq!(config.sample_groups[0].label, "sample1");
    assert_eq!(
        config.sample_groups[1].fastq_paths,
        vec![fq2.path().to_path_buf()]
    );
    assert_eq!(config.sample_groups[1].label, "sample2");
    assert_eq!(
        config.sample_groups[2].fastq_paths,
        vec![fq3.path().to_path_buf()]
    );
    assert_eq!(config.sample_groups[2].label, "sample3");
}

#[test]
fn test_req_03_single_sample_group_multiple_technical_replicates() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();
    let fq3 = NamedTempFile::new().unwrap();

    let grouped_arg = format!(
        "{},{},{}",
        fq1.path().display(),
        fq2.path().display(),
        fq3.path().display()
    );

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![grouped_arg],
        output_prefix: "out".to_string(),
        sample_label: Some(vec!["ctrl_grp".to_string()]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    let config = CountConfig::try_from(args).expect("Must succeed");
    assert_eq!(config.sample_groups.len(), 1);
    assert_eq!(
        config.sample_groups[0].fastq_paths,
        vec![
            fq1.path().to_path_buf(),
            fq2.path().to_path_buf(),
            fq3.path().to_path_buf()
        ]
    );
    assert_eq!(config.sample_groups[0].label, "ctrl_grp");
}

#[test]
fn test_req_04_multiple_sample_groups_with_multiple_technical_replicates() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();
    let fq3 = NamedTempFile::new().unwrap();
    let fq4 = NamedTempFile::new().unwrap();

    let grp1 = format!("{},{}", fq1.path().display(), fq2.path().display());
    let grp2 = format!("{},{}", fq3.path().display(), fq4.path().display());

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![grp1, grp2],
        output_prefix: "out".to_string(),
        sample_label: Some(vec!["ctrl".to_string(), "treat".to_string()]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    let config = CountConfig::try_from(args).expect("Must succeed");
    assert_eq!(config.sample_groups.len(), 2);
    assert_eq!(
        config.sample_groups[0].fastq_paths,
        vec![fq1.path().to_path_buf(), fq2.path().to_path_buf()]
    );
    assert_eq!(config.sample_groups[0].label, "ctrl");
    assert_eq!(
        config.sample_groups[1].fastq_paths,
        vec![fq3.path().to_path_buf(), fq4.path().to_path_buf()]
    );
    assert_eq!(config.sample_groups[1].label, "treat");
}

#[test]
fn test_req_05_preservation_of_sample_group_order() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq_z = NamedTempFile::new().unwrap();
    let fq_a = NamedTempFile::new().unwrap();
    let fq_m = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![
            fq_z.path().display().to_string(),
            fq_a.path().display().to_string(),
            fq_m.path().display().to_string(),
        ],
        output_prefix: "out".to_string(),
        sample_label: Some(vec![
            "grp_z".to_string(),
            "grp_a".to_string(),
            "grp_m".to_string(),
        ]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    let config = CountConfig::try_from(args).expect("Must succeed");
    assert_eq!(config.sample_groups[0].label, "grp_z");
    assert_eq!(
        config.sample_groups[0].fastq_paths,
        vec![fq_z.path().to_path_buf()]
    );
    assert_eq!(config.sample_groups[1].label, "grp_a");
    assert_eq!(
        config.sample_groups[1].fastq_paths,
        vec![fq_a.path().to_path_buf()]
    );
    assert_eq!(config.sample_groups[2].label, "grp_m");
    assert_eq!(
        config.sample_groups[2].fastq_paths,
        vec![fq_m.path().to_path_buf()]
    );
}

#[test]
fn test_req_06_preservation_of_physical_file_order_inside_groups() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq_first = NamedTempFile::new().unwrap();
    let fq_second = NamedTempFile::new().unwrap();
    let fq_third = NamedTempFile::new().unwrap();

    // Pass in custom order: third, first, second
    let grp = format!(
        "{},{},{}",
        fq_third.path().display(),
        fq_first.path().display(),
        fq_second.path().display()
    );

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![grp],
        output_prefix: "out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    let config = CountConfig::try_from(args).expect("Must succeed");
    assert_eq!(
        config.sample_groups[0].fastq_paths,
        vec![
            fq_third.path().to_path_buf(),
            fq_first.path().to_path_buf(),
            fq_second.path().to_path_buf()
        ]
    );
}

#[test]
fn test_req_07_omitted_labels_default_generation() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();
    let fq3 = NamedTempFile::new().unwrap();
    let fq4 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![
            fq1.path().display().to_string(),
            fq2.path().display().to_string(),
            fq3.path().display().to_string(),
            fq4.path().display().to_string(),
        ],
        output_prefix: "out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    let config = CountConfig::try_from(args).expect("Must succeed");
    let labels: Vec<&str> = config
        .sample_groups
        .iter()
        .map(|g| g.label.as_str())
        .collect();
    assert_eq!(labels, vec!["sample1", "sample2", "sample3", "sample4"]);
}

#[test]
fn test_req_08_explicit_labels_map_to_sample_groups() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();
    let fq3 = NamedTempFile::new().unwrap();
    let fq4 = NamedTempFile::new().unwrap();
    let fq5 = NamedTempFile::new().unwrap();

    let grp1 = format!(
        "{},{},{}",
        fq1.path().display(),
        fq2.path().display(),
        fq3.path().display()
    );
    let grp2 = format!("{},{}", fq4.path().display(), fq5.path().display());

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![grp1, grp2],
        output_prefix: "out".to_string(),
        sample_label: Some(vec!["day0".to_string(), "day14".to_string()]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    let config = CountConfig::try_from(args).expect("Must succeed");
    assert_eq!(config.sample_groups.len(), 2);
    assert_eq!(config.sample_groups[0].label, "day0");
    assert_eq!(config.sample_groups[0].fastq_paths.len(), 3);
    assert_eq!(config.sample_groups[1].label, "day14");
    assert_eq!(config.sample_groups[1].fastq_paths.len(), 2);
}

#[test]
fn test_req_09_label_count_mismatch_with_technical_replicates() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();
    let fq3 = NamedTempFile::new().unwrap();
    let fq4 = NamedTempFile::new().unwrap();

    let grp1 = format!("{},{}", fq1.path().display(), fq2.path().display());
    let grp2 = format!("{},{}", fq3.path().display(), fq4.path().display());

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![grp1, grp2],
        output_prefix: "out".to_string(),
        sample_label: Some(vec![
            "rep1".to_string(),
            "rep2".to_string(),
            "rep3".to_string(),
            "rep4".to_string(),
        ]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    match CountConfig::try_from(args) {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("does not match number of sample groups"),
                "Expected group mismatch error, got: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got: {other:?}"),
    }
}

#[test]
fn test_req_10_empty_first_label() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![
            fq1.path().display().to_string(),
            fq2.path().display().to_string(),
        ],
        output_prefix: "out".to_string(),
        sample_label: Some(vec!["".to_string(), "treatment".to_string()]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    match CountConfig::try_from(args) {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("must not be empty"),
                "Expected empty label error, got: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got: {other:?}"),
    }

    // CLI parse with --sample-label ,treatment
    let cli = Cli::try_parse_from([
        "mageckrs",
        "count",
        "-l",
        lib_file.path().to_str().unwrap(),
        "-f",
        fq1.path().to_str().unwrap(),
        fq2.path().to_str().unwrap(),
        "-n",
        "out",
        "--sample-label",
        ",treatment",
    ])
    .expect("Cli parse should succeed");
    if let Commands::Count(args) = cli.command {
        match CountConfig::try_from(args) {
            Err(MageckError::Config(msg)) => assert!(msg.contains("must not be empty")),
            other => panic!("Expected MageckError::Config, got: {other:?}"),
        }
    }
}

#[test]
fn test_req_11_empty_middle_label() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();
    let fq3 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![
            fq1.path().display().to_string(),
            fq2.path().display().to_string(),
            fq3.path().display().to_string(),
        ],
        output_prefix: "out".to_string(),
        sample_label: Some(vec![
            "control".to_string(),
            "".to_string(),
            "treatment".to_string(),
        ]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    match CountConfig::try_from(args) {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("must not be empty"),
                "Expected empty label error, got: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got: {other:?}"),
    }

    // CLI parse with --sample-label control,,treatment
    let cli = Cli::try_parse_from([
        "mageckrs",
        "count",
        "-l",
        lib_file.path().to_str().unwrap(),
        "-f",
        fq1.path().to_str().unwrap(),
        fq2.path().to_str().unwrap(),
        fq3.path().to_str().unwrap(),
        "-n",
        "out",
        "--sample-label",
        "control,,treatment",
    ])
    .expect("Cli parse should succeed");
    if let Commands::Count(args) = cli.command {
        match CountConfig::try_from(args) {
            Err(MageckError::Config(msg)) => assert!(msg.contains("must not be empty")),
            other => panic!("Expected MageckError::Config, got: {other:?}"),
        }
    }
}

#[test]
fn test_req_12_empty_final_label() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![
            fq1.path().display().to_string(),
            fq2.path().display().to_string(),
        ],
        output_prefix: "out".to_string(),
        sample_label: Some(vec!["control".to_string(), "".to_string()]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    match CountConfig::try_from(args) {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("must not be empty"),
                "Expected empty label error, got: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got: {other:?}"),
    }

    // CLI parse with --sample-label control,
    let cli = Cli::try_parse_from([
        "mageckrs",
        "count",
        "-l",
        lib_file.path().to_str().unwrap(),
        "-f",
        fq1.path().to_str().unwrap(),
        fq2.path().to_str().unwrap(),
        "-n",
        "out",
        "--sample-label",
        "control,",
    ])
    .expect("Cli parse should succeed");
    if let Commands::Count(args) = cli.command {
        match CountConfig::try_from(args) {
            Err(MageckError::Config(msg)) => assert!(msg.contains("must not be empty")),
            other => panic!("Expected MageckError::Config, got: {other:?}"),
        }
    }
}

#[test]
fn test_req_13_duplicate_labels_rejected() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![
            fq1.path().display().to_string(),
            fq2.path().display().to_string(),
        ],
        output_prefix: "out".to_string(),
        sample_label: Some(vec!["sample_a".to_string(), "sample_a".to_string()]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    match CountConfig::try_from(args) {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("Duplicate sample label 'sample_a'"),
                "Expected duplicate label error, got: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got: {other:?}"),
    }
}

#[test]
fn test_req_14_labels_containing_tab_rejected() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![fq1.path().display().to_string()],
        output_prefix: "out".to_string(),
        sample_label: Some(vec!["ctrl\tsample".to_string()]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    match CountConfig::try_from(args) {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("invalid control characters"),
                "Expected control char error, got: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got: {other:?}"),
    }
}

#[test]
fn test_req_15_labels_containing_carriage_return_rejected() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![fq1.path().display().to_string()],
        output_prefix: "out".to_string(),
        sample_label: Some(vec!["ctrl\rsample".to_string()]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    match CountConfig::try_from(args) {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("invalid control characters"),
                "Expected control char error, got: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got: {other:?}"),
    }
}

#[test]
fn test_req_16_labels_containing_newline_rejected() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![fq1.path().display().to_string()],
        output_prefix: "out".to_string(),
        sample_label: Some(vec!["ctrl\nsample".to_string()]),
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    match CountConfig::try_from(args) {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("invalid control characters"),
                "Expected control char error, got: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got: {other:?}"),
    }
}

#[test]
fn test_req_17_leading_empty_fastq_component_rejected() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![format!(",{}", fq1.path().display())],
        output_prefix: "out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    match CountConfig::try_from(args) {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("Empty FASTQ path in sample group 1"),
                "Expected empty FASTQ path error, got: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got: {other:?}"),
    }
}

#[test]
fn test_req_18_middle_empty_fastq_component_rejected() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();
    let fq2 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![format!(
            "{},,{}",
            fq1.path().display(),
            fq2.path().display()
        )],
        output_prefix: "out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    match CountConfig::try_from(args) {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("Empty FASTQ path in sample group 1"),
                "Expected empty FASTQ path error, got: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got: {other:?}"),
    }
}

#[test]
fn test_req_19_trailing_empty_fastq_component_rejected() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![format!("{},", fq1.path().display())],
        output_prefix: "out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    match CountConfig::try_from(args) {
        Err(MageckError::Config(msg)) => {
            assert!(
                msg.contains("Empty FASTQ path in sample group 1"),
                "Expected empty FASTQ path error, got: {msg}"
            );
        }
        other => panic!("Expected MageckError::Config, got: {other:?}"),
    }
}

#[test]
fn test_req_20_repeated_physical_fastq_paths_preserved() {
    let lib_file = NamedTempFile::new().unwrap();
    let fq1 = NamedTempFile::new().unwrap();

    let repeated = format!("{},{}", fq1.path().display(), fq1.path().display());
    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![repeated],
        output_prefix: "out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    let config = CountConfig::try_from(args).expect("Must succeed");
    assert_eq!(config.sample_groups.len(), 1);
    assert_eq!(
        config.sample_groups[0].fastq_paths,
        vec![fq1.path().to_path_buf(), fq1.path().to_path_buf()]
    );
    assert_eq!(config.sample_groups[0].fastq_paths.len(), 2);
}

#[test]
fn test_req_21_paths_containing_punctuation_preserved() {
    let lib_file = NamedTempFile::new().unwrap();
    let temp_dir = tempfile::tempdir().unwrap();
    let file_path = temp_dir.path().join("my-sample_1.2+3.fastq");
    std::fs::write(&file_path, b"").unwrap();

    let args = CountArgs {
        list_seq: lib_file.path().to_path_buf(),
        fastq: vec![file_path.display().to_string()],
        output_prefix: "out".to_string(),
        sample_label: None,
        norm_method: NormMethod::Median,
        control_sgrna: None,
        sgrna_len: 20,
        trim_5: mageckrs::cli::Trim5Mode::None,
        test_mismatch: 0,
        umi_start: None,
        umi_length: None,
        umi_dedup_scope: None,
        umi_max_keys: None,
        count_n: false,
        unmapped_to_file: false,
        threads: None,
    };

    let config = CountConfig::try_from(args).expect("Must succeed");
    assert_eq!(config.sample_groups.len(), 1);
    assert_eq!(config.sample_groups[0].fastq_paths[0], file_path);
}

#[test]
fn test_req_22_existing_valid_one_file_per_sample_cli() {
    let cli = Cli::try_parse_from([
        "mageckrs",
        "count",
        "-l",
        "library.tsv",
        "-f",
        "s1.fq",
        "s2.fq",
        "-n",
        "out",
        "--sample-label",
        "ctrl,trt",
    ])
    .expect("Cli parse must succeed");

    match cli.command {
        Commands::Count(args) => {
            assert_eq!(args.fastq, vec!["s1.fq", "s2.fq"]);
            assert_eq!(
                args.sample_label,
                Some(vec!["ctrl".to_string(), "trt".to_string()])
            );
        }
        _ => panic!("Expected Commands::Count"),
    }
}
