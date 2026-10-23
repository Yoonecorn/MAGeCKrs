use mageckrs::error::MageckError;

#[test]
fn cnv_error_mapping_uses_existing_public_math_category() {
    // The CNV kernel itself remains crate-private; this is its public error boundary.
    let error = MageckError::Math("CNV population is empty".to_string());
    assert_eq!(
        error.to_string(),
        "Math/Modeling error: CNV population is empty"
    );
    assert!(std::error::Error::source(&error).is_none());
}

#[test]
fn public_library_dispatch_connects_cnv_to_mle() {
    use clap::Parser;
    use mageckrs::cli::Cli;
    use std::{fs, path::Path};
    let root = tempfile::tempdir().unwrap();
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/m11a");
    let prefix = root.path().join("library");
    let cli = Cli::try_parse_from([
        "mageckrs",
        "mle",
        "-k",
        fixture.join("counts.tsv").to_str().unwrap(),
        "-d",
        fixture.join("design.tsv").to_str().unwrap(),
        "-n",
        prefix.to_str().unwrap(),
        "--cnv-norm",
        fixture.join("cnv.tsv").to_str().unwrap(),
        "--max-sgrnas-per-gene",
        "2",
    ])
    .unwrap();
    mageckrs::run_with_cli(cli).unwrap();
    let info = fs::read_to_string(root.path().join("library.cnv_run_info.txt")).unwrap();
    assert!(info.contains("M11A-CNV-ADAPTER-1"));
    assert!(info.contains("A\tA\tmatched-genes\t10\n"));
    assert_eq!(fs::read_dir(root.path()).unwrap().count(), 4);
}
