# MAGeCKrs

MAGeCKrs 0.1.0 provides CRISPR screen analysis from Rust source on macOS Apple Silicon (`aarch64-apple-darwin`). The crate, library and command are `mageckrs`.

This is a source-only release under the [MIT License](LICENSE), Copyright (c) 2026 Yoon Sing Yap. No executable or application is distributed. [Citation metadata](CITATION.cff) and [third-party notices](THIRD_PARTY_NOTICES.md) accompany the source.

## Build from source

Install Rust/Cargo and Apple command-line developer tools on an Apple Silicon Mac. The release has been validated with Rust 1.98.0; the declared minimum Rust version has not been independently qualified.

```sh
git clone --branch v0.1.0 --depth 1 https://github.com/Yoonecorn/mageckrs.git
cd mageckrs
cargo build --release --locked
PATH="$PWD/target/aarch64-apple-darwin/release:$PATH"
mageckrs --version
mageckrs --help
```

If `CARGO_TARGET_DIR` is set, use its `aarch64-apple-darwin/release` directory instead. Builds link dynamically to macOS system and Accelerate libraries and use `target-cpu=native`. Intel macOS, Linux, Windows and universal binaries are unsupported. Source builds do not require an Apple Developer ID, signing or notarization.

## Workflows

Implemented commands are COUNT, exact explicitly scoped UMI counting, ordered unmapped output, COUNT-to-test `run`, RRA/`test`, classical MLE and a beta-only CNV-to-MLE adapter. COUNT remains serial even when a positive fixed `--threads` limit is supplied. RRA, run analysis and per-gene MLE support deterministic CPU workers; `auto` and `auto:MAX` use host-dependent scheduler metadata. Fixed workers are preferred for reproducible workflows.

Runnable synthetic examples from the repository root:

```sh
OUT=$(mktemp -d)
mageckrs count -l examples/campaign-c/library.tsv -f examples/campaign-c/day0.fastq examples/campaign-c/treat-a.fastq examples/campaign-c/treat-b.fastq -n "$OUT/counts" --sample-label Day0,TreatA,TreatB
mageckrs test -k "$OUT/counts.count.txt" -c Day0 -t TreatA,TreatB -n "$OUT/comparison"
mageckrs mle -k "$OUT/counts.count.txt" -d examples/campaign-c/design.tsv -n "$OUT/model"
mageckrs run -l examples/campaign-c/library.tsv -f examples/campaign-c/day0.fastq examples/campaign-c/treat-a.fastq examples/campaign-c/treat-b.fastq -n "$OUT/combined" --sample-label Day0,TreatA,TreatB -c Day0 -t TreatA,TreatB
mageckrs mle -k "$OUT/counts.count.txt" -d examples/campaign-c/design.tsv -n "$OUT/cnv-model" --cnv-norm examples/campaign-c/cnv.tsv
```

## Interpretation and limitations

RRA p-values and FDRs are nominal diagnostics of the implemented rank model. Interpret them with effect direction, gene ranking, guide consistency, experimental replication and biological validation. Cross-tool outputs need not be identical; universal biological calibration and guaranteed false-discovery control are not claimed.

MLE computes ordinary Wald inference. Permutation-style columns remain 1.0 placeholders. CNV correction changes matched betas without recomputing covariance, standard errors, z statistics, p-values or FDR. Same-model retries and failures are reported in diagnostic sidecars; numerical convergence does not establish biological truth.

Pathway/GSEA, MLE permutation inference, supplied MLE seeds and operational GPU/Metal acceleration are deferred. Metal compile checking does not qualify GPU execution. UMI/CNV interfaces do not establish biological accuracy. Atomic publication of a complete output set across crashes is not guaranteed.

See [usage and output semantics](docs/usage.md), [UMI and combined-run details](docs/public-surface/CAMPAIGN-A.md), [numerical acceptance](docs/numerical-contract.md), and the [synthetic examples](examples/campaign-a/README.md). Manuscript material and benchmark evidence are maintained separately and are not included in this code repository. No DOI or preprint deposition is asserted.

## Validate

```sh
cargo test --all-targets --locked
cargo test --test e2e_validation --locked
cargo test --doc --locked
cargo check --all-targets --locked
cargo check --all-targets --features metal --locked
cargo clippy --all-targets --locked -- -D warnings
cargo fmt --check
```

The public suite uses synthetic examples and permitted support fixtures. Real-data regression fixtures and their dependent tests are excluded because redistribution rights were not established. No benchmark is run by the four empty Cargo benchmark targets.
