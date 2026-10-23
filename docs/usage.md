# Usage and output semantics

COUNT sample groups are space-separated; comma-separated paths within a group are technical replicates. Plain and gzip FASTQ are supported. Physical splits are summed. Fixed-offset exact or unique-one-mismatch matching and exact-only automatic offset scanning have distinct semantics; inspect `count --help`.

Exact UMI counting requires `--umi-start`, `--umi-length`, `--umi-dedup-scope file|sample`, and `--umi-max-keys` together. Coordinates address the original read and may not overlap the matched guide. Counts are molecules; QC remains physical reads. The key cap limits distinct keys, not total memory. See [details](public-surface/CAMPAIGN-A.md).

RRA selects control/treatment labels or zero-based indexes and retains control/all-zero sentinel rows according to its output policy. The default seed is fixed; metadata record the effective seed. Multi-control negative-binomial dispersion is the maximum of raw and fitted values; the separate single-control adapter uses pooled moments, centered weighted least squares and fitted-only dispersion.

MLE aligns samples with a full-rank design matrix, uses default scalar dispersion 0.01 and computes ordinary Wald inference. Default maximum guides per gene is 500. Control/oversized genes can be excluded from inference. Numerical checks and retry rules are in [the numerical contract](numerical-contract.md).

Use fresh output prefixes. Normalized counts have exactly four decimal places; gene-level columns use pipe delimiters. Staged UMI/unmapped COUNT, RRA/test, MLE and run output publication uses sequential per-file renames with best-effort rollback. Ordinary COUNT retains overwrite behavior. Crash-wide atomicity of the whole file set is not guaranteed. `run` performs COUNT then RRA/test and has no MLE or CNV stage.

Use `mageckrs --help` and subcommand `--help` for supported options. `--count-n`, pathway execution and supplied MLE seeds are rejected. MLE permutation columns are placeholders, CNV does not recompute inference, and operational Metal acceleration is unqualified. CPU execution remains available.
