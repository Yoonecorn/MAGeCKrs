# Exact UMI counting and combined run

## Exact UMI counting

Supply all four options together, with no implicit scope:

```sh
mageckrs count -l library.tsv -f lane1.fq,lane2.fq other.fq -n molecules \
  --sample-label baseline,treatment --sgrna-len 20 --trim-5 8 \
  --umi-start 0 --umi-length 8 --umi-dedup-scope sample --umi-max-keys 1000000
```

| Option | Meaning |
|---|---|
| `--umi-start` | Absolute zero-based position in the original sequence, before guide trimming |
| `--umi-length` | Fixed 1..=64 bases; checked half-open range start..start+length |
| `--umi-dedup-scope` | Required `file` or `sample`; no default |
| `--umi-max-keys` | Required positive maximum distinct `(guide, UMI)` keys in the active scope |

`file` deduplicates each physical FASTQ independently, then sums molecule counts
within its sample group. `sample` shares the deduper across every physical file
in that group and resets between sample groups. Even a repeated input path is a
separate physical-file occurrence for `file`. Equal UMIs on different guides or
samples remain independent. Exact uppercase A/C/G/T equality is the only
molecule rule: lowercase bases normalize to uppercase; no distance clustering,
quality correction, quality filtering, gene-wide collapse or error correction.
Packed keys apply through 32 bases; longer fixed UMIs use byte keys.

Guide matching follows ordinary COUNT: space-separated sample groups,
comma-separated technical replicates, existing exact/unique-one-mismatch
precedence. `--trim-5 auto` selects the first exact match at offsets 0..=40 and
remains exact-only even with `--test-mismatch 1`. A mapped UMI must be disjoint
from that selected guide window; overlap fails, rather than searching for another
guide offset. A mapped read with a short or non-ACGT UMI fails with file/record
diagnostics. Unmatched/ambiguous/short-guide records are not UMI observations.
Coordinate addition overflow and incomplete/zero/out-of-range UMI settings fail.

At capacity a duplicate remains legal. A new key beyond the active scope cap
fails before insertion, with no eviction or approximate fallback. Memory also
includes guide sets, reader buffers, matrices and allocator overhead: the key cap
is **not a complete process-memory bound**. A fatal parse/UMI/normalization error
publishes no successful final files in this opt-in staged workflow.

| Artifact | Meaning |
|---|---|
| `.count.txt` | Distinct molecules, one column per sample group, guide order unchanged |
| `.count.normalized.txt` | Existing normalization from molecule counts, exactly four decimals |
| `.countsummary.txt` | Existing six-column physical-read QC: File, Label, Reads, Mapped, Percentage, Zerocounts |
| `.umi_summary.txt` | File, Label, ValidMappedReads, NewMolecules, DuplicateMolecules in supplied order |
| `.run_info.txt` | Guide trim plus UMI coordinates, explicit scope, cap, equality and QC/count semantics |

Mapped and Zerocounts remain based on physical reads before deduplication, so
Mapped can exceed molecule counts. New molecules are attributed to their first
occurrence in supplied file order. Reordering technical replicates leaves sample
molecule totals unchanged but can change per-file NewMolecules/DuplicateMolecules
attribution. Empty files retain zero reads and all guides zero in read QC.
Header/index/paired-read UMI sources and lengths above 64 are unsupported.

## COUNT option dispositions and unmapped output

`--count-n` is explicitly retired from this bounded campaign. The compatibility
spelling is retained and always returns a configuration error before input I/O;
no wildcard assignment, N-aware library matching or ambiguity policy is implied.
COUNT accepts positive fixed worker limits but remains serial; zero is rejected.

`--unmapped-to-file` writes `<prefix>.unmapped.fastq`, an uncompressed stream in
sample-group, physical-file and record order. It contains ordinary no-match,
ambiguous and short-guide records, preserving ID, sequence and quality values.
The plus line is canonical `+`; lines use LF. No source tags are added and IDs
may repeat across inputs. Mapped UMI duplicates are mapped reads and are not
unmapped output. Invalid FASTQ/UMI aborts instead of being written as unmapped.
FASTQ is required for UMI/unmapped workflows; FASTA lacks quality and is rejected.
Plain and gzip FASTQ inputs may be mixed. Output compression is unsupported.
Auto retains the existing reader's short/no-match classification.

COUNT with UMI or unmapped output stages all requested artifacts and refuses
pre-existing destinations (including directories and dangling links). It checks
again atomically at each publication using macOS RENAME_EXCL. Ordinary COUNT
with neither option retains its earlier reader and overwrite/publication
behavior.

## Combined COUNT-to-test run

```sh
mageckrs run -l library.tsv -f lane1.fq,lane2.fq treatment.fq -n combined \
  --sample-label baseline,treatment -c baseline -t treatment \
  --norm-method median --permutation 1000 --seed 17
```

Run uses COUNT's grouped inputs and auto/none/integer trimming, including the
four UMI options and optional unmapped output. Roles use existing test label or
zero-based column-index resolution. Seed input is decimal u64; omission uses
existing fixed RRA seed 81985529216486895. Positive `--threads N` limits RRA workers while counting stays serial; omission
or 1 is serial. Retired count-n is rejected. This convenience command runs RRA/test only, with no MLE stage.
It exposes the existing permutation/seed/normalization/controls selection; other
test options retain their current defaults (BH/fdr, directional threshold 0.25,
remove-zero both and the existing strict threshold).

The shared COUNT preparation creates a validated role-neutral in-memory table;
no disk count TSV is reparsed. Raw values and nominal order are preserved.
COUNT writes its normalized artifact. Test separately derives its existing
normalization from the same **raw** matrix, exactly as `count` then `test -k
<prefix>.count.txt`; COUNT's normalized matrix is not fed into test. Both stages
retain their existing control-file grammar. COUNT resolves guide sequence or ID;
test resolves opaque guide IDs. Use IDs accepted by both for control normalization;
sequence-only controls can pass COUNT then fail the test stage. Existing error
precedence, RRA mathematics, ranks, draw order and CNV/MLE remain unchanged.

Run publishes the standard count and gene/sgRNA files under one prefix, plus any
UMI/unmapped outputs. Count metadata is named `.count_run_info.txt`, avoiding
collision with test's `.run_info.txt`. Count metadata bytes match standalone
COUNT for identical options. Test metadata retains command=test for the test
stage and records the final combined output prefix. When comparison prefixes
differ, the expected output_prefix difference is pinned separately; metadata
comparison is not mislabeled complete-file identity.

All requested artifacts stay staged until both stages succeed. Every final
path must be absent. Publication uses sequential exclusive renames; a late
conflict cannot overwrite an existing target. On failure, best-effort rollback
removes already-published owned outputs and reports cleanup failures. Temporary
stage cleanup failures are logged. This is not atomic publication of the entire
file set across crashes, and no crash-rollback guarantee is claimed.

The process tests compare full numerical output files for all four normalization
modes, both UMI scopes, grouped replicates, label/index roles and fixed/default
seeds. Unit tests inject late publication conflicts and rollback failures.
See the [synthetic example](../../examples/campaign-a/README.md).
No new Python oracle, measured UMI accuracy, real-CNV certification, GPU execution,
statistical significance calibration, speedup or release readiness is claimed.

RRA p-values/FDRs are nominal rank-model diagnostics. Interpret them with rank,
direction, guide consistency, reproducibility and biological validation; no
universal biological calibration is claimed. Experimental estimators and null models are outside the supported interface.
