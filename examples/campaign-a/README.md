# Tiny reproducible Campaign A example

Run from the repository root after `cargo build --release --locked`. Supply a
fresh output directory; the commands below use `/tmp/mageckrs-campaign-a-example`
only as an example. Existing destinations are protected by run/opt-in COUNT.
The six guides and four physical FASTQ files are synthetic, with a two-base AA
suffix UMI. The first two files form one sample; the next two are distinct samples.

```sh
mkdir /tmp/mageckrs-campaign-a-example
BIN=target/aarch64-apple-darwin/release/mageckrs
EX=examples/campaign-a
$BIN count -l "$EX/library.tsv" -f "$EX/sample-0.fastq,$EX/sample-1.fastq" "$EX/sample-2.fastq" "$EX/sample-3.fastq" -n /tmp/mageckrs-campaign-a-example/counts --sample-label C1,C2,T --sgrna-len 4 --norm-method none
$BIN test -k /tmp/mageckrs-campaign-a-example/counts.count.txt -c C1,C2 -t T -n /tmp/mageckrs-campaign-a-example/test --norm-method none --permutation 19 --seed 17
$BIN run -l "$EX/library.tsv" -f "$EX/sample-0.fastq,$EX/sample-1.fastq" "$EX/sample-2.fastq" "$EX/sample-3.fastq" -n /tmp/mageckrs-campaign-a-example/combined --sample-label C1,C2,T -c C1,C2 -t T --sgrna-len 4 --norm-method none --permutation 19 --seed 17
$BIN count -l "$EX/library.tsv" -f "$EX/sample-0.fastq,$EX/sample-1.fastq" "$EX/sample-2.fastq" "$EX/sample-3.fastq" -n /tmp/mageckrs-campaign-a-example/molecules --sample-label C1,C2,T --sgrna-len 4 --norm-method none --umi-start 4 --umi-length 2 --umi-dedup-scope sample --umi-max-keys 6 --unmapped-to-file
```

Use the binary path reported by Cargo when `CARGO_TARGET_DIR` is external.
Without deduplication, g1 counts are C1=10, C2=11, T=2. With sample-scope UMI,
each guide has exactly one molecule per sample. With file scope, C1 has two
molecules per guide because its files deduplicate independently. Read QC retains
physical reads in both cases. This toy input has no unmapped reads, so its
unmapped file is empty. Change scope explicitly; there is no default.

Combined count/numerical test artifacts match the separate commands above.
Count metadata moves from `.run_info.txt` to `.count_run_info.txt`. Test metadata
records its own final prefix and is compared with that expected difference.
The campaign's process tests execute this example with temporary destinations.
No statistical benchmark or real-data validity claim follows from this fixture.
