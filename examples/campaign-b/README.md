# Deterministic CPU examples

From the repository root, run the release binary with a fresh output prefix:

```sh
mageckrs mle -k tests/fixtures/m11a/counts.tsv \
  -d tests/fixtures/m11a/design.tsv --norm-method none \
  --threads 4 -n /tmp/mageckrs-campaign-b-mle

mageckrs test -k tests/fixtures/m11a/counts.tsv \
  -c Day0 -t TreatA,TreatB --norm-method none \
  --seed 0 --permutation 129 --threads 4 -n /tmp/mageckrs-campaign-b-test

mageckrs run -l examples/campaign-a/library.tsv \
  -f examples/campaign-a/sample-0.fastq,examples/campaign-a/sample-1.fastq \
     examples/campaign-a/sample-2.fastq examples/campaign-a/sample-3.fastq \
  --sample-label C,T1,T2 -c C -t T1,T2 --sgrna-len 4 --norm-method none \
  --permutation 129 --threads 4 -n /tmp/mageckrs-campaign-b-run
```

These are small retained fixtures for usage and parity; they do not demonstrate
real-screen acceleration. Remove `--threads` or use 1 for serial execution.
Combined counting remains serial. Output files must not already exist.
