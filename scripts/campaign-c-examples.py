#!/usr/bin/env python3
"""Execute current supported README examples with an explicitly selected binary.
Output must be a fresh directory. Standard library only; no biological claim.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import time


def identity(path):
    p = Path(path).resolve()
    return {"path": str(p), "bytes": p.stat().st_size,
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def rows(path, key):
    with Path(path).open() as stream:
        result = list(csv.DictReader(stream, delimiter="\t"))
    assert len({r[key] for r in result}) == len(result)
    return {r[key]: r for r in result}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    args = parser.parse_args()
    binary = Path(args.binary).resolve()
    root = Path(args.root).resolve()
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    report = {"binary": identity(binary), "root": str(root), "kind": "synthetic interface qualification",
              "runs": [], "assertions": [], "documents": []}

    def invoke(name, argv, expected=0):
        index = len(report["runs"])
        start = time.monotonic()
        proc = subprocess.run([str(binary), *argv], cwd=root, capture_output=True)
        stdout = out / f"{index:03}-{name}.stdout"
        stderr = out / f"{index:03}-{name}.stderr"
        stdout.write_bytes(proc.stdout)
        stderr.write_bytes(proc.stderr)
        rec = {"name": name, "argv": [str(binary), *argv], "exit": proc.returncode,
               "expected_exit": expected, "seconds": time.monotonic() - start,
               "stdout": identity(stdout), "stderr": identity(stderr)}
        report["runs"].append(rec)
        (out / "records.json").write_text(json.dumps(report, indent=2) + "\n")
        assert proc.returncode == expected, rec
        assert b"panicked" not in proc.stderr, rec
        return proc

    # Execute command text from current documents, not a separately maintained approximation.
    for doc, group in [("README.md", "root"), ("examples/campaign-a/README.md", "a"),
                       ("examples/campaign-b/README.md", "b")]:
        text = (root / doc).read_text()
        report["documents"].append(identity(root / doc))
        folder = out / group
        folder.mkdir()
        for block in re.findall(r"```sh\n(.*?)```", text, re.S):
            for line in block.replace("\\\n", " ").splitlines():
                line = line.strip()
                if not line.startswith(("mageckrs ", "$BIN ")):
                    continue
                tokens = shlex.split(line)[1:]
                argv = []
                for token in tokens:
                    token = token.replace("$OUT", str(folder))
                    token = token.replace("$EX", "examples/campaign-a")
                    token = token.replace("/tmp/mageckrs-campaign-a-example", str(folder))
                    token = token.replace("/tmp/mageckrs-campaign-b-", str(folder / "campaign-b-"))
                    argv.append(token)
                invoke(group + "-" + str(len(report["runs"])), argv)
    # Operational public-surface placeholders get concrete synthetic inputs.
    # Preserve argument text; only substitute named file/prefix placeholders.
    public = out / "public"
    public.mkdir()
    shutil_source = root / "examples/campaign-c"
    (public / "library.tsv").write_bytes((shutil_source / "library.tsv").read_bytes())
    for original, names in [("day0.fastq", ["lane1.fq", "lane2.fq"]), ("treat-a.fastq", ["other.fq", "treatment.fq"])]:
        lines = (shutil_source / original).read_text().splitlines()
        for name in names:
            with (public / name).open("w") as stream:
                for i in range(0, len(lines), 4):
                    stream.write(lines[i] + "\nAAAAAAAA" + lines[i + 1][:20] + "\n+\n" + "I" * 28 + "\n")
    for doc in ["docs/public-surface/CAMPAIGN-A.md"]:
        text = (root / doc).read_text()
        report["documents"].append(identity(root / doc))
        for block in re.findall(r"```sh\n(.*?)```", text, re.S):
            for line in block.replace("\\\n", " ").splitlines():
                if not line.strip().startswith("mageckrs "):
                    continue
                tokens = shlex.split(line)[1:]
                mapped = []
                input_folder = public
                if tokens[0] == "run":
                    input_folder = public / "run-inputs"
                    input_folder.mkdir()
                    (input_folder / "library.tsv").write_bytes((shutil_source / "library.tsv").read_bytes())
                    for original, names in [("day0.fastq", ["lane1.fq", "lane2.fq"]), ("treat-a.fastq", ["treatment.fq"])]:
                        for name in names:
                            (input_folder / name).write_bytes((shutil_source / original).read_bytes())
                for token in tokens:
                    if token in ["library.tsv", "lane1.fq,lane2.fq", "other.fq", "treatment.fq"]:
                        token = ",".join(str(input_folder / x) for x in token.split(","))
                    elif token in ["molecules", "combined"]:
                        token = str(public / token)
                    mapped.append(token)
                invoke("public-surface-a", mapped)
    invoke("public-surface-cnv", ["mle", "-k", str(root / "tests/fixtures/m11a/counts.tsv"),
                                "-d", str(root / "tests/fixtures/m11a/design.tsv"), "-n", str(public / "cnv"),
                                "--cnv-norm", str(root / "tests/fixtures/m11a/cnv.tsv")])
    for cmd in [[], ["count"], ["test"], ["mle"], ["run"], ["pathway"]]:
        invoke("help-" + (cmd[0] if cmd else "root"), [*cmd, "--help"])
    invoke("version", ["--version"])
    ex = root / "examples/campaign-c"
    count_args = ["count", "-l", str(ex / "library.tsv"), "-f", str(ex / "day0.fastq"),
                  str(ex / "treat-a.fastq"), str(ex / "treat-b.fastq"),
                  "--sample-label", "Day0,TreatA,TreatB", "--norm-method", "none"]
    for scope in ["file", "sample"]:
        invoke("umi-" + scope, [*count_args, "-n", str(out / ("umi-" + scope)),
                              "--umi-start", "20", "--umi-length", "2",
                              "--umi-dedup-scope", scope, "--umi-max-keys", "25"])
        molecule = rows(out / ("umi-" + scope + ".count.txt"), "sgRNA")
        source = rows(root / "tests/fixtures/m11a/counts.tsv", "sgRNA")
        for g, row in molecule.items():
            for sample in ["Day0", "TreatA", "TreatB"]:
                assert int(row[sample]) == int(int(source[g][sample]) > 0)
    report["assertions"].append("exact known molecule counts, both explicit scopes")
    # Unsupported options fail before analysis. Use valid inputs to distinguish these errors.
    mle_args = ["mle", "-k", str(root / "tests/fixtures/m11a/counts.tsv"),
                "-d", str(root / "tests/fixtures/m11a/design.tsv"), "-n", str(out / "rejected")]
    invoke("count-threads", [*count_args, "-n", str(out / "rejected"), "--threads", "0"], 1)
    invoke("count-serial-workers", [*count_args, "-n", str(out / "serial-workers"), "--threads", "1"], 0)
    invoke("count-n", [*count_args, "-n", str(out / "rejected"), "--count-n"], 1)
    invoke("mle-seed", [*mle_args, "--seed", "1"], 1)
    invoke("mle-zero-threads", [*mle_args, "--threads", "0"], 1)
    invoke("reserved-pathway", ["pathway", "-g", "absent", "--gmt-file", "absent", "-n", str(out / "rejected")], 1)
    invoke("run-cnv-unsupported", ["run", "--cnv-norm", str(ex / "cnv.tsv")], 2)
    assert not list(out.glob("rejected.*"))
    report["assertions"].append("unsupported options: exact exits, no outputs, no panics")
    # Literal counts, normalized precision, combined vs separate numerical output.
    source = rows(root / "tests/fixtures/m11a/counts.tsv", "sgRNA")
    emitted = rows(out / "root/counts.count.txt", "sgRNA")
    for g, row in emitted.items():
        for sample in ["Day0", "TreatA", "TreatB"]:
            assert row[sample] == source[g][sample]
    normalized = rows(out / "root/counts.count.normalized.txt", "sgRNA")
    assert all(re.fullmatch(r"[0-9]+\.[0-9]{4}", r[s]) for r in normalized.values()
               for s in ["Day0", "TreatA", "TreatB"])
    for group, separate in [("root", "comparison"), ("a", "test")]:
        for suffix in ["gene_summary.txt", "sgrna_summary.txt"]:
            assert (out / group / (separate + "." + suffix)).read_bytes() == (out / group / ("combined." + suffix)).read_bytes()
        assert (out / group / "counts.count.txt").read_bytes() == (out / group / "combined.count.txt").read_bytes()
    report["assertions"].append("known read counts; four-decimal output; separate/combined scientific bytes")
    normal = rows(out / "root/model.gene_summary.txt", "Gene")
    corrected = rows(out / "root/cnv-model.gene_summary.txt", "Gene")
    assert normal.keys() == corrected.keys()
    changed = 0
    for gene, row in normal.items():
        for col, value in row.items():
            if col.endswith("|beta"):
                changed += value != corrected[gene][col]
            else:
                assert value == corrected[gene][col], (gene, col)
            if col.endswith("|p-value") or col.endswith("|fdr"):
                assert value == "1.0", (gene, col, value)
    assert changed > 0
    assert (out / "root/model.sgrna_summary.txt").read_bytes() == (out / "root/cnv-model.sgrna_summary.txt").read_bytes()
    report["assertions"].append("CNV changes betas while preserving ordinary inference/eff; permutation columns remain 1")
    report["outputs"] = [identity(p) for p in sorted(out.rglob("*")) if p.is_file() and p.name != "records.json"]
    report["status"] = "PASS"
    (out / "records.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "executions": len(report["runs"]), "assertions": report["assertions"]}, indent=2))


if __name__ == "__main__":
    main()
