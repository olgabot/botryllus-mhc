#!/usr/bin/env python3
"""Run one kmerseek 0.4.0 arm of notebook 233: Botryllus proteins against human GENCODE v49
canonical proteins with one alphabet at its bit-matched k (a row of 233_arms.tsv).

Steps, each skipped when its output already exists, so a rerun picks up where it stopped:

1. Index the human proteome: scaled 1, low-complexity k-mers removed, the alphabet's own
   extension penalty and X-drop, and a Karlin-Altschul fit on 200 database sequences
   (the settings of notebook 232's hp_lehninger2 run, PR 10).
2. If the fit is refused, refit once with 1,000 calibration sequences
   (`kmerseek calibrate`). If that is refused too, the arm has no E-values: every match
   gets E = inf, and summary.json says so.
3. Search all Botryllus proteins with the filters of notebook 232's run: a match is kept
   when the whole-query p-value is below 0.05 or the best region's score is at least 1.3.
4. Copy the matches with E <= 10 or E = inf to a smaller gzipped CSV, the only file the
   notebook reads. E = inf means kmerseek found shared k-mers but did not extend them.

The index and the full search CSV go to WORK_DIR, node-local disk on Sherlock (RocksDB
opens and closes many files, which overloads $SCRATCH's metadata servers), and are
deleted when the arm is done. Everything kept goes to OUT_DIR.

Standard library only (Sherlock's python module has no polars).

Usage:
    233_kmerseek_alphabet_arm.py <row number in 233_arms.tsv, 0 = first arm>
Environment: KMERSEEK (binary), HUMAN (target FASTA), QUERY (Botryllus FASTA), OUT_DIR,
WORK_DIR (scratch space for the index and the full search CSV), and optionally ARMS_TSV
(the table of arms; default 233_arms.tsv next to this script).
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# notebook 235 runs the same arms from its own table (235_arms.tsv)
ARMS_TSV = Path(os.environ.get("ARMS_TSV") or Path(__file__).with_name("233_arms.tsv"))
KMERSEEK = os.environ["KMERSEEK"]
HUMAN = os.environ["HUMAN"]
QUERY = os.environ["QUERY"]
OUT_DIR = Path(os.environ["OUT_DIR"])
WORK_DIR = Path(os.environ["WORK_DIR"])

KA_QUERIES, KA_SEED, RETRY_KA_QUERIES = 200, 1, 1_000
MAX_QUERY_PVALUE, MIN_REGION_SCORE = 0.05, 1.3
MAX_EVALUE = 10.0
csv.field_size_limit(sys.maxsize)


def run(cmd: list[str], log: Path) -> tuple[int, float]:
    """Run a command, writing its output to a log; return (exit code, seconds)."""
    t0 = time.time()
    with open(log, "w") as fh:
        fh.write(" ".join(cmd) + "\n\n")
        fh.flush()
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
    return rc, round(time.time() - t0, 1)


def fitted(survival_csv: Path) -> bool | None:
    """The `fitted` flag kmerseek 0.4.0 writes in its survival CSV; None if missing."""
    if not survival_csv.exists() or survival_csv.stat().st_size == 0:
        return None
    with open(survival_csv) as fh:
        row = next(csv.DictReader(fh), None)
    return None if row is None else row.get("fitted", "").strip().lower() == "true"


def main() -> None:
    row_number = int(sys.argv[1])
    with open(ARMS_TSV) as fh:
        arm = list(csv.DictReader(fh, delimiter="\t"))[row_number]
    alphabet, k, penalty, xdrop = arm["alphabet"], arm["ksize"], arm["penalty"], arm["xdrop"]
    tag = f"{alphabet}.k{k}.scaled1.lcremoved"
    for sub in ("ka_survival", "evalue_le_10", "logs"):
        (OUT_DIR / sub).mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    kept = OUT_DIR / "evalue_le_10" / f"{tag}.evalue_le_10_or_inf.csv.gz"
    if kept.exists():
        print(f"{tag}: already done, {kept} exists")
        return
    idx = WORK_DIR / f"human_gencode_v49_canonical.{tag}.rocksdb"
    survival = OUT_DIR / "ka_survival" / f"{tag}.csv"
    retry_survival = OUT_DIR / "ka_survival" / f"{tag}.retry.csv"
    raw = WORK_DIR / f"{tag}.search.csv"
    logs = OUT_DIR / "logs"
    summary = {**arm, "tag": tag, "kmerseek": KMERSEEK}
    extension = ["--extend-mismatch-penalty", penalty, "--extend-xdrop", xdrop]

    if not (idx / "CURRENT").exists():
        summary["index_rc"], summary["index_s"] = run(
            [KMERSEEK, "index", "--input", HUMAN, "--output", str(idx), "--alphabet", alphabet,
             "--ksize", k, "--scaled", "1", "--remove-low-complexity", *extension,
             "--ka-queries", str(KA_QUERIES), "--ka-seed", str(KA_SEED),
             "--ka-survival-out", str(survival)],
            logs / f"{tag}.index.log")
        if summary["index_rc"] != 0:
            sys.exit(f"{tag}: index failed, see {logs / f'{tag}.index.log'}")

    summary["fitted_first"] = fitted(survival)
    fit = summary["fitted_first"]
    if fit is False:
        if fitted(retry_survival) is None:
            summary["calibrate_rc"], summary["calibrate_s"] = run(
                [KMERSEEK, "calibrate", "--target", str(idx), *extension,
                 "--ka-queries", str(RETRY_KA_QUERIES), "--ka-seed", str(KA_SEED),
                 "--ka-survival-out", str(retry_survival)],
                logs / f"{tag}.calibrate.log")
        fit = fitted(retry_survival)
    summary["fitted"] = fit

    if not raw.exists():
        summary["search_rc"], summary["search_s"] = run(
            [KMERSEEK, "search", "--query", QUERY, "--target", str(idx), "--alphabet", alphabet,
             "--ksize", k, *extension, "--ka-queries", str(KA_QUERIES), "--ka-seed", str(KA_SEED),
             "--max-query-pvalue", str(MAX_QUERY_PVALUE), "--min-region-score", str(MIN_REGION_SCORE),
             "--output", str(raw)],
            logs / f"{tag}.search.log")
        if summary["search_rc"] != 0:
            sys.exit(f"{tag}: search failed, see {logs / f'{tag}.search.log'}")

    n_all = n_kept = n_inf = 0
    tmp = kept.with_suffix(".tmp")
    with open(raw, newline="") as fin, gzip.open(tmp, "wt", newline="") as fout:
        reader = csv.reader(fin)
        header = next(reader)
        col = header.index("region_evalue")
        writer = csv.writer(fout)
        writer.writerow(header)
        for row in reader:
            n_all += 1
            try:
                e = float(row[col])
            except ValueError:
                continue
            if math.isinf(e):
                n_inf += 1
            if e <= MAX_EVALUE or math.isinf(e):
                n_kept += 1
                writer.writerow(row)
    tmp.rename(kept)
    raw_bytes = raw.stat().st_size
    raw.unlink()
    shutil.rmtree(idx)
    summary.update(n_matches=n_all, n_evalue_le_10_or_inf=n_kept, n_evalue_inf=n_inf,
                   raw_csv_bytes=raw_bytes)
    (OUT_DIR / "evalue_le_10" / f"{tag}.summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
