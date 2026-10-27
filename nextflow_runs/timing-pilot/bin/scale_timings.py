"""Scale the 500-protein pilot's measured times to the 45,339-protein Botryllus proteome.

Every step's time grows with the number of residues it reads, plus a fixed cost per task
(loading a model or an HMM library). The fixed cost comes from the 10-protein Sherlock run
of 2026-10-25 (fixed_costs_mini10.tsv), where almost all the time is that cost. So:

    whole proteome = fixed + (pilot - fixed) * (proteome residues / pilot residues)

proteinannotator drops proteins longer than 5,000 aa (its max_seq_length), so its residue
count leaves them out. ESMFold is scaled protein by protein instead: each of the 500 stands
for the 1/500 of the proteome in its length decile, so its fold time counts 45,339 / 500
times.

Usage:
    python scale_timings.py <raw dir> <pilot lengths.tsv> <proteome.fasta> <out dir>
"""

import re
import sys
from datetime import datetime
from pathlib import Path

import polars as pl

PA_MAX_LEN = 5000


def residues(fasta, max_len=None):
    n, lengths, cur = 0, [], 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    lengths.append(cur)
                n, cur = n + 1, 0
            else:
                cur += len(line.strip().replace("*", ""))
    lengths.append(cur)
    if max_len:
        lengths = [x for x in lengths if 30 <= x <= max_len]
    return len(lengths), sum(lengths)


def duration_s(text):
    """Nextflow trace durations: '1h 2m 3s', '41.5s', '850ms', '-'."""
    if text in ("-", ""):
        return None
    total = 0.0
    for value, unit in re.findall(r"([\d.]+)\s*(ms|d|h|m|s)", text):
        total += float(value) * {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}[unit]
    return total


def proteinannotator(raw, ratio):
    traces = sorted(raw.glob("proteinannotator/results/pipeline_info/execution_trace_*.txt"))
    if not traces:
        return None
    trace = pl.read_csv(traces[-1], separator="\t", infer_schema_length=0)
    trace = trace.filter(pl.col("status") == "COMPLETED").with_columns(
        step=pl.col("name").str.extract(r":([A-Z0-9_]+) \(").fill_null(pl.col("name")),
        realtime_s=pl.col("realtime").map_elements(duration_s, return_dtype=pl.Float64),
        cpu_frac=pl.col("%cpu").str.strip_chars("%").cast(pl.Float64, strict=False) / 100,
    )
    fixed = pl.read_csv(Path(__file__).parent.parent / "fixed_costs_mini10.tsv", separator="\t")
    steps = (
        trace.group_by("step")
        .agg(
            n_tasks=pl.len(),
            pilot_wall_s=pl.col("realtime_s").max(),
            pilot_cpu_h=(pl.col("realtime_s") * pl.col("cpu_frac")).sum() / 3600,
        )
        .join(fixed, on="step", how="left")
        .with_columns(pl.col("fixed_s").fill_null(0.0))
        .with_columns(
            proteome_cpu_h=pl.col("pilot_cpu_h") * ratio,
            # one task over the whole proteome (hmmsearch, s4pred run as one task)
            proteome_wall_h=(pl.col("fixed_s")
                             + (pl.col("pilot_wall_s") - pl.col("fixed_s")).clip(0) * ratio) / 3600,
        )
        .sort("proteome_wall_h", descending=True)
    )
    return steps


def epoch_wall_s(folder):
    start, end = folder / "start_epoch.txt", folder / "end_epoch.txt"
    if not (start.exists() and end.exists()):
        return None
    return int(end.read_text()) - int(start.read_text())


LINE_TS = r"^(\d\d/\d\d/\d\d \d\d:\d\d:\d\d)"


def esmfold(raw, n_proteome, n_pilot):
    log = raw / "esmfold" / "esm-fold.log"
    if not log.exists():
        return None, None
    folds, failed, stamps = [], [], {}
    for line in log.read_text().splitlines():
        ts = re.match(LINE_TS, line)
        if ts and "Loading model" in line:
            stamps["load"] = datetime.strptime(ts.group(1), "%y/%m/%d %H:%M:%S")
        if ts and "Starting Predictions" in line:
            stamps["start"] = datetime.strptime(ts.group(1), "%y/%m/%d %H:%M:%S")
        m = re.search(r"Predicted structure for (\S+) with length (\d+), pLDDT ([\d.]+).* in ([\d.]+)s", line)
        if m:
            folds.append((m.group(1), int(m.group(2)), float(m.group(3)), float(m.group(4))))
        m = re.search(r"Failed \((.*?)\) on sequence (\S+) of length (\d+)", line)
        if m:
            failed.append((m.group(2), int(m.group(3)), m.group(1)))
    per_protein = pl.DataFrame(folds, schema=["protein", "length_aa", "mean_plddt", "fold_s"], orient="row")
    load_s = (stamps["start"] - stamps["load"]).total_seconds() if len(stamps) == 2 else None
    weight = n_proteome / n_pilot
    gpu = (raw / "esmfold" / "gpu.txt").read_text().strip() if (raw / "esmfold" / "gpu.txt").exists() else "?"
    summary = {
        "gpu": gpu,
        "n_folded": per_protein.height,
        "n_failed": len(failed),
        "failed": "; ".join(f"{p} {n} aa ({why})" for p, n, why in failed),
        "model_load_s": load_s,
        "pilot_fold_h": per_protein["fold_s"].sum() / 3600,
        "proteome_fold_gpu_h": per_protein["fold_s"].sum() * weight / 3600,
        # nf-core/proteinfold runs one task per protein, so it loads the model every time
        "proteome_pipeline_gpu_h": (per_protein["fold_s"].sum() * weight
                                    + (load_s or 0) * n_proteome) / 3600,
        "proteome_failed_estimate": round(len(failed) * weight),
    }
    return per_protein, summary


def main():
    raw, lengths_tsv, proteome, out = (Path(a) for a in sys.argv[1:5])
    out.mkdir(parents=True, exist_ok=True)
    pilot = pl.read_csv(lengths_tsv, separator="\t")
    n_pilot, res_pilot = pilot.height, pilot["length_aa"].sum()
    n_prot, res_prot = residues(proteome)
    _, res_prot_pa = residues(proteome, PA_MAX_LEN)
    pilot_pa = pilot.filter(pl.col("length_aa").is_between(30, PA_MAX_LEN))["length_aa"].sum()
    rows = []

    steps = proteinannotator(raw, res_prot_pa / pilot_pa)
    if steps is not None:
        steps.write_csv(out / "proteinannotator_steps.tsv", separator="\t")
        for r in steps.iter_rows(named=True):
            rows.append(("proteinannotator", r["step"], r["n_tasks"], r["pilot_wall_s"] / 3600,
                         r["pilot_cpu_h"], r["proteome_cpu_h"], r["proteome_wall_h"]))

    for name, folder, cpus in [("interproscan_full", raw / "interproscan_full", 16),
                               ("unifire_lite", raw / "unifire", 4)]:
        wall = epoch_wall_s(folder)
        if wall is not None:
            ratio = res_prot / res_pilot
            rows.append((name, name, 1, wall / 3600, wall * cpus / 3600,
                         wall * cpus * ratio / 3600, wall * ratio / 3600))

    per_protein, esm = esmfold(raw, n_prot, n_pilot)
    if per_protein is not None:
        per_protein.write_csv(out / "esmfold_per_protein.tsv", separator="\t")
        pl.DataFrame([esm]).write_csv(out / "esmfold_summary.tsv", separator="\t")
        rows.append(("proteinfold_esmfold", "RUN_ESMFOLD", esm["n_folded"], esm["pilot_fold_h"],
                     None, None, esm["proteome_pipeline_gpu_h"]))

    table = pl.DataFrame(rows, orient="row", schema=[
        "pipeline", "step", "pilot_tasks", "pilot_wall_h", "pilot_cpu_h",
        "proteome_cpu_h", "proteome_wall_h_one_task"])
    table.write_csv(out / "timing_scaled_to_proteome.tsv", separator="\t")
    print(f"proteome {n_prot} proteins {res_prot} residues; pilot {n_pilot} proteins {res_pilot} residues")
    print(f"residue ratio {res_prot / res_pilot:.2f} (proteinannotator, <= {PA_MAX_LEN} aa: {res_prot_pa / pilot_pa:.2f})")
    with pl.Config(tbl_rows=50, tbl_cols=20, fmt_str_lengths=60):
        print(table)
        if esm:
            print(esm)


if __name__ == "__main__":
    main()
