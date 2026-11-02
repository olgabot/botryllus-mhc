"""BHF vs human: kmerseek main (fd549f9) against release-0.4.0 (0aed5ca), same settings.

Inputs are the CSVs from scripts/232_bhf_kmerseek_main_vs_v040.sbatch, copied from
Sherlock to results/233_sherlock/, plus the September 0.4.0 results for BHF.

1. Control: the v0.4.0 rerun must give back the September BHF rows.
2. Per gene, the best E-value under each build, for both low-complexity arms.
3. A figure: percent change in each region's E-value from 0.4.0 to main.

    python notebooks/233_bhf_kmerseek_main_vs_v040.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

import pubfig as pf

REPO = Path(__file__).resolve().parents[1]
RUN_DIR = REPO / "results/233_sherlock"
SEPT_DIR = Path.home() / "data/botryllus/kmerseek-0.4"
OUT_TABLE = REPO / "results/233_bhf_kmerseek_main_vs_v040.csv"
OUT_FIG = REPO / "figures/233_bhf_kmerseek_main_vs_v040"

FIG_2021 = ["ZNF292", "RSF1", "TSHZ1", "TSHZ2", "TSHZ3", "RNMT", "SFI1", "TRAPPC10", "NDNF"]
HITS_2021 = FIG_2021 + [
    "ADGB", "CCDC96", "DHX8", "HMGB3", "IFT122", "JMJD1C", "MYRFL", "ORC2", "PDIA5", "PNN",
    "SLC16A7", "TNPO2", "ZNF713",
]
ARMS = ["lckept", "lcremoved"]
KEY = ["target_name", "region_start", "region_end"]


def load(path):
    df = pl.read_csv(path, infer_schema_length=100000)
    df = df.filter(pl.col("query_name").str.strip_chars() == "BHF")
    return df.with_columns(pl.col("target_name").str.split("|").list.get(6).alias("gene"))


# 1. Control: v0.4.0 rerun on Sherlock (BHF as the only query) against September.
for arm in ARMS:
    sept = load(SEPT_DIR / f"BHF.hp_lehninger2.k24.scaled1.{arm}.v0.4.0.query_subset.csv")
    rerun = load(RUN_DIR / f"v040.{arm}.csv")
    joined = sept.select(KEY + ["region_evalue"]).join(
        rerun.select(KEY + ["region_evalue"]), on=KEY, how="full", suffix="_rerun", coalesce=True
    )
    n_both = joined.filter(pl.col("region_evalue").is_not_null() & pl.col("region_evalue_rerun").is_not_null()).height
    max_rel = (
        joined.filter(pl.col("region_evalue").is_finite() & pl.col("region_evalue_rerun").is_finite())
        .select(((pl.col("region_evalue") - pl.col("region_evalue_rerun")).abs() / pl.col("region_evalue")).max())
        .item()
    )
    print(f"control {arm}: September {sept.height} rows, rerun {rerun.height}, in both {n_both}, "
          f"largest relative E-value difference {max_rel}")

# 2. Best E-value per gene under each build.
rows = []
for arm in ARMS:
    for build in ["v040", "main"]:
        d = load(RUN_DIR / f"{build}.{arm}.csv")
        print(f"{build} {arm}: {d.height} rows, {d['gene'].n_unique()} genes, "
              f"min E {d['region_evalue'].min():.3g}, E < 1: {(d['region_evalue'] < 1).sum()}")
        rows.append(d.group_by("gene").agg(pl.col("region_evalue").min().alias("best_evalue"))
                    .with_columns(pl.lit(arm).alias("arm"), pl.lit(build).alias("build")))
long = pl.concat(rows)
wide = long.pivot(on="build", index=["arm", "gene"], values="best_evalue").with_columns(
    pl.col("gene").is_in(HITS_2021).alias("found_2021"),
    pl.col("gene").is_in(FIG_2021).alias("in_2021_figure"),
).sort("arm", "main", nulls_last=True)
OUT_TABLE.parent.mkdir(exist_ok=True)
wide.write_csv(OUT_TABLE)
pl.Config.set_tbl_rows(200)
print(wide.filter(pl.col("found_2021") | (pl.col("main") < 1) | (pl.col("v040") < 1)))

# 3. Figure: per matched region, the change in E-value from 0.4.0 to main.
pf.use_style()
fig, axs = plt.subplots(1, 2, figsize=(pf.TWO_COLUMN_MM * pf.MM, 70 * pf.MM), sharey=True)
groups = [
    ("not found in 2021", ~pl.col("gene").is_in(HITS_2021), pf.GREY, 6),
    ("found in 2021, not in the 2021 figure", pl.col("gene").is_in(HITS_2021) & ~pl.col("gene").is_in(FIG_2021), pf.OKABE_ITO["orange"], 10),
    ("in the 2021 figure", pl.col("gene").is_in(FIG_2021), pf.OKABE_ITO["blue"], 10),
]
for ax, arm in zip(axs, ARMS):
    a = load(RUN_DIR / f"v040.{arm}.csv").select(KEY + ["gene", "region_evalue"])
    b = load(RUN_DIR / f"main.{arm}.csv").select(KEY + ["region_evalue"])
    d = a.join(b, on=KEY, suffix="_main").with_columns(
        ((pl.col("region_evalue_main") / pl.col("region_evalue") - 1) * 100).alias("pct")
    )
    assert d.height == a.height == b.height, "main and 0.4.0 matched different regions"
    for label, mask, colour, size in groups:
        sub = d.filter(mask)
        ax.scatter(sub["region_evalue"], sub["pct"], s=size, color=colour, label=label, zorder=3)
    for r in d.filter(pl.col("gene") == "DHX8").iter_rows(named=True):
        ax.annotate("DHX8", (r["region_evalue"], r["pct"]), xytext=(3, 3), textcoords="offset points", fontsize=5)
    ax.axhline(0, color="black", lw=0.5, label="no change")
    ax.axvline(1, color="black", lw=0.5, ls="--", label="E-value = 1")
    ax.set_xscale("log")
    ax.set_xticks([0.3, 1, 3, 10, 30, 100, 300], ["0.3", "1", "3", "10", "30", "100", "300"])
    ax.set_xlim(0.3, 600)
    ax.set_ylim(-1.2, 0.2)
    ax.set_xlabel("E-value, kmerseek 0.4.0")
    ax.text(0.98, 0.04, "low-complexity k-mers " + ("kept" if arm == "lckept" else "removed"),
            transform=ax.transAxes, fontsize=6, ha="right")
axs[0].set_ylabel("Change in E-value, main vs 0.4.0 (%)")
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, ncol=5, loc="upper center", frameon=False)
fig.subplots_adjust(top=0.85, wspace=0.12)
for ax, letter in zip(axs, "ab"):
    pf.panel_label(ax, letter)
OUT_FIG.parent.mkdir(exist_ok=True)
pf.save(fig, str(OUT_FIG))
