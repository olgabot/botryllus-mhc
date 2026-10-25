"""Figure: where each Botryllus protein's best hit to an MHC protein, B2M or
another C1-set Ig-domain protein ranks among all of that protein's Swiss-Prot
hits.

One dot per Botryllus protein and target class. x = rank by bit score among
the protein's hits in the full search (1 = its best hit of all). A Botryllus
MHC homolog would have its MHC hit at rank 1.

  python plot_mhc_hits.py results/best_hit_per_protein_and_class.tsv results/<name>
"""

import sys
from pathlib import Path

import numpy as np
from matplotlib.ticker import FuncFormatter
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "notebooks"))
import pubfig as pf  # noqa: E402

# Bottom to top on the y axis.
CLASSES = {
    "other C1-set": ("other C1-set\nIg-domain protein", pf.GREY),
    "B2M": ("B2M (MHC class I\nlight chain)", pf.OKABE_ITO["bluish_green"]),
    "MHC class I or I-like": ("MHC class I\nor class I-like", pf.OKABE_ITO["blue"]),
    "MHC class II": ("MHC class II", pf.OKABE_ITO["vermillion"]),
}
TOOLS = ["MMseqs2", "Foldseek (ProstT5 3Di)"]


def main(table: Path, out_stem: str) -> None:
    df = pl.read_csv(table, separator="\t", infer_schema_length=0).with_columns(
        pl.col("rank_by_bits").cast(pl.Int64)
    )
    rng = np.random.default_rng(0)
    pf.use_style()
    fig, axs = pf.figure(pf.TWO_COLUMN_MM, 70, ncols=2, sharey=True)
    xmax = df["rank_by_bits"].max() * 1.5
    for ax, tool, letter in zip(axs, TOOLS, "ab"):
        sub = df.filter(pl.col("tool") == tool)
        ax.axvline(1, color="black", lw=0.5, zorder=0)
        for y, (cls, (label, colour)) in enumerate(CLASSES.items()):
            s = sub.filter(pl.col("target_class") == cls)
            ax.scatter(s["rank_by_bits"], y + rng.uniform(-0.25, 0.25, s.height),
                       s=6, color=colour, alpha=0.7, linewidths=0)
            ax.text(xmax, y + 0.32, f"{s.height} Botryllus protein{'' if s.height == 1 else 's'}", va="bottom",
                    ha="right", fontsize=6)
        ax.set_xscale("log")
        ax.set_xlim(0.7, xmax)
        ax.set_ylim(-0.5, len(CLASSES) - 0.35)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_yticks(range(len(CLASSES)), [v[0] for v in CLASSES.values()])
        ax.set_xlabel("Rank of the best hit in this class among\n"
                      "the Botryllus protein's Swiss-Prot hits (1 = best hit)")
        ax.set_title(tool, loc="left")
        pf.panel_label(ax, letter)
    pf.save(fig, out_stem)


if __name__ == "__main__":
    main(Path(sys.argv[1]), sys.argv[2])
