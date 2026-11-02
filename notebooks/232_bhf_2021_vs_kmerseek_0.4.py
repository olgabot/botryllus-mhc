"""BHF vs human, 2021 sourmash run compared with kmerseek 0.4.0 (Lehninger HP, k=24).

Reads the kmerseek 0.4.0 BHF results (8 arms: scaled 1/2/5/10, low-complexity
k-mers kept/removed) and the two protein FASTA files, then

1. writes the zinc finger (TSHZ1/2/3) and TRAPPC10 rows to a table, and
2. draws the two regions residue by residue, with the E-value of every hit.

Run with the 2025-kmerseek-analysis env:
    python notebooks/232_bhf_2021_vs_kmerseek_0.4.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

import pubfig as pf

KMERSEEK_DIR = Path.home() / "data/botryllus/kmerseek-0.4"
HITS_CSV = KMERSEEK_DIR / "BHF.all_8_arms.kmerseek-0.4.0.hp_lehninger2.k24.csv"
BOTRYLLUS_FA = Path.home() / "data/botryllus/Bs_proteins.fa"
HUMAN_FA = Path.home() / "data/gencode/human/v49/gencode.v49.pc_translations.canonical.fa"

REPO = Path(__file__).resolve().parents[1]
OUT_TABLE = REPO / "results/232_bhf_zinc_finger_trappc10_kmerseek_v040.csv"
OUT_FIG = REPO / "figures/232_bhf_zinc_finger_trappc10_kmerseek_v040"

# Lehninger split, the same as sourmash's aa_to_hp that the 2021 run used.
HYDROPHOBIC = set("AFGILMPVWY")
POLAR = set("CDEHKNQRST")
FOCUS = ["TSHZ1", "TSHZ2", "TSHZ3", "TRAPPC10"]
ARM = "scaled1.lckept"  # every k-mer indexed, low-complexity k-mers kept: closest to 2021


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def hp(seq):
    return "".join("h" if aa in HYDROPHOBIC else "p" if aa in POLAR else "x" for aa in seq)


hits = pl.read_csv(HITS_CSV, infer_schema_length=10000)
bhf = next(s for n, s in read_fasta(BOTRYLLUS_FA) if n.strip() == "BHF")
targets = set(hits.filter(pl.col("gene").is_in(FOCUS))["target_name"])
human = {n: s for n, s in read_fasta(HUMAN_FA) if n in targets}

# Coordinates in the CSV are 0-based, end excluded; check them against the FASTA files.
for r in hits.filter(pl.col("gene").is_in(FOCUS)).iter_rows(named=True):
    assert bhf[r["region_start"] : r["region_end"]] == r["region_subseq"], r["gene"]
    assert human[r["target_name"]][r["target_start"] : r["target_end"]] == r["target_subseq"], r["gene"]

table = (
    hits.filter(pl.col("gene").is_in(FOCUS))
    .select(
        "arm", "gene", "region_start", "region_end", "target_start", "target_end",
        "region_n_shared_kmers", "region_n_mismatches", "region_expected_shared_kmers",
        "region_tail_probability", "query_poisson_pvalue", "region_ka_bits", "region_evalue",
        "region_tfidf", "region_subseq", "target_subseq",
    )
    .sort("arm", "gene", "region_start")
)
OUT_TABLE.parent.mkdir(exist_ok=True)
table.write_csv(OUT_TABLE)
print(table.filter(pl.col("arm") == ARM).drop("region_subseq", "target_subseq"))


def draw_pair(ax, row, motif_positions, title):
    """Four rows per column: BHF residue, BHF HP letter, human HP letter, human residue."""
    q, t = row["region_subseq"], row["target_subseq"]
    qh, th = hp(q), hp(t)
    ys = {"BHF": 3, "BHF HP": 2, "human HP": 1, "human": 0}
    for i, (a, b, ah, bh) in enumerate(zip(q, t, qh, th)):
        if ah != bh:
            ax.add_patch(plt.Rectangle((i - 0.5, 0.5), 1, 2, color=pf.OKABE_ITO["vermillion"], alpha=0.35, lw=0))
        if a == b:
            ax.add_patch(plt.Rectangle((i - 0.5, -0.5), 1, 4, fill=False, ec="black", lw=0.6))
        for txt, y in ((a, 3), (ah, 2), (bh, 1), (b, 0)):
            ax.text(i, y, txt, ha="center", va="center", family="monospace")
    q0, t0 = row["region_start"], row["target_start"]
    for p in motif_positions:
        ax.plot(p - q0, 3.75, marker="v", color="black", ms=3)
    for i, b in enumerate(t):
        if b in "CH":
            ax.plot(i, -0.75, marker="^", color="black", ms=3)
    ax.set_yticks(list(ys.values()), list(ys.keys()))
    ax.set_ylim(-1.2, 4.2)
    ax.set_xlim(-0.8, len(q) - 0.2)
    ticks = [i for i in range(len(q)) if (q0 + i + 1) % 5 == 0]
    ax.set_xticks(ticks, [f"{q0 + i + 1}\n{t0 + i + 1}" for i in ticks])
    ax.tick_params(axis="x", length=2)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(f"Residue number: BHF (top), {title} (bottom)")


pf.use_style()
fig = plt.figure(figsize=(pf.TWO_COLUMN_MM * pf.MM, 125 * pf.MM))
gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1.6], hspace=0.9)
ax_a, ax_b, ax_c = (fig.add_subplot(gs[i]) for i in range(3))

arm = hits.filter(pl.col("arm") == ARM)
zf = arm.filter((pl.col("gene") == "TSHZ2") & (pl.col("region_start") == 9)).row(0, named=True)
tr = arm.filter(pl.col("gene") == "TRAPPC10").row(0, named=True)
# BHF C20, C33, H36, C40 (1-based) are the C-x12-C-x2-H-x3-C residues in the 2021 figure.
draw_pair(ax_a, zf, [19, 32, 35, 39], "TSHZ2")
draw_pair(ax_b, tr, [], "TRAPPC10")

# Legend for panels a and b, placed above panel a.
handles = [
    plt.Rectangle((0, 0), 1, 1, color=pf.OKABE_ITO["vermillion"], alpha=0.35, lw=0),
    plt.Rectangle((0, 0), 1, 1, fill=False, ec="black", lw=0.6),
    plt.Line2D([], [], marker="v", ls="", color="black", ms=3),
    plt.Line2D([], [], marker="^", ls="", color="black", ms=3),
]
labels = [
    "HP letters differ (h = hydrophobic, p = polar)",
    "same amino acid",
    "BHF C or H in the 2021 zinc-finger pattern",
    "human C or H",
]
ax_a.legend(handles, labels, ncol=4, loc="lower left", bbox_to_anchor=(0, 1.15), frameon=False)

# Panel c: E-value of every BHF hit in this arm, focus genes in blue.
d = arm.sort("region_evalue").with_row_index("rank")
is_focus = d["gene"].is_in(FOCUS).to_numpy()
ax_c.scatter(d["rank"].to_numpy()[~is_focus], d["region_evalue"].to_numpy()[~is_focus], s=6, color=pf.GREY, label="other human proteins")
ax_c.scatter(d["rank"].to_numpy()[is_focus], d["region_evalue"].to_numpy()[is_focus], s=10, color=pf.OKABE_ITO["blue"], label="TSHZ1, TSHZ2, TSHZ3, TRAPPC10")
ax_c.axhline(1, color="black", lw=0.5, ls="--", label="E-value = 1")
for r in d.filter(pl.col("gene").is_in(FOCUS + ["DHX8"])).iter_rows(named=True):
    label = r["gene"] if r["gene"] != "TSHZ2" else f'TSHZ2 (BHF {r["region_start"] + 1}-{r["region_end"]})'
    ax_c.annotate(label, (r["rank"], r["region_evalue"]), xytext=(2, -9) if r["region_start"] == 9 and r["gene"] == "TSHZ2" else (2, 4), textcoords="offset points", fontsize=5)
ax_c.set_yscale("log")
ax_c.set_yticks([1, 3, 10, 30, 100, 300], ["1", "3", "10", "30", "100", "300"])
ax_c.set_ylim(0.7, 500)
ax_c.set_xlabel("Human hits to BHF, ranked by E-value")
ax_c.set_ylabel("E-value")
ax_c.legend(loc="upper left", frameon=False)

for ax, letter in zip((ax_a, ax_b, ax_c), "abc"):
    pf.panel_label(ax, letter)
OUT_FIG.parent.mkdir(exist_ok=True)
pf.save(fig, str(OUT_FIG))
