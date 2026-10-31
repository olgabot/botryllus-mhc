#!/usr/bin/env python3
"""Figure: the k-mer size searched for each alphabet, against the range notebook 274 of
seanome/2024-kmerseek-analysis says each alphabet can use.

Reads assets/pairs.tsv (bin/plan_pairs.py) and assets/kmerseek_ksize_range_per_alphabet.tsv
(k_min and k_max columns of tables/274_ksizes_to_test_per_alphabet_human_swissprot.csv).
Writes figures/kmerseek_botryllus_swissprot_ksizes.{pdf,svg,png}.

Usage: python bin/plot_pairs.py   (from this folder, in the 2025-kmerseek-analysis env)
"""

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pubfig as pf  # noqa: E402

pairs = list(csv.DictReader(open(HERE.parent / "assets/pairs.tsv"), delimiter="\t"))
ranges = {r["alphabet"]: r for r in csv.DictReader(
    open(HERE.parent / "assets/kmerseek_ksize_range_per_alphabet.tsv"), delimiter="\t")}

# alphabets ordered by bits per letter, most bits at the top
alphabets = sorted({p["alphabet"]: float(p["bits_per_position"]) for p in pairs}.items(),
                   key=lambda x: x[1])
order = [a for a, _ in alphabets]
y = {a: i for i, a in enumerate(order)}

pf.use_style()
fig, ax = pf.figure(pf.ONE_COLUMN_MM, 80)
for a in order:
    r = ranges[a]
    ax.plot([int(r["k_min"]), int(r["k_max"])], [y[a], y[a]], color=pf.GREY, lw=3,
            solid_capstyle="butt", alpha=0.5, zorder=1,
            label="k-mer sizes notebook 274 can use" if a == order[0] else None)
for ref, colour, marker in (("6", pf.OKABE_ITO["blue"], "o"), ("7", pf.OKABE_ITO["vermillion"], "s")):
    sel = [p for p in pairs if p["protein20_k"] == ref]
    bits = sel[0]["bits_per_kmer"] if sel[0]["alphabet"] == "protein20" else ""
    ax.scatter([int(p["ksize"]) for p in sel], [y[p["alphabet"]] for p in sel], color=colour,
               marker=marker, s=12, zorder=2,
               label=f"searched: as many bits as protein20 k = {ref} ({bits} bits)")
ax.set_yticks(range(len(order)), [f"{a} ({b:.2f})" for a, b in alphabets])
ax.set_ylabel("Alphabet (bits per letter)")
ax.set_xlabel("k-mer size")
ax.set_xlim(0, 41)
pf.shared_legend(fig, ncol=1)
out = HERE.parent / "figures"
out.mkdir(exist_ok=True)
pf.save(fig, str(out / "kmerseek_botryllus_swissprot_ksizes"))
