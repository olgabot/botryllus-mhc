#!/usr/bin/env python3
"""Write notebooks/235_arms.tsv: one kmerseek arm per alphabet, each carrying at least
the bits of protein20 at k = 7, with bits measured from Swiss-Prot residue frequencies.

Notebook 233 set the bar at hp_lehninger2 k = 24 (24.0 bits, Shannon entropy of the
human proteome's class shares), which put protein20 at k = 6. This table raises the bar
to protein20 k = 7 and changes how bits are counted, to match notebook 274 in
seanome/2024-kmerseek-analysis:

* Bits per position = -log2(sum of q_c squared), where q_c is the share of Swiss-Prot
  residues in letter class c. sum q_c^2 is the chance that two residues drawn at random
  fall in the same class, so these bits say how rarely two unrelated k-mers match.
* Swiss-Prot composition: UniProtKB/Swiss-Prot release 2026_03, "6.1 Composition in
  percent for the complete database" (https://web.expasy.org/docs/relnotes/relstat.html),
  the numbers notebook 274 used.
* For each alphabet k is the smallest k whose bits reach 7 x protein20's bits (28.45).

Class partitions and extension settings are notebook 233's (233_plan_arms.py).

Usage:
    python3 notebooks/235_plan_arms.py
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "235_arms.tsv"

REFERENCE_ALPHABET, REFERENCE_K = "protein20", 7

SWISSPROT_PERCENT = {
    "A": 8.25, "Q": 3.93, "L": 9.64, "S": 6.66,
    "R": 5.52, "E": 6.71, "K": 5.79, "T": 5.36,
    "N": 4.06, "G": 7.07, "M": 2.41, "W": 1.10,
    "D": 5.46, "H": 2.27, "F": 3.86, "Y": 2.92,
    "C": 1.38, "I": 5.90, "P": 4.75, "V": 6.85,
}  # fmt: skip


def load_plan_233():
    """Notebook 233's planning module, for its CLUSTERS and EXTENSION tables."""
    spec = importlib.util.spec_from_file_location("plan_233", HERE / "233_plan_arms.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def swissprot_bits(clusters: list[str]) -> float:
    """Bits per position: -log2 of the chance two Swiss-Prot residues share a class."""
    total = sum(SWISSPROT_PERCENT.values())
    shares = [sum(SWISSPROT_PERCENT[aa] for aa in c) / total for c in clusters]
    return -math.log2(sum(q * q for q in shares))


def main() -> None:
    plan = load_plan_233()
    bits = {name: swissprot_bits(cl) for name, cl in plan.CLUSTERS.items()}
    reference_bits = REFERENCE_K * bits[REFERENCE_ALPHABET]
    rows = ["alphabet\tn_classes\tbits_per_position\tksize\tbits_per_kmer\tpenalty\txdrop"]
    for name, cl in plan.CLUSTERS.items():
        # the 1e-9 keeps k = 7 for protein20 itself, whose quotient is exactly 7
        k = math.ceil(reference_bits / bits[name] - 1e-9)
        penalty, xdrop = plan.EXTENSION[name]
        rows.append(f"{name}\t{len(cl)}\t{bits[name]:.4f}\t{k}\t{k * bits[name]:.2f}\t{penalty}\t{xdrop}")
    OUT.write_text("\n".join(rows) + "\n")
    print(f"reference: {REFERENCE_ALPHABET} k = {REFERENCE_K} = {reference_bits:.2f} bits (Swiss-Prot 2026_03)")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
