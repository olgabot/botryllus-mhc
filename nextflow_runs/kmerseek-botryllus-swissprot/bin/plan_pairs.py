#!/usr/bin/env python3
"""Write assets/pairs.tsv: for each of kmerseek's 19 alphabets, the k-mer sizes that carry
at least the bits of protein20 at k = 6 and at k = 7.

The rule is notebook 235's (notebooks/235_plan_arms.py, branch
nb232-kmerseek-botryllus-human-immune-hits), which follows notebook 274 of
seanome/2024-kmerseek-analysis:

* Bits per position = -log2(sum of q_c squared), where q_c is the share of Swiss-Prot
  residues in letter class c. sum q_c^2 is the chance that two residues drawn at random
  fall in the same class, so these bits say how rarely two unrelated k-mers match.
* Swiss-Prot composition: UniProtKB/Swiss-Prot release 2026_03, "6.1 Composition in
  percent for the complete database" (https://web.expasy.org/docs/relnotes/relstat.html).
* For each alphabet, k is the smallest k whose bits reach protein20's bits at the
  reference k (6 x 4.0639 = 24.38 bits, 7 x 4.0639 = 28.45 bits).

Class partitions are kmerseek v0.4.0's (src/rust/alphabets.rs); the extension mismatch
penalty and X-drop per alphabet are the ones notebooks 233-235 ran (notebook 241's plan).
Standard library only.

Usage:
    python3 bin/plan_pairs.py            # writes assets/pairs.tsv next to bin/
"""

from __future__ import annotations

import math
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "assets" / "pairs.tsv"
REFERENCE_ALPHABET, REFERENCE_KS = "protein20", (6, 7)

SWISSPROT_PERCENT = {
    "A": 8.25, "Q": 3.93, "L": 9.64, "S": 6.66,
    "R": 5.52, "E": 6.71, "K": 5.79, "T": 5.36,
    "N": 4.06, "G": 7.07, "M": 2.41, "W": 1.10,
    "D": 5.46, "H": 2.27, "F": 3.86, "Y": 2.92,
    "C": 1.38, "I": 5.90, "P": 4.75, "V": 6.85,
}  # fmt: skip

CLUSTERS: dict[str, list[str]] = {
    "protein20": list("ACDEFGHIKLMNPQRSTVWY"),
    "dayhoff6": ["C", "AGPST", "DENQ", "HKR", "ILMV", "FWY"],
    "hp_lehninger2": ["AFGILMPVWY", "CDEHKNQRST"],
    "hp_thomas_dill2": ["ACFILMVWY", "DEGHKNPQRST"],
    "hp_kyte_doolittle2": ["ACFILMV", "DEGHKNPQRSTWY"],
    "hp_thomas_dill_no_c2": ["AFILMVWY", "CDEGHKNPQRST"],
    "hp_lehninger_c_nonpolar2": ["ACFGILMPVWY", "DEHKNQRST"],
    "hp_lehninger_hpc3": ["AFGILMPVWY", "DEHKNQRST", "C"],
    "hp_pbotc_1st_ed2": ["ACFILMPVWY", "DEGHKNQRST"],
    "gbmr4": ["ADKERNTSQ", "YFLIVMCWH", "G", "P"],
    "polarity4": ["GAVLIFWMP", "STCYNQ", "DE", "HKR"],
    "wwmj5": ["CMFILVWY", "ATH", "GP", "DE", "SNQRK"],
    "gbmr7": ["DN", "AEFIKLMQRVWY", "CH", "T", "S", "G", "P"],
    "funcgroups8": ["GVALI", "ST", "CM", "FY", "WHP", "NQ", "DE", "KR"],
    "sdm12": ["A", "D", "KER", "N", "TSQ", "YF", "LIVM", "C", "W", "H", "G", "P"],
    "mmseqs12": ["AST", "LM", "IV", "KR", "EQ", "ND", "FY", "C", "G", "H", "P", "W"],
    "wass14": ["WM", "DI", "P", "C", "AV", "K", "T", "RE", "G", "L", "Y", "SH", "F", "NQ"],
    "hsdm17": ["A", "D", "KE", "R", "N", "T", "S", "Q", "Y", "F", "LIV", "M", "C", "W", "H", "G", "P"],
    "uniprot18": ["A", "R", "N", "D", "C", "Q", "EP", "G", "HL", "I", "K", "M", "F", "S", "T", "W", "Y", "V"],
}  # fmt: skip

# alphabet -> (extension mismatch penalty, X-drop)
EXTENSION: dict[str, tuple[str, str]] = {
    "protein20": ("0.14", "0.57"), "dayhoff6": ("0.41", "1.65"),
    "hp_lehninger2": ("1.51", "6.03"), "hp_thomas_dill2": ("1.63", "6.53"),
    "hp_kyte_doolittle2": ("1.57", "6.28"), "hp_thomas_dill_no_c2": ("1.62", "6.47"),
    "hp_lehninger_c_nonpolar2": ("1.52", "6.08"), "hp_lehninger_hpc3": ("0.87", "3.46"),
    "hp_pbotc_1st_ed2": ("1.59", "6.35"), "gbmr4": ("0.67", "2.7"),
    "polarity4": ("0.54", "2.17"), "wwmj5": ("0.47", "1.89"), "gbmr7": ("0.32", "1.27"),
    "funcgroups8": ("2", "8"), "sdm12": ("0.24", "0.94"), "mmseqs12": ("0.22", "0.89"),
    "wass14": ("0.17", "0.7"), "hsdm17": ("0.18", "0.7"), "uniprot18": ("0.15", "0.6"),
}  # fmt: skip


def swissprot_bits(clusters: list[str]) -> float:
    """Bits per position: -log2 of the chance two Swiss-Prot residues share a class."""
    total = sum(SWISSPROT_PERCENT.values())
    shares = [sum(SWISSPROT_PERCENT[aa] for aa in c) / total for c in clusters]
    return -math.log2(sum(q * q for q in shares))


def main() -> None:
    bits = {name: swissprot_bits(cl) for name, cl in CLUSTERS.items()}
    rows = ["alphabet\tn_classes\tbits_per_position\tprotein20_k\tksize\tbits_per_kmer\tpenalty\txdrop"]
    for ref_k in REFERENCE_KS:
        target = ref_k * bits[REFERENCE_ALPHABET]
        for name, cl in CLUSTERS.items():
            # the 1e-9 keeps protein20 at the reference k, whose quotient is exact
            k = math.ceil(target / bits[name] - 1e-9)
            penalty, xdrop = EXTENSION[name]
            rows.append(f"{name}\t{len(cl)}\t{bits[name]:.4f}\t{ref_k}\t{k}\t{k * bits[name]:.2f}\t{penalty}\t{xdrop}")
    OUT.write_text("\n".join(rows) + "\n")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
