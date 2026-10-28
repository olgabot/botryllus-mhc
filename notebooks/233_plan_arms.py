#!/usr/bin/env python3
"""Write notebooks/233_arms.tsv: one kmerseek arm per alphabet, each carrying at least
the bits of hp_lehninger2 at k = 24.

Bits per k-mer = k x H, where H is the Shannon entropy of the alphabet's class shares
measured on the human target proteome (not log2 of the class count: gbmr7 has 7 classes
but only 2.0 bits per position, because one class holds most residues). Notebook 241 in
seanome/2024-kmerseek-analysis uses the same measure. For each alphabet k is the smallest
k whose bits reach hp_lehninger2 k = 24, so no arm is easier to hit by chance than the HP
arm of notebook 232.

The mismatch penalty and X-drop of the region extension are each alphabet's own, from
its measured copy rate (kappa_by_alphabet.tsv of the dark-set pipeline, the values
notebook 241 ran with; funcgroups8 has no measured kappa and takes the fallback 2 / 8).

Usage:
    python3 notebooks/233_plan_arms.py [human.fa]
"""

from __future__ import annotations

import collections
import math
import sys
from pathlib import Path

HUMAN = Path(
    sys.argv[1]
    if len(sys.argv) > 1
    else "/Users/olga/data/gencode/human/v49/gencode.v49.pc_translations.canonical.fa"
)
OUT = Path(__file__).with_name("233_arms.tsv")

REFERENCE_ALPHABET, REFERENCE_K = "hp_lehninger2", 24

# Class partitions, copied from kmerseek v0.4.0 src/rust/alphabets.rs (protein20 and
# dayhoff6 from the sourmash encoder it uses).
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
}

# alphabet -> (mismatch penalty, X-drop), as notebook 241's plan.json ran them.
EXTENSION: dict[str, tuple[str, str]] = {
    "protein20": ("0.14", "0.57"),
    "dayhoff6": ("0.41", "1.65"),
    "hp_lehninger2": ("1.51", "6.03"),
    "hp_thomas_dill2": ("1.63", "6.53"),
    "hp_kyte_doolittle2": ("1.57", "6.28"),
    "hp_thomas_dill_no_c2": ("1.62", "6.47"),
    "hp_lehninger_c_nonpolar2": ("1.52", "6.08"),
    "hp_lehninger_hpc3": ("0.87", "3.46"),
    "hp_pbotc_1st_ed2": ("1.59", "6.35"),
    "gbmr4": ("0.67", "2.7"),
    "polarity4": ("0.54", "2.17"),
    "wwmj5": ("0.47", "1.89"),
    "gbmr7": ("0.32", "1.27"),
    "funcgroups8": ("2", "8"),
    "sdm12": ("0.24", "0.94"),
    "mmseqs12": ("0.22", "0.89"),
    "wass14": ("0.17", "0.7"),
    "hsdm17": ("0.18", "0.7"),
    "uniprot18": ("0.15", "0.6"),
}


def residue_counts(fasta: Path) -> collections.Counter:
    """Counts of the 20 standard residues in a protein FASTA."""
    counts: collections.Counter = collections.Counter()
    with open(fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                counts.update(line.strip().upper())
    return collections.Counter({aa: counts[aa] for aa in "ACDEFGHIKLMNPQRSTVWY"})


def entropy(clusters: list[str], counts: collections.Counter) -> float:
    """Bits per position: Shannon entropy of the class shares."""
    total = sum(counts.values())
    shares = [sum(counts[aa] for aa in c) / total for c in clusters]
    return -sum(p * math.log2(p) for p in shares if p > 0)


def main() -> None:
    counts = residue_counts(HUMAN)
    for name, cl in CLUSTERS.items():
        assert sorted("".join(cl)) == sorted("ACDEFGHIKLMNPQRSTVWY"), name
    h = {name: entropy(cl, counts) for name, cl in CLUSTERS.items()}
    reference_bits = REFERENCE_K * h[REFERENCE_ALPHABET]
    rows = ["alphabet\tn_classes\tbits_per_position\tksize\tbits_per_kmer\tpenalty\txdrop"]
    for name, cl in CLUSTERS.items():
        # the 1e-9 keeps k = 24 for the reference itself, whose quotient is exactly 24
        k = math.ceil(reference_bits / h[name] - 1e-9)
        penalty, xdrop = EXTENSION[name]
        rows.append(f"{name}\t{len(cl)}\t{h[name]:.4f}\t{k}\t{k * h[name]:.2f}\t{penalty}\t{xdrop}")
    OUT.write_text("\n".join(rows) + "\n")
    print(f"reference: {REFERENCE_ALPHABET} k = {REFERENCE_K} = {reference_bits:.2f} bits")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
