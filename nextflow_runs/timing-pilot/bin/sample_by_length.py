"""Pick proteins whose lengths follow the whole proteome's length spread.

The proteome is sorted by length and cut into N_BINS bins of equal protein count
(deciles for 10 bins). The same number of proteins is drawn at random from each bin,
so the sample has the proteome's length distribution, its longest proteins included.
Stop codons (*) are removed, because InterProScan rejects them.

Usage:
    python3 sample_by_length.py <proteome.fasta> <n> <out.fasta> <out_lengths.tsv>
"""

import random
import sys

N_BINS = 10
SEED = 20261027


def read_fasta(path):
    records = []
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq).replace("*", "")))
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        records.append((name, "".join(seq).replace("*", "")))
    return records


def main():
    fasta, n, out_fasta, out_tsv = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
    records = sorted(read_fasta(fasta), key=lambda r: len(r[1]))
    if n % N_BINS:
        sys.exit(f"n={n} must be a multiple of {N_BINS}")
    rng = random.Random(SEED)
    per_bin = n // N_BINS
    picked = []
    for b in range(N_BINS):
        lo = b * len(records) // N_BINS
        hi = (b + 1) * len(records) // N_BINS
        for name, seq in rng.sample(records[lo:hi], per_bin):
            picked.append((b, name, seq))

    with open(out_fasta, "w") as fa, open(out_tsv, "w") as tsv:
        tsv.write("protein\tlength_aa\tlength_bin\n")
        for b, name, seq in picked:
            fa.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fa.write(seq[i : i + 60] + "\n")
            tsv.write(f"{name}\t{len(seq)}\t{b}\n")

    total = sum(len(s) for _, s in records)
    sample = sum(len(s) for _, _, s in picked)
    print(f"proteome: {len(records)} proteins, {total} residues")
    print(f"sample:   {len(picked)} proteins, {sample} residues")


if __name__ == "__main__":
    main()
