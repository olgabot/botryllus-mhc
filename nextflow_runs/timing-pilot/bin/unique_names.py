"""Give every protein in a FASTA its own name, numbering isoforms .1, .2, ...

The Botryllus proteome reuses 2,019 names across 4,613 records (up to 6 records share a
name, and they are different proteins: isoforms of one gene model). nf-core/proteinannotator
removes duplicates by name, so it would keep one record per name and drop 2,594 proteins
without an error. Every record of a shared name gets an isoform number in file order
(FUN039215_FUN039215.1, FUN039215_FUN039215.2); a name used once is left as it is. Stop
codons (*) are removed, because InterProScan rejects them.

Usage:
    python3 unique_names.py <in.fasta> <out.fasta> <renamed.tsv>
"""

import sys
from collections import Counter


def names(path):
    with open(path) as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">")]


def main():
    fasta, out_fasta, out_tsv = sys.argv[1:4]
    total = Counter(names(fasta))
    shared = {name for name, n in total.items() if n > 1}
    clash = {f"{name}.{i}" for name in shared for i in range(1, total[name] + 1)} & set(total)
    if clash:
        sys.exit(f"isoform names already used in the input: {sorted(clash)[:5]}")
    seen, renamed = Counter(), []
    with open(fasta) as fh, open(out_fasta, "w") as out:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seen[name] += 1
                new = f"{name}.{seen[name]}" if name in shared else name
                if new != name:
                    renamed.append((name, new))
                out.write(f">{new}\n")
            else:
                out.write(line.replace("*", ""))
    with open(out_tsv, "w") as tsv:
        tsv.write("original\tunique\n")
        for old, new in renamed:
            tsv.write(f"{old}\t{new}\n")
    print(f"{total.total()} proteins, {len(renamed)} renamed ({len(shared)} names were shared)")


if __name__ == "__main__":
    main()
