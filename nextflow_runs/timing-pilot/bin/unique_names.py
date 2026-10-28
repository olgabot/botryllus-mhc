"""Give every protein in a FASTA its own name.

The Botryllus proteome reuses 2,019 names across 4,613 records (up to 6 records share a
name, and they are different proteins). nf-core/proteinannotator removes duplicates by
name, so it would keep one record per name and drop 2,594 proteins without an error. The
first record keeps its name; later ones get _2, _3, ... in file order. Stop codons (*) are
removed, because InterProScan rejects them.

Usage:
    python3 unique_names.py <in.fasta> <out.fasta> <renamed.tsv>
"""

import sys


def main():
    fasta, out_fasta, out_tsv = sys.argv[1:4]
    seen, renamed = {}, []
    with open(fasta) as fh, open(out_fasta, "w") as out:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seen[name] = seen.get(name, 0) + 1
                new = name if seen[name] == 1 else f"{name}_{seen[name]}"
                if new != name:
                    renamed.append((name, new))
                out.write(f">{new}\n")
            else:
                out.write(line.replace("*", ""))
    clash = set(seen) & {new for _, new in renamed}
    if clash:
        sys.exit(f"new names already used in the input: {sorted(clash)[:5]}")
    with open(out_tsv, "w") as tsv:
        tsv.write("original\tunique\n")
        for old, new in renamed:
            tsv.write(f"{old}\t{new}\n")
    print(f"{sum(seen.values())} proteins, {len(renamed)} renamed "
          f"({sum(1 for v in seen.values() if v > 1)} names were shared)")


if __name__ == "__main__":
    main()
