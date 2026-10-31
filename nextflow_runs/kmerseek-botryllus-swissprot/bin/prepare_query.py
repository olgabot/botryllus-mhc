#!/usr/bin/env python3
"""

from __future__ import annotationsMake the Botryllus proteome safe to search with kmerseek and to read as CSV, and split it
into chunks.

Input is the renamed proteome (botryllus_new_genome_protein_unique_names.fasta, 45,339
proteins), where every protein already has its own name. Two changes, both recorded in the
renaming table so a match can be traced back to the input record (the same rules as
notebook 234's 234_prepare_fun_proteome.py):

* 10 headers contain a comma (FUN001900_ACTIN,). kmerseek writes the header into the
  query_name column of a CSV, so a comma would split that column. Commas become '-'.
* One sequence (FUN029508_FUN029508) has a '.' inside it. It becomes 'X', the letter for
  an unknown residue.

The script stops if two records share a name, so the original FASTA (2,019 repeated names)
cannot be used by mistake.

Usage: prepare_query.py <in.fasta> <chunk size> <out dir>
Writes <out dir>/chunk_0000.fasta, chunk_0001.fasta, ... and <out dir>/query_renaming.tsv.
"""

import sys
from collections import Counter
from pathlib import Path


def read_fasta(path: str) -> list[tuple[str, str]]:
    """Return (header without '>', sequence) pairs in file order."""
    records: list[tuple[str, str]] = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                records.append((line[1:].strip(), ""))
            elif records:
                records[-1] = (records[-1][0], records[-1][1] + line.strip())
    return records


def main(in_path: str, chunk_size: int, out_dir: Path) -> None:
    records = read_fasta(in_path)
    repeated = [name for name, n in Counter(h for h, _ in records).items() if n > 1]
    if repeated:
        sys.exit(f"{len(repeated)} names occur more than once (e.g. {repeated[0]}): "
                 "use the renamed proteome, not the original FASTA")
    out_dir.mkdir(parents=True, exist_ok=True)
    new_names: set[str] = set()
    with open(out_dir / "query_renaming.tsv", "w") as renaming:
        renaming.write("original_header\tnew_header\tchanges\n")
        for i in range(0, len(records), chunk_size):
            with open(out_dir / f"chunk_{i // chunk_size:04d}.fasta", "w") as out:
                for header, sequence in records[i:i + chunk_size]:
                    changes = []
                    new_header = header.replace(",", "-")
                    if new_header != header:
                        changes.append("comma")
                    if "." in sequence:
                        sequence = sequence.replace(".", "X")
                        changes.append("dot_to_X")
                    if changes:
                        renaming.write(f"{header}\t{new_header}\t{','.join(changes)}\n")
                    new_names.add(new_header)
                    out.write(f">{new_header}\n")
                    for start in range(0, len(sequence), 60):
                        out.write(sequence[start:start + 60] + "\n")
    if len(new_names) != len(records):
        sys.exit("replacing commas made two names equal")
    print(f"{len(records)} proteins in {-(-len(records) // chunk_size)} chunks of up to {chunk_size}")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), Path(sys.argv[3]))
