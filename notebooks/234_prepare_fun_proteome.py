"""Make the FUN Botryllus proteome safe to search with kmerseek and to filter as CSV.

Input: botryllus_new_genome_protein.fasta (45,339 proteins, 23 Oct 2026). Three changes,
each recorded in the renaming table so a match can be traced to the original record:

* 2,019 IDs occur 2 to 6 times with different sequences (isoforms whose transcript suffix
  was lost). Each copy of a repeated ID gets -T1, -T2, ... in file order; IDs that occur
  once keep their name.
* 10 headers contain a comma (FUN001900_ACTIN,), which would split the query_name column
  of kmerseek's CSV output. Commas become '-'.
* One sequence (FUN029508_FUN029508) has a '.' inside it. It becomes 'X', the letter for
  an unknown residue that the old proteome (Bs_proteins.fa) already uses.

Usage: python3 234_prepare_fun_proteome.py <in.fasta> <out.fasta> <renaming.tsv>
"""

import sys
from collections import Counter


def read_fasta(path: str) -> list[tuple[str, str]]:
    """Return (header without '>', sequence) pairs in file order."""
    records: list[tuple[str, str]] = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                records.append((line[1:], ""))
            elif records:
                records[-1] = (records[-1][0], records[-1][1] + line.strip())
    return records


def main(in_path: str, out_path: str, renaming_path: str) -> None:
    records = read_fasta(in_path)
    n_copies = Counter(header for header, _ in records)
    seen: Counter[str] = Counter()
    n_renamed = n_comma = n_dot = 0
    with open(out_path, "w") as out, open(renaming_path, "w") as renaming:
        renaming.write("original_header\tnew_header\tchanges\n")
        for header, sequence in records:
            changes = []
            new_header = header.replace(",", "-")
            if new_header != header:
                changes.append("comma")
                n_comma += 1
            if n_copies[header] > 1:
                seen[header] += 1
                new_header = f"{new_header}-T{seen[header]}"
                changes.append("repeated_id")
                n_renamed += 1
            if "." in sequence:
                sequence = sequence.replace(".", "X")
                changes.append("dot_to_X")
                n_dot += 1
            if changes:
                renaming.write(f"{header}\t{new_header}\t{','.join(changes)}\n")
            out.write(f">{new_header}\n")
            for start in range(0, len(sequence), 60):
                out.write(sequence[start:start + 60] + "\n")
    print(f"{len(records)} proteins; {n_renamed} copies of repeated IDs renamed, "
          f"{n_comma} headers with a comma, {n_dot} sequences with '.'")


if __name__ == "__main__":
    main(*sys.argv[1:4])
