#!/usr/bin/env python3
"""

from __future__ import annotationsRewrite each UniProtKB/Swiss-Prot FASTA header to its bare accession, and save the rest
of the header as a table.

kmerseek writes the whole FASTA header into the target_name column of its CSV output, and
Swiss-Prot descriptions contain commas, which would split that column. With the accession
as the whole header, target_name is the accession, and the table gives the entry name,
protein name, organism and taxon ID for each one.

Usage: prepare_swissprot.py <uniprot_sprot.fasta.gz> <out.fasta> <out.metadata.tsv>
"""

import gzip
import re
import sys

HEADER = re.compile(r"^>sp\|(?P<acc>[^|]+)\|(?P<entry>\S+) (?P<rest>.*)$")
FIELD = re.compile(r" (OS|OX|GN|PE|SV)=")


def parse(rest: str) -> dict[str, str]:
    """Split 'Protein name OS=... OX=... GN=... PE=... SV=...' into its fields."""
    parts = FIELD.split(rest)
    out = {"protein_name": parts[0].strip()}
    for key, value in zip(parts[1::2], parts[2::2]):
        out[key] = value.strip()
    return out


def main(in_path: str, out_fasta: str, out_tsv: str) -> None:
    n = 0
    with gzip.open(in_path, "rt") as fin, open(out_fasta, "w") as fa, open(out_tsv, "w") as tsv:
        tsv.write("accession\tentry_name\tprotein_name\torganism\ttaxon_id\tgene_name\n")
        for line in fin:
            if not line.startswith(">"):
                fa.write(line)
                continue
            m = HEADER.match(line.rstrip("\n"))
            if not m:
                sys.exit(f"not a Swiss-Prot header: {line.strip()}")
            f = parse(m["rest"])
            fa.write(f">{m['acc']}\n")
            tsv.write("\t".join([m["acc"], m["entry"], f["protein_name"], f.get("OS", ""),
                                 f.get("OX", ""), f.get("GN", "")]) + "\n")
            n += 1
    print(f"{n} Swiss-Prot entries")


if __name__ == "__main__":
    main(*sys.argv[1:4])
