#!/usr/bin/env python3
"""

Summarise the kmerseek Botryllus vs Swiss-Prot search across alphabet-ksize pairs.

Reads one <chunk>.<alphabet>.k<k>.best.tsv per search task (one row per Botryllus protein
with a kept match: number of Swiss-Prot proteins matched, number at the strong cutoff or
below, and the best match) and writes three files:

<prefix>.pairs.tsv
    one row per alphabet-ksize pair: whether the E-value fit worked, which number the
    matches were kept on (region_evalue, or region_tail_probability without a fit), and
    how many of the Botryllus proteins have a kept match and a strong match.
<prefix>.best_per_protein_per_pair.tsv.gz
    every Botryllus protein x pair with a kept match, with its best Swiss-Prot match.
<prefix>.per_protein.tsv.gz
    one row per Botryllus protein (all of them, matched or not): in how many of the fitted
    pairs it has a strong match, and its best match over the fitted pairs, with the
    Swiss-Prot protein name and organism.

Standard library only (Sherlock's python module has no polars).

Usage: collect_summary.py --pairs pairs.tsv --best-dir best --fit-status fit_status.tsv
           --metadata swissprot.metadata.tsv --queries chunk_*.fasta --strong 1e-3
           --out-prefix kmerseek_botryllus_swissprot
"""

import argparse
import csv
import gzip
from pathlib import Path


def query_names(paths: list[str]) -> list[str]:
    names = []
    for p in sorted(paths):
        with open(p) as fh:
            names += [line[1:].strip() for line in fh if line.startswith(">")]
    return names


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--best-dir", required=True)
    ap.add_argument("--fit-status", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--queries", nargs="+", required=True)
    ap.add_argument("--strong", type=float, required=True)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    pairs = list(csv.DictReader(open(a.pairs), delimiter="\t"))
    fit = {}
    for line in open(a.fit_status):
        if line.strip():
            alphabet, k, status = line.rstrip("\n").split("\t")
            fit[(alphabet, k)] = status
    meta = {r["accession"]: r for r in csv.DictReader(open(a.metadata), delimiter="\t")}
    queries = query_names(a.queries)
    n_queries = len(queries)

    # (alphabet, k) -> query -> row
    best: dict[tuple[str, str], dict[str, dict]] = {(p["alphabet"], p["ksize"]): {} for p in pairs}
    for f in sorted(Path(a.best_dir).glob("*.best.tsv")):
        # chunk_0000.<alphabet>.k<k>.best.tsv
        _, alphabet, kpart, _, _ = f.name.split(".")
        rows = best[(alphabet, kpart[1:])]
        for r in csv.DictReader(open(f), delimiter="\t"):
            rows[r["query_name"]] = r

    with open(f"{a.out_prefix}.pairs.tsv", "w") as out:
        out.write("alphabet\tksize\tprotein20_k\tbits_per_kmer\tfit\tstatistic\t"
                  f"n_botryllus_proteins\tn_with_match\tn_with_strong_match\tshare_with_strong_match\treference_alphabet\treference_k\treference_bits\trounding\n")
        for p in pairs:
            key = (p["alphabet"], p["ksize"])
            rows = best[key]
            stat = {"fitted": "region_evalue", "unfitted": "region_tail_probability"}.get(fit.get(key), "NA")
            n_strong = sum(1 for r in rows.values() if int(r["n_targets_strong"]) > 0)
            out.write(f"{p['alphabet']}\t{p['ksize']}\t{p['protein20_k']}\t{p['bits_per_kmer']}\t"
                      f"{fit.get(key, 'missing')}\t{stat}\t{n_queries}\t{len(rows)}\t{n_strong}\t"
                      f"{n_strong / n_queries:.4f}\t{p.get('reference_alphabet', 'protein20')}\t"
                      f"{p.get('reference_k', p['protein20_k'])}\t{p.get('reference_bits', '')}\t"
                      f"{p.get('rounding', 'ceil')}\n")

    with gzip.open(f"{a.out_prefix}.best_per_protein_per_pair.tsv.gz", "wt") as out:
        out.write("query_name\talphabet\tksize\tprotein20_k\tstatistic\tn_targets\t"
                  "n_targets_strong\tbest_target\tbest_value\n")
        for p in pairs:
            for q, r in sorted(best[(p["alphabet"], p["ksize"])].items()):
                out.write(f"{q}\t{p['alphabet']}\t{p['ksize']}\t{p['protein20_k']}\t{r['statistic']}\t"
                          f"{r['n_targets']}\t{r['n_targets_strong']}\t{r['best_target']}\t{r['best_value']}\n")

    fitted_pairs = [p for p in pairs if fit.get((p["alphabet"], p["ksize"])) == "fitted"]
    with gzip.open(f"{a.out_prefix}.per_protein.tsv.gz", "wt") as out:
        out.write(f"query_name\tn_fitted_pairs\tn_pairs_with_strong_match\tbest_alphabet\tbest_ksize\t"
                  "best_evalue\tbest_target\tbest_target_entry\tbest_target_protein_name\t"
                  "best_target_organism\n")
        for q in queries:
            n_strong, top = 0, None
            for p in fitted_pairs:
                r = best[(p["alphabet"], p["ksize"])].get(q)
                if r is None:
                    continue
                n_strong += int(r["n_targets_strong"]) > 0
                v = float(r["best_value"])
                if top is None or v < top[0]:
                    top = (v, p["alphabet"], p["ksize"], r["best_target"])
            if top is None:
                out.write(f"{q}\t{len(fitted_pairs)}\t0\t\t\t\t\t\t\t\n")
                continue
            m = meta.get(top[3], {})
            out.write(f"{q}\t{len(fitted_pairs)}\t{n_strong}\t{top[1]}\t{top[2]}\t{top[0]:g}\t{top[3]}\t"
                      f"{m.get('entry_name', '')}\t{m.get('protein_name', '')}\t{m.get('organism', '')}\n")
    print(f"{n_queries} Botryllus proteins, {len(pairs)} pairs, {len(fitted_pairs)} with an E-value fit")


if __name__ == "__main__":
    main()
