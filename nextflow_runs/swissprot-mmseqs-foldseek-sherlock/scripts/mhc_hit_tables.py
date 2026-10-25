"""Build the MHC hit tables from the pre-filtered MMseqs2 and Foldseek results.

Inputs (from `make pull`):
  prefiltered/<run>.<tool>.mhc_c1set_hits.tsv  hits to Swiss-Prot entries with an
                                              MHC or C1-set Pfam domain
  prefiltered/<run>.<tool>.best_hit.tsv        each Botryllus protein's best hit
  prefiltered/<run>.<tool>.n_hits.tsv          number of hits per Botryllus protein
  prefiltered/mhc_exhaustive.<tool>.mhc_c1set_hits.tsv
                                              every-pair search against the 384
                                              MHC/C1-set entries only
  prefiltered/swissprot_headers.tsv.gz         accession -> Swiss-Prot FASTA header
  reference/swissprot_mhc_c1set_pfam.tsv       UniProt query result (see README.md)

Outputs in --outdir:
  mhc_hits.tsv                 every hit to an MHC class I or II protein
  mhc_best_hit_per_protein.tsv one row per (tool, Botryllus protein): its best MHC hit
  b2m_and_c1set_summary.tsv    number of Botryllus proteins hitting each target class, per tool
  b2m_hits.tsv                 every hit to B2M, the MHC class I light chain
  exhaustive_check_summary.tsv how the exhaustive search against the 384
                               MHC/C1-set entries compares with the full search
  best_hit_per_protein_and_class.tsv  each protein's best hit in each target class

groove_start..groove_end is the target's peptide-binding region: from the end
of the signal peptide to the start of the first Ig-like domain (UniProt
features). aligned_aa_in_peptide_binding_region counts the aligned target
residues inside it, gaps included; 0 means the hit lies only in the Ig-like
part. aligned_pairs_in_peptide_binding_region counts only target residues
aligned to a Botryllus residue (no gap), from the exhaustive search's aligned
strings (qaln, taln). Both are left
empty for AlphaFold models of isoforms (target_is_canonical false), whose
residue numbers do not match the UniProt coordinates.
"""

import argparse
from pathlib import Path

import polars as pl

TOOLS = {
    "mmseqs_swissprot": "MMseqs2",
    "foldseek_afdb_swissprot": "Foldseek (ProstT5 3Di)",
}

# Pfam families that mark an MHC class I or class I-like heavy chain, and the
# class II alpha and beta chains.
CLASS_I_PFAM = {
    "PF00129": "MHC_I",
    "PF06623": "MHC_I_C",
    "PF16497": "MHC_I_3",  # CD1, EPCR
}
CLASS_II_PFAM = {"PF00993": "MHC_II_alpha", "PF00969": "MHC_II_beta"}
C1_SET = "PF07654"

# Named in the task, in case an entry lacks the Pfam annotation.
MHC_GENES = {
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-H",
    "CD1A", "CD1B", "CD1C", "CD1D", "CD1E", "MR1", "FCGRT", "AZGP1",
    "MICA", "MICB", "HFE", "PROCR",
}
VERTEBRATA_TAXID = "7742"

NUMERIC_COLUMNS = [
    "pident", "fident", "alnlen", "qstart", "qend", "qlen", "tstart", "tend",
    "tlen", "evalue", "bits", "prob", "qcov", "tcov",
]


def hits_columns(stem: Path) -> list[str]:
    with open(f"{stem}.mhc_c1set_hits.tsv") as f:
        return f.readline().rstrip("\n").split("\t")


def classify(pfam: str, genes: str) -> str:
    pf = set(filter(None, (pfam or "").split(";")))
    gene_set = set((genes or "").split())
    if pf & set(CLASS_II_PFAM):
        return "MHC class II"
    if pf & set(CLASS_I_PFAM) or gene_set & MHC_GENES or any(
        g.startswith("HLA-D") for g in gene_set
    ):
        return "MHC class I or I-like"
    if "B2M" in gene_set:
        return "B2M"
    if C1_SET in pf:
        return "other C1-set"
    return "none"


def load_reference(path: Path) -> pl.DataFrame:
    ref = pl.read_csv(path, separator="\t", infer_schema_length=0).rename(
        {
            "Entry": "target_accession",
            "Entry Name": "target_entry",
            "Gene Names": "target_genes",
            "Protein names": "target_protein",
            "Organism": "target_organism",
            "Organism (ID)": "target_taxid",
            "Pfam": "target_pfam",
            "Taxonomic lineage (Ids)": "lineage",
            "Signal peptide": "signal",
            "Domain [FT]": "domains",
        }
    )
    # The peptide-binding part of an MHC chain (class I alpha1-alpha2, class
    # II alpha1 or beta1) lies between the signal peptide and the first
    # Ig-like domain, which UniProt annotates; the Ig-like domain itself is
    # shared with thousands of non-MHC proteins.
    ref = ref.with_columns(
        (
            pl.col("signal").str.extract(r"SIGNAL \d+\.\.(\d+)", 1).cast(pl.Int64).fill_null(0)
            + 1
        ).alias("groove_start"),
        (
            pl.col("domains").str.extract(r"DOMAIN (\d+)\.\.\d+; /note=\"Ig-like", 1).cast(pl.Int64)
            - 1
        ).alias("groove_end"),
    )
    return ref.with_columns(
        pl.struct("target_pfam", "target_genes")
        .map_elements(
            lambda r: classify(r["target_pfam"], r["target_genes"]),
            return_dtype=pl.String,
        )
        .alias("target_class"),
        pl.col("lineage")
        .str.contains(rf"(^|, ){VERTEBRATA_TAXID} \(")
        .alias("target_is_vertebrate"),
    ).drop("lineage", "signal", "domains")


def load_headers(path: Path) -> pl.DataFrame:
    return pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=["accession", "header"],
        quote_char=None,
    )


def accession_expr(col: str) -> pl.Expr:
    """sp|P01892|HLA_A2_HUMAN -> P01892; AF-P01892-F1-model_v4 -> P01892."""
    c = pl.col(col)
    return (
        pl.when(c.str.starts_with("AF-"))
        .then(c.str.split("-").list.get(1, null_on_oob=True))
        .when(c.str.contains(r"\|"))
        .then(c.str.split("|").list.get(1, null_on_oob=True))
        .otherwise(c)
    )


def read_hits(path: Path) -> pl.DataFrame:
    cols = hits_columns(Path(str(path).removesuffix(".mhc_c1set_hits.tsv")))
    return pl.read_csv(
        path, separator="\t", quote_char=None, infer_schema_length=10000,
    ).with_columns(
        pl.col(c).cast(pl.Float64, strict=False) for c in NUMERIC_COLUMNS if c in cols
    ).with_columns(pl.col("hit_rank").cast(pl.Int64))


def pairs_in_region(taln: str, qaln: str, tstart: int, lo: int, hi: int) -> tuple[int, int]:
    """Count aligned residue pairs (no gap on either side) whose target residue
    lies in lo..hi, and how many of those are identical."""
    pos, n, same = tstart, 0, 0
    for t, q in zip(taln, qaln):
        if t != "-":
            if q != "-" and lo <= pos <= hi:
                n += 1
                same += t == q
            pos += 1
    return n, same


def rank_by_bits(hits: pl.DataFrame, bits_path: Path) -> pl.DataFrame:
    """1 + the number of the protein's full-search hits with a higher bit
    score. Works for hits only the exhaustive check found, which have no
    position in the full search's list."""
    full = pl.read_csv(bits_path, separator="\t", has_header=False,
                       new_columns=["query", "full_bits"])
    ranks = (
        hits.select("query", "bits").unique()
        .join(full, on="query", how="left")
        .group_by("query", "bits")
        .agg((pl.col("full_bits") > pl.col("bits")).sum().alias("n_higher"))
        .with_columns((pl.col("n_higher") + 1).cast(pl.Int64).alias("rank_by_bits"))
        .drop("n_higher")
    )
    return hits.join(ranks, on=["query", "bits"], how="left")


def add_exhaustive_hits(
    hits: pl.DataFrame, exh: pl.DataFrame, label: str
) -> tuple[pl.DataFrame, dict]:
    """Add hits that only the exhaustive search against the 384 MHC/C1-set
    entries found (the full search keeps 2,000 candidates per protein).

    The exhaustive search ran against a smaller database, so its e-values are
    smaller. They are put on the full-database scale by the median ratio
    full/exhaustive over the pairs both searches found."""
    pair = ["query", "target"]
    shared = exh.join(
        hits.select(*pair, pl.col("evalue").alias("evalue_full")), on=pair
    ).with_columns((pl.col("evalue_full") / pl.col("evalue")).alias("ratio"))
    factor = shared["ratio"].median()
    new = (
        exh.join(hits, on=pair, how="anti")
        .with_columns(
            (pl.col("evalue") * factor).alias("evalue"),
            pl.lit(None, dtype=pl.Int64).alias("hit_rank"),
        )
        .filter(pl.col("evalue") < 1e-3)
    )
    summary = {
        "tool": label,
        "pairs_full_search": hits.height,
        "pairs_exhaustive": exh.height,
        "pairs_in_both": shared.height,
        "pairs_full_only": hits.join(exh, on=pair, how="anti").height,
        "evalue_scale_factor_median": factor,
        "evalue_scale_factor_min": shared["ratio"].min(),
        "evalue_scale_factor_max": shared["ratio"].max(),
        "pairs_exhaustive_only_e_below_1e-3": new.height,
        "proteins_exhaustive_only_e_below_1e-3": new["query"].n_unique(),
    }
    # The exhaustive search wrote the aligned residues; carry them over to the
    # pairs the full search also found.
    if "qaln" in exh.columns:
        hits = hits.join(exh.select(*pair, "qaln", "taln"), on=pair, how="left")
    merged = pl.concat(
        [
            hits.with_columns(pl.lit("full search").alias("found_by")),
            new.with_columns(pl.lit("exhaustive check only").alias("found_by")),
        ],
        how="diagonal_relaxed",
    )
    return merged, summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prefiltered", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--run", default="full")
    ap.add_argument("--exhaustive-run", default="mhc_exhaustive")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    ref = load_reference(args.reference)
    headers = load_headers(args.prefiltered / "swissprot_headers.tsv.gz")

    hit_tables, summaries, checks = [], [], []
    for tool, label in TOOLS.items():
        stem = args.prefiltered / f"{args.run}.{tool}"
        if not Path(f"{stem}.mhc_c1set_hits.tsv").exists():
            print(f"skipping {label}: no {stem}.mhc_c1set_hits.tsv")
            continue
        hits = read_hits(Path(f"{stem}.mhc_c1set_hits.tsv"))
        exh_path = args.prefiltered / f"{args.exhaustive_run}.{tool}.mhc_c1set_hits.tsv"
        if exh_path.exists():
            hits, check = add_exhaustive_hits(hits, read_hits(exh_path), label)
            checks.append(check)
            hits = rank_by_bits(hits, Path(f"{stem}.bits_for_ranking.tsv"))
        else:
            print(f"no exhaustive check for {label}: {exh_path} missing")
        best = (
            pl.read_csv(f"{stem}.best_hit.tsv", separator="\t", quote_char=None,
                        infer_schema_length=10000)
            .with_columns(accession_expr("target").alias("best_hit_accession"))
            .join(headers, left_on="best_hit_accession", right_on="accession", how="left")
            .select(
                "query",
                "best_hit_accession",
                pl.col("header").alias("best_hit_header"),
                pl.col("evalue").alias("best_hit_evalue"),
                pl.col("bits").alias("best_hit_bits"),
            )
        )
        n_hits = pl.read_csv(
            f"{stem}.n_hits.tsv", separator="\t", has_header=False,
            new_columns=["query", "n_hits"],
        )
        hits = (
            hits.join(ref, on="target_accession", how="left")
            .join(n_hits, on="query", how="left")
            .join(best, on="query", how="left")
            .with_columns(
                pl.lit(label).alias("tool"),
                (pl.col("rank_by_bits") == 1).alias("mhc_is_best_hit")
                if "rank_by_bits" in hits.columns
                else (pl.col("hit_rank") == 1).alias("mhc_is_best_hit"),
                # The full search keeps at most 2,000 hits per protein, so for
                # a protein at the limit the rank is a lower bound.
                (pl.col("n_hits") >= 2000).alias("protein_at_2000_hit_limit"),
                pl.col("target").alias("target_id"),
                # AlphaFold DB also holds models of isoforms (AF-Q30201-6-F1),
                # whose residue numbers differ from the canonical sequence
                # that the UniProt domain coordinates refer to.
                (~pl.col("target").str.contains(r"^AF-[^-]+-\d+-F")).alias(
                    "target_is_canonical"
                ),
            ).with_columns(
                pl.when(pl.col("target_is_canonical"))
                .then(
                    (
                        # No Ig-like domain (ULBP, RAE-1, EPCR): the whole chain.
                        pl.min_horizontal(
                            "tend", pl.col("groove_end").fill_null(pl.col("tlen"))
                        )
                        - pl.max_horizontal("tstart", "groove_start")
                        + 1
                    ).clip(lower_bound=0)
                )
                .alias("aligned_aa_in_peptide_binding_region"),
            )
        )
        if "taln" in hits.columns:
            counts = [
                pairs_in_region(r["taln"], r["qaln"], int(r["tstart"]),
                                r["groove_start"], r["groove_end"] or int(r["tlen"]))
                if r["taln"] is not None and r["target_is_canonical"]
                and r["groove_start"] is not None else (None, None)
                for r in hits.select(
                    "taln", "qaln", "tstart", "tlen", "groove_start", "groove_end",
                    "target_is_canonical",
                ).iter_rows(named=True)
            ]
            hits = hits.with_columns(
                pl.Series("aligned_pairs_in_peptide_binding_region",
                          [c[0] for c in counts], dtype=pl.Int64),
                pl.Series("identical_pairs_in_peptide_binding_region",
                          [c[1] for c in counts], dtype=pl.Int64),
            )
        hit_tables.append(hits)
        summaries.append(
            hits.group_by("tool", "target_class").agg(
                pl.col("query").n_unique().alias("n_botryllus_proteins"),
                pl.len().alias("n_hits"),
                pl.col("evalue").min().alias("best_evalue"),
                pl.col("bits").max().alias("best_bits"),
            )
        )

    all_hits = pl.concat(hit_tables, how="diagonal_relaxed")
    keep = [
        "tool", "query", "target_id", "target_is_canonical", "target_accession", "target_entry", "target_genes",
        "target_protein", "target_organism", "target_is_vertebrate",
        "target_class", "target_pfam", "evalue", "bits", "pident", "fident",
        "prob", "alnlen", "qstart", "qend", "qlen", "tstart", "tend", "tlen",
        "qcov", "tcov", "groove_start", "groove_end",
        "aligned_aa_in_peptide_binding_region",
        "aligned_pairs_in_peptide_binding_region",
        "identical_pairs_in_peptide_binding_region",
        "hit_rank", "rank_by_bits", "n_hits", "protein_at_2000_hit_limit",
        "mhc_is_best_hit", "found_by",
        "best_hit_accession", "best_hit_header", "best_hit_evalue", "best_hit_bits",
        "qaln", "taln",
    ]
    keep = [c for c in keep if c in all_hits.columns]
    mhc = (
        all_hits.filter(pl.col("target_class").str.starts_with("MHC"))
        .select(keep)
        .sort("tool", "evalue")
    )
    mhc.write_csv(args.outdir / "mhc_hits.tsv", separator="\t")

    mhc.sort("evalue").group_by("tool", "query", maintain_order=True).first().write_csv(
        args.outdir / "mhc_best_hit_per_protein.tsv", separator="\t"
    )

    # Each Botryllus protein's best hit in each target class (MHC class I, II,
    # B2M, other C1-set), for the figure.
    all_hits.select(keep).sort("evalue").group_by(
        "tool", "query", "target_class", maintain_order=True
    ).first().sort("tool", "target_class", "evalue").write_csv(
        args.outdir / "best_hit_per_protein_and_class.tsv", separator="\t"
    )

    pl.concat(summaries).sort("tool", "target_class").write_csv(
        args.outdir / "b2m_and_c1set_summary.tsv", separator="\t"
    )
    # B2M hits in full, since B2M is the MHC class I light chain.
    all_hits.filter(pl.col("target_class") == "B2M").select(keep).sort(
        "tool", "evalue"
    ).write_csv(args.outdir / "b2m_hits.tsv", separator="\t")

    if checks:
        pl.DataFrame(checks).write_csv(
            args.outdir / "exhaustive_check_summary.tsv", separator="\t"
        )
        print(pl.DataFrame(checks).transpose(include_header=True))
    print(pl.concat(summaries).sort("tool", "target_class"))
    print(mhc.head(20))


if __name__ == "__main__":
    main()
