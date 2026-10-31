"""Count UniFIRE's whole-proteome predictions: proteins with any prediction per rule set,
and how many predictions of each annotation type.

Usage:
    python3 summarise_unifire.py <unifire_whole_proteome dir> <n proteins sent to UniFIRE>

Writes unifire_summary.tsv and unifire_annotation_types.tsv in that folder, and the figure
unifire_whole_proteome.{pdf,svg,png} next to them.
"""

import sys
from pathlib import Path

import polars as pl

import pubfig as pf

RULE_SETS = {"arba": "ARBA", "unirule": "UniRule", "unirule-pirsr": "PIRSR"}


def read(folder, key):
    return pl.read_csv(
        folder / f"predictions_{key}.out.gz", separator="\t", quote_char=None,
        infer_schema=False,
    ).with_columns(rule_set=pl.lit(RULE_SETS[key]))


def main():
    folder, n_total = Path(sys.argv[1]), int(sys.argv[2])
    preds = pl.concat([read(folder, k) for k in RULE_SETS])
    per_set = (
        preds.group_by("rule_set")
        .agg(proteins=pl.col("ProteinId").n_unique(), predictions=pl.len())
        .vstack(pl.DataFrame({
            "rule_set": ["any"], "proteins": [preds["ProteinId"].n_unique()],
            "predictions": [preds.height]}, schema_overrides={"proteins": pl.UInt32, "predictions": pl.UInt32}))
        .with_columns(proteins_sent=pl.lit(n_total),
                      percent=(100 * pl.col("proteins") / n_total).round(1))
        .sort("proteins", descending=True)
    )
    print(per_set)
    per_set.write_csv(folder / "unifire_summary.tsv", separator="\t")
    types = (
        preds.with_columns(kind=pl.col("AnnotationType").str.split(".").list.first())
        .group_by("rule_set", "kind").agg(predictions=pl.len())
        .sort("rule_set", "predictions", descending=[False, True])
    )
    print(types)
    types.write_csv(folder / "unifire_annotation_types.tsv", separator="\t")

    pf.use_style()
    fig, (a, b) = pf.figure(pf.TWO_COLUMN_MM, 55, ncols=2, width_ratios=[1, 1.6])
    order = ["any", "ARBA", "UniRule", "PIRSR"]
    rows = {r["rule_set"]: r for r in per_set.iter_rows(named=True)}
    colours = {"any": pf.GREY, "ARBA": pf.OKABE_ITO["blue"], "UniRule": pf.OKABE_ITO["vermillion"], "PIRSR": pf.OKABE_ITO["bluish_green"]}
    y = range(len(order))
    a.barh(list(y), [rows[k]["proteins"] for k in order], color=[colours[k] for k in order])
    for i, k in enumerate(order):
        a.text(rows[k]["proteins"], i, f" {rows[k]['proteins']:,} ({rows[k]['percent']}%)", va="center")
    a.set_yticks(list(y), order)
    a.invert_yaxis()
    a.set_xlim(0, n_total)
    a.set_xlabel(f"Proteins with a prediction (of {n_total:,})")

    kinds = types.group_by("kind").agg(pl.col("predictions").sum()).sort("predictions", descending=True)["kind"].to_list()
    left = [0] * len(kinds)
    for rs in ["ARBA", "UniRule", "PIRSR"]:
        sub = dict(types.filter(pl.col("rule_set") == rs).select("kind", "predictions").iter_rows())
        vals = [sub.get(k, 0) for k in kinds]
        b.barh(range(len(kinds)), vals, left=left, color=colours[rs], label=rs)
        left = [l + v for l, v in zip(left, vals)]
    names = {"keyword": "UniProt keyword", "xref": "GO term", "comment": "function, location, activity",
             "feature": "site or region", "protein": "protein name", "gene": "gene name"}
    b.set_yticks(range(len(kinds)), [names.get(k, k) for k in kinds])
    b.invert_yaxis()
    b.set_xlabel("Predictions")
    b.legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3)
    pf.panel_label(a, "a")
    pf.panel_label(b, "b")
    pf.save(fig, str(folder / "unifire_whole_proteome"))


if __name__ == "__main__":
    main()
