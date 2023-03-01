import os
import re

import pandas as pd
from kmer_utils import get_encoded_kmer_hashvals
from sourmash.logging import notify, error, print_results, set_quiet


def sanitize(x):
    """Clean a gene name so it is a nice filename"""
    return x.replace(" ", "_").replace(".", "-dot-").replace(':', '-colon-')


GENE_SYMBOL_PATTERN = re.compile("gene_name=([\w\d]+)")

SYMBOL_SEPARATOR = "---"


def annotate_found_sig(
        found_sig, found_seq, query_sig, query_kmer_hashvals, scaled,
        gene_symbol_pattern=GENE_SYMBOL_PATTERN
):
    """Get the containment and overlapping k-mers for each found signature"""
    containment = found_sig.contained_by(query_sig)
    symbol = re.findall(gene_symbol_pattern, found_sig.name)[0]

    found_kmer_hashvals = get_encoded_kmer_hashvals(
        found_seq, found_sig.name, sigobj=query_sig
    )
    contained_kmer_hashvals = query_kmer_hashvals.merge(
        found_kmer_hashvals, suffixes=("_query", "_found"), on=("hashval", "kmer_hp")
    )
    contained_kmer_hashvals["n_kmers"] = len(contained_kmer_hashvals)
    contained_kmer_hashvals["intersect_bp"] = (
            scaled * contained_kmer_hashvals["n_kmers"]
    )
    contained_kmer_hashvals["containment"] = containment
    contained_kmer_hashvals["symbol"] = symbol
    return contained_kmer_hashvals


def add_homology_information_to_matches(
        query_kmer_matches,
        human_mouse_homologs,
        class_col="DB Class Key",
        homolog_group_col="homolog_group",
        coord_col="Genomic Coordinates (mouse: , human: )",
        symbol_separator=SYMBOL_SEPARATOR,
):
    query_kmer_matches[homolog_group_col] = None

    # Annotate any genes that are human-mouse homologs
    for symbol, df in query_kmer_matches.groupby("symbol"):
        found_homolog_subset = human_mouse_homologs.query("Symbol == @symbol")

        if found_homolog_subset.empty:
            # No matches found in mouse-human homologs, continue on to next one
            continue

        found_homolog_rows = human_mouse_homologs[class_col].isin(
            found_homolog_subset[class_col]
        )
        found_homolog_groups = human_mouse_homologs.loc[found_homolog_rows]
        homolog_group_name = symbol_separator.join(
            sorted(found_homolog_groups["Symbol"])
        )

        all_homologs_found = (
            found_homolog_groups["Symbol"].isin(query_kmer_matches.symbol).all()
        )

        # Assign values
        query_kmer_matches.loc[df.index, "genomic_coord"] = found_homolog_subset[
            coord_col
        ].values[0]
        query_kmer_matches.loc[df.index, homolog_group_col] = homolog_group_name
        query_kmer_matches.loc[df.index, "all_homologs_found"] = all_homologs_found
    return query_kmer_matches


def single_gather_multi_db(query_sig, query_seq, dbs, sequences, human_mouse_homologs,
                           outdir, threshold_bp=10, scaled=5, force=False):
    csv = os.path.join(outdir, f"{sanitize(query_sig.name)}.csv")

    # Early exit if this query was already done
    if os.path.exists(csv) and not force:
        return

    # Early exit if this query is empty
    if len(query_sig.minhash.hashes) == 0:
        return

    dfs = []
    query_kmer_hashvals = get_encoded_kmer_hashvals(query_seq, query_sig.name,
                                                    sigobj=query_sig)

    for species, db in dbs.items():
        species_seqs = sequences[species]
        try:
            counter = db.counter_gather(query_sig, threshold_bp=threshold_bp)
        except ValueError:
            # Expected: ValueError: requested threshold_bp is unattainable with this query 
            # The query is too small
            return

        # Get hash and k-mer information for each matching signature
        for i, found_sig in enumerate(counter.siglist):
            found_seq = species_seqs[found_sig.name]
            df = annotate_found_sig(found_sig, found_seq, query_sig,
                                    query_kmer_hashvals, scaled)
            df["species"] = species
            df["found_i"] = i

            dfs.append(df)

    if dfs:
        query_kmer_matches = pd.concat(dfs, ignore_index=True)
        query_kmer_matches = add_homology_information_to_matches(
            query_kmer_matches, human_mouse_homologs
        )

        query_kmer_matches.to_csv(csv, index=False)
    else:
        notify(f"{query_sig.name} had no matches")
