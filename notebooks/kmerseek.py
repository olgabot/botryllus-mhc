import csv
import os

import click
import sourmash_plugin_branchwater
import sourmash
import pandas as pd

from sig2kmer import get_matching_hashes_in_file

NOTIFY_EVERY_BP = int(1e5)


protein_moltypes = "protein", "dayhoff", "hp"


def sketch(fasta, moltype, ksize, scale):
    param_string = f"{moltype},k={ksize},scale={scale}"
    sigfile = f"{fasta}.{moltype}.k{ksize}.scale{scale}.sig.zip"
    sourmash_plugin_branchwater.do_manysketch(
        fasta,
        param_string,
        output=sigfile,
        singleton=True,
        force=True,
    )
    return sigfile


def _setup_kmer_writer(sig):
    output_kmers = f"{sig}.kmers.csv"

    kmerout_fp = open(output_kmers, "wt")
    kmerout_w = csv.writer(kmerout_fp)
    kmerout_w.writerow(
        [
            "kmer_in_sequence",
            "kmer_in_alphabet",
            "hashval",
            "start",
            "read_name",
            "filename",
        ]
    )
    return output_kmers, kmerout_fp, kmerout_w


def get_kmers(sig, fasta, moltype, ksize, scale):
    # Scale is not used, but is accepted for ease of keyword argument sharing

    output_kmers, kmerout_fp, kmerout_w = _setup_kmer_writer(sig)

    # first, load the signature and extract the hashvals
    sigobj = list(
        sourmash.load_file_as_signatures(sig, ksize=ksize, select_moltype=moltype)
    ).pop()
    query_hashvals = set(sigobj.minhash.hashes.keys())
    query_ksize = sigobj.minhash.ksize

    # TODO: Simplify and remove all the tracking numbers, e.g. number residues loaded, kmers found, etc.
    # I really don't like how this relies on so many parameters for tracking how far we've
    # progressed through the files. I wish this function was massively simpler, not accepting any of the
    # keyword arguments initialized with a 0
    get_matching_hashes_in_file(
        fasta,
        query_ksize,
        moltype,
        input_is_protein=True,
        hashes=query_hashvals,
        found_kmers=0,
        m=0,  # bp in found sequences
        n=0,  # bp loaded
        n_seq=0,
        seqout_fp=None,
        kmerout_w=kmerout_w,
        watermark=int(NOTIFY_EVERY_BP),
    )

    kmerout_fp.close()
    return output_kmers


@click.command()
def main(
    query_fasta,
    target_sig,
    target_kmers,
    ksize=24,
    moltype="hp",
    scale=5,
    search_output="search.csv",
):
    sketch_kwargs = dict(ksize=ksize, moltype=moltype, scale=scale)

    # Set defaults for search
    search_params = (
        dict(sketch_kwargs)
        .update("threshold", 0)
        .update("prob_significant_overlap", True)
        .update("output_path", search_output)
    )
    # Can this be async?
    query_sig = sketch(query_fasta, **sketch_kwargs)
    query_kmers = get_kmers(query_sig, query_fasta, **sketch_kwargs)

    sourmash_plugin_branchwater.do_multisearch(query_sig, target_sig, **search_params)

    results = pd.read_csv(search_output)

    results_with_kmers = results.join(query_kmers).join(target_kmers)

    results_per_gene = results_with_kmers.groupby("gene")

    for gene, match in results_per_gene:
        domains = get_domains_in_gene(gene, match.target_start, match.target_end)

        print(
            f"Found {match.query_name}:{match.query_start}-{match.query_end} "
            f"in {match.target_name}:{match.target_start}-{match.target_end}"
        )
        if domains:
            for domain in domains:
                print(f"Found: {domain.name} in {domain.start}-{domain.end}")
        print(match.query_seq)
        print(match.hp_kmer)
        print(match.target_seq)
