
import sourmash
from orpheum.sequence_encodings import encode_peptide


import pandas as pd

def get_encoded_kmer_hashvals(sequence, name, encoding="hp", K=24, sigobj=None):
    lines = []
    for i in range(0, len(sequence) - K + 1):
        kmer = sequence[i : i + K]
        kmer_encoded = encode_peptide(kmer, encoding)
        hashvals = sigobj.minhash.seq_to_hashes(kmer, is_protein=True)

        # if len(hashval) == 1:
        #     hashval = hashval[0]
        # else:
        #     print(f"More than one hashval found for {kmer}")
        for hashval in hashvals:
            line = [i, kmer, kmer_encoded, hashval]
            lines.append(line)
    kmer_to_hashes = pd.DataFrame(
        lines, columns=["i", "kmer", f"kmer_{encoding}", "hashval"]
    )
    kmer_to_hashes["name"] = name
    return kmer_to_hashes


def get_matching_kmer_subsequence(sequence, kmers, gene_symbol):
    i_min = len(sequence)
    i_max = 0
    for kmer_human in kmers:
        i_kmer = sequence.find(kmer_human)
        # if i_kmer =< 0:
        #     continue

        if i_kmer < i_min:
            i_min = i_kmer
        if i_kmer > i_max:
            i_max = i_kmer
    j = i_max + len(kmer_human)

    #     print(f'i_min: {i_min}')
    #     print(f'i_max: {i_max}')
    #     print(f'j: {j}')

    
    print(f">{gene_symbol}")
    print(f"before match: {sequence[:i_min]}\n")
    print(f"matching: {sequence[i_min:j]}\n")
    print(f"after match: {sequence[j:]}")
    print(f"Match range (1-based): {i_min+1}-{j}")
    return sequence[i_min:j]


def subset_gene_kmers(merged_kmers, col, gene_symbol, other='human'):
    gene_subset = merged_kmers.query(f"{col} == @gene_symbol")
    
    kmer_other = f"kmer__{other}"

    tidy = pd.concat(
        [
            gene_subset[[kmer_other, "kmer_hp"]],
            gene_subset[["kmer__botryllus", "kmer_hp"]],
        ]
    )
    tidy = tidy.drop_duplicates()
    tidy["species"] = tidy[kmer_other].map(
        lambda x: "botryllus" if pd.isnull(x) else other
    )
    tidy["kmer_seq"] = tidy.apply(
        lambda x: x[kmer_other]
        if pd.isnull(x["kmer__botryllus"])
        else x["kmer__botryllus"],
        axis=1,
    )
    pivoted = tidy.pivot(columns="kmer_hp", index="species", values="kmer_seq")
    return pivoted
