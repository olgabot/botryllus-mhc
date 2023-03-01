import glob
import os

import pandas as pd


def read_csvs(folder, csv_glob='*.csv'):
    dfs = []
    csvs = glob.glob(os.path.join(folder, csv_glob))
    for csv in csvs:
        df = pd.read_csv(csv)
        dfs.append(df)
    concatenated = pd.concat(dfs, ignore_index=True)
    return concatenated




def extract_pfam_metadata(df):
    pfam_metadata_split1 = df.name.str.rstrip(";").str.split(";", expand=True)
    pfam_metadata_split2 = pfam_metadata_split1[0].str.split(expand=True)

    pfam_name = pfam_metadata_split1[1]
    pfam_name.name = "pfam_name"
    pfam_metadata = pd.concat([pfam_name, pfam_metadata_split2], axis=1)
    pfam_metadata = pfam_metadata.rename(columns={2: "pfam_id_full"})
    pfam_metadata["pfam_id"] = pfam_metadata.pfam_id_full.str.split(".").str[0]
    return pfam_metadata


def reorder_cols(df, first_cols=["similarity", "query_name"]):
    columns_reordered = first_cols + df.columns.difference(first_cols).tolist()
    df = df[columns_reordered]
    return df


def add_pfam_metadata(df, first_cols=["pfam_name", "similarity", "query_name"]):
    metadata = extract_pfam_metadata(df)
    df_with_metadata = pd.concat([df, metadata], axis=1)
    df_with_metadata = reorder_cols(
        df_with_metadata, first_cols
    )
    return df_with_metadata

def extract_gencode_gene_symbol(df, symbol_regex="gene_symbol:([\w\d\-]+)"):
    df["symbol"] = df["name"].str.extract(symbol_regex)
    return df

def read_gencode_folder(folder, csv_glob, first_cols, symbol_regex):
    df = read_csvs(folder, csv_glob)
    df = extract_gencode_gene_symbol(df, symbol_regex)
    df = reorder_cols(df, first_cols)
    return df
