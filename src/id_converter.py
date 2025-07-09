"""
This script processes gene orthology data to find Uloborus diversus genes
and their Drosophila melanogaster orthologs, and merges this information with
a gene ID converter to provide additional gene descriptions."""

import argparse
import random
import os
import pandas as pd

scripts = os.path.dirname(__file__)
assets = os.path.join(scripts, '..', 'assets')

def make_id_converter():
    """Creates a DataFrame that maps Uloborus diversus transcript IDs to 
    LOCs and nucleotide accession IDs."""

    id_converter_df = pd.read_csv(
        f"{assets}/id_converter.tsv",
        sep="\t",
        dtype=str,
        names=[
            "LOC",
            "txpt",
            "udiv_genes",
            "percent_id",
            "e_value",
            "gene_desc",
            "GO_terms",
        ],
    )

    id_converter_df["desc_index_start"] = (
        id_converter_df["gene_desc"].str.find("[Description: ") + 14
    )
    id_converter_df["desc_index_end"] = id_converter_df["gene_desc"].str.find(
        "] [Gene Type:"
    )

    id_converter_df["Description"] = id_converter_df.apply(
        lambda row: row.gene_desc[row.desc_index_start : row.desc_index_end], axis=1
    )
    id_converter_df = id_converter_df[
        ["LOC", "txpt", "udiv_genes", "GO_terms", "Description"]
    ]

    return id_converter_df

def get_udiv_dmel_genes(NX_tsv, HOG_df, ortholog_tsv, one_random_gene=False):
    """Retrieves Uloborus diversus genes and their Drosophila melanogaster orthologs
    for a given list of HOGs. The orthogroup file for the hierarchical orthogroup node
    of interest must be provided."""

    NX_df = pd.read_csv(NX_tsv, sep="\t", index_col="HOG", dtype=str)

    ortholog_df = pd.read_csv(ortholog_tsv, sep="\t", dtype=str)
    ortholog_df["Uloborus_diversus"] = ortholog_df["Uloborus_diversus"].str.split(", ")
    ortholog_df = ortholog_df.explode("Uloborus_diversus")

    # Attempt to read the HOG DataFrame, if it fails, assume it's already a DataFrame
    # and assign empty columns for udiv_genes and dmel_orthologs
    try:
        df = pd.read_csv(HOG_df, index_col="HOG", dtype=str)
        df = df.assign(udiv_genes="", dmel_orthologs="")
    except:
        df = HOG_df
        df = df.assign(udiv_genes="", dmel_orthologs="")

    # Iterate through each HOG and find corresponding Uloborus diversus genes
    # and their Drosophila melanogaster orthologs
    # If one_random_gene is True, select one random gene per orthogroup
    # Otherwise, list all genes and their orthologs
    for HOG in df.index.to_list():
        try:
            ud_genes = NX_df.at[HOG, "Uloborus_diversus"].split(", ")

            # use if i just want to find 1 U.div gene per OG
            if one_random_gene:
                df.at[HOG, "udiv_genes"] = random.choice(ud_genes)
            else:
                df.at[HOG, "udiv_genes"] = ud_genes
                # print(ud_genes)
                dmel_OLs = []
                for i in range(len(ud_genes)):
                    q = ortholog_df[ortholog_df["Uloborus_diversus"] == ud_genes[i]][
                        "Drosophila_melanogaster"
                    ]
                    try:
                        dmel_genes = q.iloc[0].split(", ")
                        # print(dmel_genes)
                        dmel_OLs.extend(dmel_genes)
                    except:
                        pass

                df.at[HOG, "dmel_orthologs"] = dmel_OLs
        except:
            pass

    return df

def main(res_df, NX_tsv, id_converter):
    """Main function to process the results DataFrame and merge it with
    Uloborus diversus genes and their Drosophila melanogaster orthologs."""
    res_with_udiv_df = get_udiv_dmel_genes(
        NX_tsv,
        res_df,
        "/Users/calvin/kelvin-scratch/data/Uloborus_diversus__v__Drosophila_melanogaster.tsv",
    )

    res_with_udiv_df = res_with_udiv_df.explode("udiv_genes")
    res_with_udiv_df["udiv_genes"] = res_with_udiv_df["udiv_genes"].apply(
        lambda x: str(x).rsplit(".", 1)[0]
    )
    res_with_udiv_df = res_with_udiv_df.reset_index()

    merged_df = pd.merge(res_with_udiv_df, id_converter, how="left", on="udiv_genes")

    return merged_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("res_csv", help="Path to the results CSV file containing HOGs of interest")
    parser.add_argument("NX_tsv", help="Path to the hierarchical orthogroup file")
    args = parser.parse_args()

    id_converter = make_id_converter()

    main(args.res_csv, args.NX_tsv, id_converter)
