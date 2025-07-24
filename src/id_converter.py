"""
This script processes gene orthology data to find Uloborus diversus genes
and their Drosophila melanogaster orthologs, and merges this information with
a gene ID converter to provide additional gene descriptions."""

import argparse
import random
import os
import pandas as pd

src_path = os.path.dirname(__file__)
data = os.path.join(src_path, "..", "data")


def make_id_converter():
    """Creates a DataFrame that maps Uloborus diversus transcript IDs to
    LOCs and nucleotide accession IDs."""

    id_df = pd.read_csv(
        f"{data}/id_converter.tsv",
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

    id_df["desc_index_start"] = id_df["gene_desc"].str.find("[Description: ") + 14
    id_df["desc_index_end"] = id_df["gene_desc"].str.find("] [Gene Type:")

    id_df["Description"] = id_df.apply(
        lambda row: row.gene_desc[row.desc_index_start : row.desc_index_end], axis=1
    )
    id_df = id_df[["LOC", "txpt", "udiv_genes", "GO_terms", "Description"]]

    return id_df


def get_udiv_dmel_genes(hog_node_genes_tsv, hogs_of_interest, ortholog_tsv, one_random_gene=False):
    """Retrieves Uloborus diversus genes and their Drosophila melanogaster orthologs
    for a given list (csv or df) of hogs. The orthogroup file for the hierarchical orthogroup node
    of interest must be provided."""

    hog_node_df = pd.read_csv(hog_node_genes_tsv, sep="\t", index_col="HOG", dtype=str)

    ortholog_df = pd.read_csv(ortholog_tsv, sep="\t", dtype=str)
    ortholog_df["Uloborus_diversus"] = ortholog_df["Uloborus_diversus"].str.split(", ")
    ortholog_df = ortholog_df.explode("Uloborus_diversus")

    # Attempt to read the hog DataFrame, if it fails, assume it's already a DataFrame
    # and assign empty columns for udiv_genes and dmel_orthologs
    try:
        df = pd.read_csv(hogs_of_interest, index_col="HOG", dtype=str)
        df = df.assign(udiv_genes="", dmel_orthologs="")
    except TypeError:
        df = hogs_of_interest
        df = df.assign(udiv_genes="", dmel_orthologs="")

    # Iterate through each hog and find corresponding Uloborus diversus genes
    # and their Drosophila melanogaster orthologs
    # If one_random_gene is True, select one random gene per orthogroup
    # Otherwise, list all genes and their orthologs
    for hog in df.index.to_list():
        try:
            ud_genes = hog_node_df.at[hog, "Uloborus_diversus"].split(", ")

            # use if i just want to find 1 U.div gene per OG
            if one_random_gene:
                df.at[hog, "udiv_genes"] = random.choice(ud_genes)
            else:
                df.at[hog, "udiv_genes"] = ud_genes
                # print(ud_genes)
                dmel_orthologs = []
                for gene in ud_genes:
                    q = ortholog_df[ortholog_df["Uloborus_diversus"] == gene][
                        "Drosophila_melanogaster"
                    ]
                    try:
                        dmel_genes = q.iloc[0].split(", ")
                        # print(dmel_genes)
                        dmel_orthologs.extend(dmel_genes)
                    except (IndexError, AttributeError):
                        pass

                df.at[hog, "dmel_orthologs"] = dmel_orthologs
        except KeyError:
            pass

    return df


def main(res_df, hog_node_genes_tsv):
    """Main function to process the results DataFrame and merge it with
    Uloborus diversus genes and their Drosophila melanogaster orthologs."""

    id_converter_df = make_id_converter()

    res_with_udiv_df = get_udiv_dmel_genes(
        hog_node_genes_tsv,
        res_df,
        f"{data}/Uloborus_diversus__v__Drosophila_melanogaster.tsv"
    )

    res_with_udiv_df = res_with_udiv_df.explode("udiv_genes")
    res_with_udiv_df["udiv_genes"] = res_with_udiv_df["udiv_genes"].apply(
        lambda x: str(x).rsplit(".", 1)[0]
    )
    res_with_udiv_df = res_with_udiv_df.reset_index()

    merged_df = pd.merge(res_with_udiv_df, id_converter_df, how="left", on="udiv_genes")

    return merged_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "results", 
        help="Path to the results CSV file containing hogs of interest, \
            OR a DataFrame with HOGs as index"
    )
    parser.add_argument(
        "hog_node_genes_tsv", 
        help="Path to the hierarchical orthogroup file")
    args = parser.parse_args()

    main(args.results, args.hog_node_genes_tsv)
