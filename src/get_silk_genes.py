#!/usr/bin/env python3

"""Get silk gland genes from a list of LOCs."""

import os
import argparse
import pandas as pd
from id_converter import make_id_converter

scripts = os.path.dirname(__file__)
data_dir = os.path.join(scripts, "..", "..", "data")

def get_all_silk_genes(top10=False) -> list:
    """Returns all silk gland genes from the id_converter."""
    
    id_converter = make_id_converter()
    silk_genes_udiv_blast_df = pd.read_csv(
        f"{data_dir}/silk_gland_genes_udiv_tblastn.csv",
        usecols=["query acc.ver", " subject acc.ver", " % identity", " evalue"],
    )

    # If top10 is False, keep only the best blast hit for each query
    if not top10:
        silk_genes_udiv_blast_df = silk_genes_udiv_blast_df.drop_duplicates(
            subset=["query acc.ver"]
        )

    # Extract the transcript ID from the subject accession version
    # and remove the version number (e.g., .1)
    silk_genes_udiv_blast_df["txpt"] = silk_genes_udiv_blast_df[
        " subject acc.ver"
    ].apply(lambda x: str(x).replace(".1", ""))

    silk_genes_udiv = id_converter["txpt"].isin(silk_genes_udiv_blast_df["txpt"])
    silk_locs = id_converter[silk_genes_udiv]["LOC"].dropna().unique()

    return list(silk_locs)

def main(loc_list, top10=False):
    """Returns silk gland genes from a list of LOCs, prints the number of overlapping LOCs."""

    with open(loc_list, encoding="utf-8") as f:
        locs = f.read().splitlines()

    silk_locs = get_all_silk_genes(top10=top10)

    intersect_list = list(set(silk_locs) & set(locs))

    print(f"{len(intersect_list)} overlapping LOCs")

    return intersect_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "loc_list", type=str, help="List of LOCs to intersect with silk genes"
    )
    parser.add_argument(
        "--top10",
        action="store_true",
        help="Whether to use the single best blast hit or top 10 hits",
    )
    args = parser.parse_args()
    main(args.loc_list, args.top10)
