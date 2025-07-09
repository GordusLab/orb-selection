#!/usr/bin/env python3

"""Get silk gland genes from a list of LOCs."""

import os
import argparse
import pandas as pd
from src.id_converter import make_id_converter


def main(LOC_list, top10=False):
    """Returns silk gland genes from a list of LOCs, prints the number of overlapping LOCs."""

    scripts = os.path.dirname(__file__)
    assets = os.path.join(scripts, "..", "assets")

    with open(LOC_list) as f:
        LOCs = f.read().splitlines()

    silk_genes_udiv_blast_df = pd.read_csv(
        f"{assets}/silk_gland_genes_udiv_tblastn.csv",
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

    id_converter = make_id_converter()

    silk_genes_udiv = id_converter["txpt"].isin(silk_genes_udiv_blast_df["txpt"])
    silk_LOCs = id_converter[silk_genes_udiv]["LOC"].dropna().unique()

    intersect_list = list(set(silk_LOCs) & set(LOCs))

    print(f"{len(intersect_list)} overlapping LOCs")

    return intersect_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "LOC_list", type=str, help="List of LOCs to intersect with silk genes"
    )
    parser.add_argument(
        "--top10",
        action="store_true",
        help="Whether to use the single best blast hit or top 10 hits",
    )
    args = parser.parse_args()
    main(args.LOC_list, args.top10)
