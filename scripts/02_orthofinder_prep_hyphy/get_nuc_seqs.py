#!/usr/bin/env python3

"""Generate one nucleotide FASTA file for a given HOG.

Usage:
	get_nuc_seqs.py <N5.tsv> <HOG_ID> <TRANSDECODER_DIR>
"""

import pandas as pd
import subprocess
import sys

if len(sys.argv) != 4 or sys.argv[1] == "-h" or sys.argv[1] == "-help" or sys.argv[1] == "--help":
    print("Usage: get_nuc_seqs_parallel.py <NX.tsv> <HOG> <CDs-DIR>")
    sys.exit()

# get df of gene names
NX_df = pd.read_csv(sys.argv[1], sep='\t', index_col="HOG", dtype=str)
NX_df = NX_df.drop(columns=['OG', 'Gene Tree Parent Clade'])

HOG = sys.argv[2]

directory = sys.argv[3]


# make a fasta file for each orthogroup in desired list of OGs
def make_fasta(hog_id):
    """Build a FASTA file for all species-associated genes in one HOG."""

    with open(f"{hog_id}.temp.fasta", "w", encoding="utf-8") as file:
        # For each species, make a list of all gene IDs for that species.
        for species in NX_df.columns:
            genes = str(NX_df.at[hog_id, species])

            # Ignore species with no genes in the HOG.
            if genes == "nan":
                continue

            # Run seqkit grep in the CDS file for that species on all genes in the HOG
            # and add their sequences to the new FASTA file.
            genes = genes.split(", ")

            for gene in genes:
                cds_path = f"{directory}/{species}.cd-hit-est.transdecoder.cds"
                seq = subprocess.run(
                    ["seqkit", "grep", "-p", f"{gene}", cds_path],
                    capture_output=True,
                    text=True,
                    check=False,
                )

                if ">" not in seq.stdout:
                    continue

                # Remove the leading '>' so we can prepend species name.
                seq_no_carat = seq.stdout.split(">", 1)[1]

                # Add species name as prefix.
                file.writelines(">" + species + "|")

                # Write gene name and sequence.
                file.writelines(seq_no_carat)

    # Clean up gene names.
    subprocess.run(
        [f'cut -d " " -f1 {hog_id}.temp.fasta > {hog_id}.fasta'],
        shell=True,
        check=True,
    )
    subprocess.run(["rm", f"{hog_id}.temp.fasta"], check=True)

    print(f"HOG {hog_id}\tDONE")
 
if __name__ == "__main__":
    make_fasta(HOG)


