#!/usr/bin/env python3

# Run in output directory

import numpy as np
import pandas as pd
import subprocess
import sys

if len(sys.argv) != 4 or sys.argv[1] == "-h" or sys.argv[1] == "-help" or sys.argv[1] == "--help":
    print("Usage: get_nuc_seqs_parallel.py <NX.tsv> <HOG> <CDs-DIR>")
    sys.exit()

# get df of gene names
NX_df = pd.read_csv(sys.argv[1], sep='\t', index_col="HOG", dtype=str)
NX_df = NX_df.drop(columns=['OG', 'Gene Tree Parent Clade'])

HOG=sys.argv[2]

directory = sys.argv[3]


# make a fasta file for each orthogroup in desired list of OGs
def make_fasta(HOG):
	file = open("%s.temp.fasta" % HOG, "w")

	# for each species, make a list of all gene IDs for that species
	for species in NX_df.columns:
		genes = str(NX_df.at[HOG, species])

		# ignore species with no genes in the OG
		if genes == 'nan':
			pass
		else: 
			# print(species)

			# run seqkit grep in the CDS file for that species on all genes in the 
			# OG and add their sequences to the new fasta file
			genes=genes.split(", ")
			
			for gene in genes: 
				# print(gene)
				
				seq = subprocess.run(['seqkit', 'grep', '-p', '%s' % gene, '/%s/%s.cd-hit-est.transdecoder.cds' % (directory, species)], capture_output=True, text=True)

				# print(seq.stdout)

				# get rid of the > so i can add it before the species prefix
				seq_no_carat = seq.stdout.split(">")[1]

				# print(seq_no_carat)

				# add species name as prefix
				file.writelines(">" + species + "|")

				# write gene name and sequence
				file.writelines(seq_no_carat)

	file.close()

	# clean up gene names
	subprocess.run(['cut -d " " -f1 %s.temp.fasta > %s.fasta' % (HOG, HOG)], shell=True)
	subprocess.run(['rm', '%s.temp.fasta' % HOG])

	print("HOG %s\tDONE" % HOG)
 
if __name__ == "__main__":
	make_fasta(HOG)


