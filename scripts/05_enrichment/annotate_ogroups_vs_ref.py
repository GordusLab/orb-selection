#!/usr/bin/env python3

from Bio.Blast.Applications import NcbitblastnCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML 
import pandas as pd 
import sys 


if len(sys.argv) < 3 : 
    print( 'annotate_ogroups_vs_ref.py  N0.tsv  ortho_dir/ ref_genes.fasta \n ortho_dir contains fastas for all species in N0. Also, do this first: \n makeblastdb  -dbtype nucl -in reference_genes.fasta ')
    exit()

df  = pd.read_csv(sys.argv[1], sep="\t")
out = pd.read_csv(sys.argv[1], sep="\t", usecols =["HOG", "OG"] )

out["Best_hit"] = "NO_HIT"
out["E_val"]    = "NO_HIT"
out["iden"]     = "NO_HIT"
out["len"]      = "NO_HIT"
out["Coords"]   = "NO_HIT"

# Go through relevant groups of Orthofinder N0
for row in range(0, len(df)):

    # Orthogroup in question
	HOG = df.iloc[row,0]

    # Detect species where HOG exists   
	cols = df.iloc[row,3:].notna()
	idxs =[ i for i, val in enumerate(cols) if val == True ]
	cols = list(cols[idxs].index)

    # Prepare loops
    i = 0     # column iterator
    hit = False
    
	while hit == False:
	
		if i < len(cols):
        
			# Change species while there is no blast hit 
			species = cols[i]

			# Get first gene of that species in Ogroup
			gene = df.loc[row,species].split(",")[0]
        
            ### Extract gene from species fasta #############
			infasta = f"{sys.argv[2]}/{species}.fasta"
			for seq_record in SeqIO.parse(infasta, "fasta"):

				if seq_record.id == gene:
                
					# Extract gene from corresponding species
					SeqIO.write(seq_record, "tmp.fasta", "fasta")

					# Blast gene against reference db (needs makeblastdb!)
					blastx_cline = NcbitblastnCommandline(query="tmp.fasta", db = sys.argv[3],
                                         evalue=1e-5,  outfmt=5, out="my_blast.xml") 

					stdout, stderr = blastx_cline()
  
					if len (stderr) > 0: print("BLAST ERROR: ", stderr)

					for record in NCBIXML.parse(open("my_blast.xml")): 
						if record.alignments : #skip queries with no matches 

							# record.alignments is list of alignments ordered by e-val. I only want the 1st elem now
							align = record.alignments[0]
							# .hsps[0] gets alignment and .expect the evalue
							hsps = align.hsps[0] 
							e_val= hsps.expect

							# Multiply each prot aminoacid by 3 to get numbers in  nucleotides
							iden = hsps.identities * 3
							leng = hsps.align_length * 3 


							out.loc[row, "Best_hit"] = align.title
							out.loc[row, "E_val"] = e_val
							out.loc[row, "iden"] = iden
							out.loc[row, "len"] =  leng 
							out.loc[row, "Coords"] = str( [hsps.sbjct_start, hsps.sbjct_end] )
                        
							hit = True   # End while loop
			# Update column iterator for next while loop
			i += 1

		else: hit = True                        
                    
# Output will have same name as N0 input plus "_blasted.tsv" suffix
outname = sys.argv[1][:-4] + "_blasted.tsv"               
out.to_csv(outname, sep="\t", index= False)


