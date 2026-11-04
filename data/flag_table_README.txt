Flag table for suspected cross-sample contamination in PRJDB7399 (built 2026-11-02, rebuilt by the scripts in scripts/04_permulation_loss_dup/flag_table/)

How it was made
  Queries: all 190,288 proteins of the nine "suspect samples" (Araneus inustus, Ozyptila nongae, Pseudomicrargus latitegulatus, Phrurolithus coreanus,
           Drapetisca socialis, Zelotes iriomotensis, Floronia bucculenta, Arcuphantes paiki, Latrodectus hasseltii).
  Search:  DIAMOND blastp (--faster, >=97% identity, >=70% query coverage, e<=1e-10) against all 2.04 million proteins of the 102 transcriptomes.
           Queries with no hit were re-scanned with low-complexity masking off (masking hides Gly/Ala-rich proteins).
  Kept:    matches to a species in a DIFFERENT family, protein identity >=98%, alignment >=100 aa. n_species_sharing counts the species with such a
           match to the gene (itself included); the best partner (highest identity, then longest alignment) is the one reported.
  Check:   for the best partner of each gene, coding-sequence identity from a codon alignment (protein alignment mapped back to the CDS).

flag_table_genes.csv  (one row per flagged suspect-sample gene; 16,063 rows)
  species, gene_id, HOG            the suspect sample's gene and its N5 HOG
  partner_species/gene/family      its best near-identical match in another family; partner_is_suspect = the partner is also one of the nine suspect samples
  protein_identity, aligned_aa     from DIAMOND
  nt_identity, aligned_codons      coding-sequence identity over aligned codons
  n_species_sharing                number of species (incl. this one) with a >=98% match to this gene
  strength                         strong: nt_identity >= 99 and n_species_sharing <= 20 (8,852 genes)
                                   moderate: nt_identity 97-99 and <= 20 species, or nt_identity >= 99 and > 20 species (1,735 genes)
                                   weak: everything else (5,476 genes)
                                   (Above about 20 species the medians fall to ordinary conserved-gene levels, so wide sharing is not treated as strong evidence.)

flag_table_cells.csv  (one row per HOG x suspect species that has >=1 flagged gene; 12,350 rows)
  n_copies                         genes of that species in the HOG
  n_flagged_strong / _or_moderate  flagged genes at each level
  remaining_if_..._removed         copies left if flagged genes are removed
  cell_status_strong / _strong_moderate   missing = no copies left (treat as MISSING, not absent) | reduced = fewer copies | ok
  (strong flags alone: 6,675 cells in 2,301 HOGs, 5,920 of them missing)

Caveats
  Flags are not proof of contamination; the direction of transfer is unknown (only the suspect-side gene is flagged; the partner may be the source).
  Species that are not suspect samples are only flagged through being a partner (see partner_species), not scanned themselves.
  Thresholds (98% protein, 99%/97% nucleotide, <=20 species) are first choices, not yet calibrated per species pair with BUSCO orthologs.
