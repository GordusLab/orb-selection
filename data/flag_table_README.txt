Flag table for suspected cross-sample contamination in PRJDB7399 (built 2026-11-02)

How it was made
  Queries: all 190,288 proteins of the nine "hub" samples (Araneus inustus, Ozyptila nongae, Pseudomicrargus latitegulatus, Phrurolithus coreanus,
           Drapetisca socialis, Zelotes iriomotensis, Floronia bucculenta, Arcuphantes paiki, Latrodectus hasseltii).
  Search:  DIAMOND blastp (--faster, >=97% identity, >=70% query coverage, e<=1e-10) against all 2.04 million proteins of the 102 transcriptomes.
  Kept:    matches to a species in a DIFFERENT family, protein identity >=98%, alignment >=100 aa, and the sequence shared by <=6 species in all
           (so proteins shared by many species, as real ultra-conserved genes are, are excluded).
  Check:   for the best partner of each gene, coding-sequence identity from a codon alignment (protein alignment mapped back to the CDS).

flag_table_genes.csv  (one row per flagged hub gene; 7,822 rows)
  species, gene_id, HOG            the hub sample's gene and its N5 HOG
  partner_species/gene/family      its best near-identical match in another family; partner_is_hub = the partner is also one of the nine hubs
  protein_identity, aligned_aa     from DIAMOND
  nt_identity, aligned_codons      coding-sequence identity over aligned codons
  n_species_sharing                number of species (incl. this one) with a >=98% match to this gene
  strength                         strong: nt_identity >= 99 | moderate: 97-99 | weak: < 97

flag_table_cells.csv  (one row per HOG x hub species that has >=1 flagged gene; 5,998 rows)
  n_copies                         genes of that species in the HOG
  n_flagged_strong / _or_moderate  flagged genes at each level
  remaining_if_..._removed         copies left if flagged genes are removed
  cell_status_strong / _strong_moderate   missing = no copies left (treat as MISSING, not absent) | reduced = fewer copies | ok

ort_run3_hog_summary.csv           per-HOG results of the 1,000,000-permulation run (copied from the repo, read-only) with the two significance labels

Caveats
  Flags are not proof of contamination; the direction of transfer is unknown (only the hub-side gene is flagged; the partner may be the source).
  Species that are not hubs are only flagged through being a partner (see partner_species), not scanned themselves.
  Thresholds (98% protein, 99%/97% nucleotide, <=6 species) are first choices, not yet calibrated per species pair with BUSCO orthologs.
