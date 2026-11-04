# Data directory contents

|Filename|Description|
|---|---|
|`GO.gmt`| Human GO GMT used as template for creating spider GO GMT files. |
|`N5.GeneCount.tsv`| Output of OrthoFinder's [orthogroup_gene_count.py](https://github.com/OrthoFinder/OrthoFinder/blob/main/tools/orthogroup_gene_count.py) tool run on the [N5.tsv](data/N5.tsv) OrthoFinder output (below) |
|`N5.tsv`| Orthogroups inferred by OrthoFinder at node 5, which includes only the species in [species_list_N5.txt](data/species_list_N5.txt).|
|`N5.udiv.o75_list.txt`| List of N5 orthogroups with occupancy of at least 75 and the presence of at least one gene from the spider _Uloborus diversus_. These 4756 HOGs were analysed in the HyPhy selection tests. |
|`N5_occ_30.tsv`| Orthogroups inferred by OrthoFinder at node 5, filtered to occupancy of at least 30 species. |
|`N5_occ_30_blasted_ptep.tsv`| Best BLAST hit in the *P. tepidariorum* genome for each orthogroup with occupancy ≥ 30; used for ontology enrichment. |
|`N5_occ_30_blasted_udiv.tsv`| Best BLAST hit in the *U. diversus* genome for each orthogroup with occupancy ≥ 30; used for ontology enrichment. |
|`SpeciesTree_full_brlen.nwk`| OrthoFinder-generated species tree with branch lengths of the 103 species used in this study ([species_list_all.txt](data/species_list_all.txt)).|
|`Uloborus_diversus__v__Drosophila_melanogaster.tsv`| 1:1 orthologs between _U. diversus_ and _D. melanogaster_ used in the helper script [src/id_converter.py](src/id_converter.py) for annotation of significant orthogroups.|
|`buscos.csv`| BUSCO results used to scale gene counts according to transcriptome redundancy and completeness in the [log odds ratio test](scripts/04_permulation_loss_dup).|
|`combined_GO.gmt`| GMT file containing GO terms for all *U. diversus* and *P. tepidariorum* genes, used as the reference for creating the network plot in [Figure 5](figures/figure_5); created by combining the `ptep_GO.gmt` and `udiv_GO.gmt` files below. |
|`id_converter.tsv`| Data compiled from _Uloborus diversus_ transcriptomes, genome and annotation files (Miller et al. 2023) used by the [src/id_converter.py](src/id_converter.py) to identify annotations for significant orthogroups.|
|`interesting_hits.txt`| List of HOGs with potential relevance to neural processes used to generate [Supplementary Table 6](results/Supplementary_Table_6_OddsRatioTest_interesting_hits.xlsx).|
|`non-orb-weavers-list.txt`| List of spiders designated as non-orb-weavers for this study.|
|`orbweavers-list.txt`| List of spiders designated as orb-weavers for this study.|
|`orthorun_list_fams_busco90.csv`| List of species of accessed transcriptomes filtered to those with at least 90% complete BUSCOs and no more than 30% duplicated BUSCOs.|
|`perms10000.RData`| RData object from the permulation run used in the log odds ratio test analysis.|
|`perms_tip_values.csv`| 10000 permulated phenotype vectors used in the log odds ratio test analysis.|
|`perms_tip_values_100k.csv.gz`, `perms_tip_values_500k.csv.gz`, `perms_tip_values_1M.csv.gz`| 100,000 / 500,000 / 1,000,000 permulated phenotype vectors (generated in batches by [`permulations_batch.R`](../scripts/04_permulation_loss_dup/permulations_batch.R)); the 500k and 1M sets extend the 100k set. Gzipped (the plain 1M CSV is 203 MB, over GitHub's file limit); `run_odds_ratio_test.py`, `power_analysis.py` and the other readers take the `.csv.gz` path directly. The plain `.csv` versions are git-ignored.|
|`perms_tip_values_nohubs_100k.csv.gz`, `perms_tip_values_nohubs_500k.csv.gz`| 100,000 / 500,000 permulated phenotype vectors on the tree without the species in `hub_species_to_exclude.txt` (`permulations_batch.R` with an exclusion file); used for the hub-removal odds ratio test runs `Run4`/`Run5` in `results/odds_ratio_test/Results_Nov02/`.|
|`hub_species_to_exclude.txt`| Nine species dropped in the species-removal robustness analysis of the log odds ratio test.|
|`orbweavers-list_nohubs.txt`, `non-orb-weavers-list_nohubs.txt`| The orb-weaver and non-orb-weaver lists without the species in `hub_species_to_exclude.txt`.|
|`flag_table_genes.csv`, `flag_table_cells.csv`| Suspected cross-sample contamination in the nine hub samples of PRJDB7399 (rules and caveats in [`flag_table_README.txt`](data/flag_table_README.txt)). The gene-level table has one row per flagged gene with its best near-identical match in another family and a `strength` label (strong, moderate or weak, from coding-sequence identity); the cell-level table has one row per HOG and hub species with the copies flagged and the copies remaining. Passed to the [log odds ratio test](scripts/04_permulation_loss_dup) with `--flagged-cells` (the gene-level table is the input; `--flag-strengths` picks the labels used).|
|`flag_table_README.txt`| How the flag tables were made, what each column means, and their caveats (thresholds are first choices; the direction of transfer is unknown).|
|`udiv_go_term_gene_sets.tsv`| GO term to *U. diversus* gene sets used by the permulation-based enrichment ([`permulation_enrichment.py`](../scripts/06_enrichment/permulation_enrichment.py)).|
|`ptep_GO.gmt`| GMT file containing GO terms for all *P. tepidariorum* genes; combined with `udiv_GO.gmt` for use as the reference for creating the network plot in [Figure 5](figures/figure_5). |
|`ptep_gene_ontology.gaf`| Gene annotation file from _P. tepidariorum_ (Zhu et al. 2023) used to annotate orthogroups. |
|`ptep_go_annots.all.tsv`| GO annotations for _P. tepidariorum_ gene IDs used in the ontology enrichment pipeline [scripts/06_enrichment](scripts/06_enrichment).|
|`silk_gland_genes_udiv_tblastn.csv`| Best BLAST hit for silk-gland assocated genes from _Argiope argentata_ (Chaw et al. 2021) in the _U. diversus_ genome.|
|`species_list_N5.txt`| List of 98 entelegyne spider species used in most analyses in this study.|
|`species_list_all.txt`| List of 102 spiders and outgroup species used in the OrthoFinder orthology search.|
|`spidroins_LOCs.tsv`| Gene IDs of _U. diversus_ silk genes.|
|`txptome_accessions.csv`| List of transcriptome accession IDs for ~900 spiders, some of which were not included in this study.|
|`udiv_GO.gmt`| GMT file containing GO terms for all *U. diversus* genes; combined with `ptep_GO.gmt` for use as the reference for creating the network plot in [Figure 5](figures/figure_5). |
|`udiv_gene_ontology.gaf`| Gene annotation file from _U. diversus_ (Miller et al. 2023) used to annotate orthogroups. |
|`udiv_go_annots.all.tsv`| GO annotations for _U. diversus_ gene IDs used in the ontology enrichment pipeline [scripts/06_enrichment](scripts/06_enrichment). |