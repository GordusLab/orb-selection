# Data directory contents

|Filename|Description|
|---|---|
|`N5.GeneCount.tsv`| Output of OrthoFinder's [orthogroup_gene_count.py](https://github.com/OrthoFinder/OrthoFinder/blob/main/tools/orthogroup_gene_count.py) tool run on the [N5.tsv](data/N5.tsv) OrthoFinder output (below) |
|`N5.tsv`| Orthogroups inferred by OrthoFinder at node 5, which includes only the species in [species_list_N5.txt](data/species_list_N5.txt).|
|`N5.udiv.o75_list.txt`| List of N5 orthogroups with occupancy of at least 75 and the presence of at least one gene from the spider _Uloborus diversus_. These 4756 HOGs were analysed in the HyPhy selection tests. |
|`N5_blasted.tsv`| [N5.tsv](data/N5.tsv) orthogroups best BLAST hits from the _Parasteatoda tepidariorum_ genome (Zhu et al. 2023).|
|`SpeciesTree_full_brlen.nwk`| OrthoFinder-generated species tree with branch lengths of the 103 species used in this study ([species_list_all.txt](data/species_list_all.txt)).|
|`Uloborus_diversus__v__Drosophila_melanogaster.tsv`| 1:1 orthologs between _U. diversus_ and _D. melanogaster_ used in the helper script [src/id_converter.py](src/id_converter.py) for annotation of significant orthogroups.|
|`buscos.csv`| BUSCO results used to scale gene counts according to transcriptome redundancy and completeness in the [log odds ratio test](scripts/04_permulation_loss_dup).|
|`id_converter.tsv`| Data compiled from _Uloborus diversus_ transcriptomes, genome and annotation files (Miller et al. 2023) used by the [src/id_converter.py](src/id_converter.py) to identify annotations for significant orthogroups.|
|`interesting_hits.txt`| List of HOGs with potential relevance to neural processes used to generate [Supplementary Table 6](results/Supplementary_Table_6_OddsRatioTest_interesting_hits.xlsx).|
|`non-orb-weavers-list.txt`| List of spiders designated as non-orb-weavers for this study.|
|`orb-weavers-list.txt`| List of spiders designated as orb-weavers for this study.|
|`orthorun_list_fams_busco90.csv`| List of species of accessed transcriptomes filtered to those with at least 90% complete BUSCOs and no more than 30% duplicated BUSCOs.|
|`perms10000.RData`| RData object from the permulation run used in the log odds ratio test analysis.|
|`perms_tip_values.csv`| 10000 permulated phenotype vectors used in the log odds ratio test analysis.|
|`ptep_gene_ontology.gaf`| Gene annotation file from _P. tepidariorum_ (Zhu et al. 2023) used to annotate orthogroups.|
|`ptep_go_annots.all.tsv`| GO annotations for _P. tepidariorum_ gene IDs used in the ontology enrichment pipeline [scripts/06_enrichment](scripts/06_enrichment).|
|`silk_gland_genes_udiv_tblastn.csv`| Best BLAST hit for silk-gland assocated genes from _Argiope argentata_ (Chaw et al. 2021) in the _U. diversus_ genome.|
|`species_list_N5.txt`| List of 98 entelegyne spider species used in most analyses in this study.|
|`species_list_all.txt`| List of 102 spiders and outgroup species used in the OrthoFinder orthology search.|
|`spidroins_LOCs.tsv`| Gene IDs of _U. diversus_ silk genes.|
|`txptome_accessions.csv`| List of transcriptome accession IDs for ~900 spiders, some of which were not included in this study.|
|`udiv_go_annots.all.tsv`| GO annotations for _U. diversus_ gene IDs used in the ontology enrichment pipeline [scripts/06_enrichment](scripts/06_enrichment). |