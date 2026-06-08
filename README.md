# orb-selection

Analysis repository for "Comparative transcriptomic analysis reveals signatures of selection
for orb-weaving behavior in spiders," Runnels et al. 2026

## Repository contents

|Directory / Filename|Description|
|---|---|
|`data/`| Raw and intermediate data files used in the pipeline, _e.g._ OrthoFinder outputs, BUSCO scores, lists of species belonging to different categories, lists of HOGs tested in the HyPhy analyses, permulation-generated phenotype designations, and resources used to annotate results.|
|`figures/`| PDF files of all figures and figure elements output by the figure-generating notebooks in [`scripts/07_figures_tables`](scripts/07_figures_tables).|
|`results/`| Results from the analyses including all Supplementary Tables ([source code](scripts/06_figures_tables/Supplementary%20Data%20Tables.ipynb)), lists of _U. diversus_ or _P. tepidariorum_ gene IDs for significant HOGs and complete GO enrichment of these genes (Fig. 3-5, [source code](scripts/05_enrichment)), and all outputs from the Log Odds Ratio test (Fig. 5, [source code](scripts/04_permulation_loss_dup)) and phyloGLM. |
|`scripts/`| Full analysis pipeline divided into stages following the paper's methods section. See below for a complete description of the steps required for each stage of the analysis.|
|`src/` | Contains helper modules used in various stages of the analysis. |
|`README.md`| Top-level project overview, pipeline stage documentation, and usage notes.|
|`CHANGELOG.md`|Repository reorganization details. Reorganization support was performed with GitHub Copilot (GPT-5.3-Codex). |
|`pyproject.toml`| Python packaging and project metadata configuration (build system, dependencies, package discovery).|
|`.Rprofile`|Project-level R startup settings.|
|`.gitignore`|Version control exclusion rules.|
|`environment.yml`|Python dependencies.|
|`renv.lock`|R dependencies.|

## Important Data Availability Notes

Some required inputs and outputs are not available in this repository.

- Raw and processed transcriptome FASTA inputs for OrthoFinder and the HyPhy pipeline are external. `TODO_ADD_TXTOME_DATA_ARCHIVE_LOCATION`
- HyPhy per-gene JSON outputs are external. `TODO_ADD_HYPHY_JSON_ARCHIVE_LOCATION`
- Some cache pickle files for HyPhy and Log Odds Ratio test results used by downstream scripts are external. `TODO_ADD_PKL_ARCHIVE_LOCATION`

## Stage 01: Transcriptome data collection and pre-processing

The scripts in [`scripts/01_pre_processing`](scripts/01_pre_processing) prepare the raw sequence set used in the orthology search.

### Steps: 
1. Download transcriptomes ([Supplementary Table 1](results/Supplementary_Table_1_SpiderAccessions_BUSCOs.xlsx)) from NCBI GenBank.
2. Cluster with CD-HIT: [`process_fsas_cd-hit.sh`](scripts/01_pre_processing/process_fsas_cd-hit.sh)
3. Identify open reading frames with TransDecoder.LongOrfs: [`process_fsas.TD-LO.sh`](scripts/01_pre_processing/process_fsas.TD-LO.sh)
4. Identify homology with BLAST-P: [`process_fsas.blastp.sh`](scripts/01_pre_processing/process_fsas.blastp.sh)
5. Predict coding sequences with TransDecoder.Predict: [`process_fsas.TD-P.sh`](scripts/01_pre_processing/process_fsas.TD-P.sh)
6. Analyze transcriptome quality using BUSCO: [`process_fsas.busco.sh`](scripts/01_pre_processing/process_fsas.busco.sh)

## Stage 02: Orthology search and pre-testing pipeline

The scripts in [`scripts/02_orthofinder_prep_hyphy`](scripts/02_orthofinder_prep_hyphy) process OrthoFinder results and corresponding sequences to prepare them for selection testing. OrthoFinder was re-run several times due to errors related to the number of files and the need to use an edited species tree; a record of the options used each time is in the [stage 2 README.md](scripts/02_orthofinder_prep_hyphy/README.md).

### Steps: 
1. Run OrthoFinder (see [stage 2 README.md](scripts/02_orthofinder_prep_hyphy/README.md))
2. Filter orthogroups to N5 HOGs with occupancy ≥ 75 and get nucleotide sequences: [`get_nuc_seqs.sh`](scripts/02_orthofinder_prep_hyphy/get_nuc_seqs.sh)
3. Preparation for HyPhy testing: [`prep_for_hyphy.sh`](scripts/02_orthofinder_prep_hyphy/prep_for_hyphy.sh). Includes:
   - Pre-alignment quality filtering with PREQUAL
   - Alignment with MACSE 
   - Gene tree generation with IQ-TREE 
   - Error-filtering alignments and trees with BUSTED + hyphy error-filter
4. Foreground branch labeling: [`label_trees.sh`](scripts/02_orthofinder_prep_hyphy/label_trees.sh)

Data note: Stage 02 uses processed FASTA files not tracked in GitHub due to size.

## Stage 03: Testing for positive and relaxed selection

The scripts/modules in [`scripts/03_selection_tests`](scripts/03_selection_tests/) run the HyPhy analyses on the 4,576 N5 orthogroups and parse the results. 

### Steps: 

1. Run RELAX: [`relax.sh`](scripts/03_selection_tests/relax.sh)
2. Run BUSTED-PH with orb-weavers as foreground: [`busted_ph.sh`](scripts/03_selection_tests/busted_ph.sh)
3. Run BUSTED-PH with non-orb-weavers as foreground: [`busted_ph_switch_fg.sh`](scripts/03_selection_tests/busted_ph_switch_fg.sh)
4. Parse the results: [`parse_hyphy_results.py`](scripts/03_selection_tests/parse_hyphy_results.py)

Data note: Stage 03 uses thousands of JSON result files not tracked in GitHub due to size.

## Stage 04: Gene loss and duplication analysis

The scripts and notebooks in [`scripts/04_permulation_loss_dup`](scripts/04_permulation_loss_dup) run the odds-ratio permulation workflow and evaluate distribution shape assumptions for test statistics.

### Steps:
1. Generate and export permulation tip assignments: [`permulations.R`](scripts/04_permulation_loss_dup/permulations.R)
2. Run loss/duplication odds-ratio test module [`odds_ratio_test.py`](scripts/04_permulation_loss_dup/odds_ratio_test.py) using the workflow in the [`Odds Ratio Permulation Test.ipynb`](scripts/04_permulation_loss_dup/Odds%20Ratio%20Permulation%20Test.ipynb) notebook

Data note: Stage 04 uses some cached pickle inputs/outputs that are external to this repository due to size.

## Stage 05: Phylogenetic regression

The script and notebook in [`scripts/05_phyloglm`](scripts/05_phyloglm) fit a phylogenetic generalized linear model to each gene and evaluate the results. 

### Steps: 
1. Run phyloglm analysis in parallel on ~12000 genes: [`phyloglm.R`](scripts/05_phyloglm/phyloglm.R)
2. Inspect significant results: [`PhyloGLM Analysis.ipynb`](scripts/05_phyloglm/PhyloGLM%20Analysis.R)

## Stage 06: Ontology enrichment analysis

The scripts in [`scripts/06_enrichment`](scripts/06_enrichment) create significant gene ID lists from HyPhy and `odds_ratio_test.py` results and run GO enrichment summaries.

### Steps:
1. Make BLAST db from the [_P. tepidariorum_ genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_043381705.1/)
2. Run [`annotate_ogroups_vs_ref.py`](scripts/06_enrichment/annotate_ogroups_vs_ref.py) to determine best BLAST hit for each orthogroup from the _P. tepidariorum_ for enrichment of significant gene sets for orthogroups more commonly found in non-orb-weavers
3. Helper functions to generate significant ID list files for downstream enrichment: [`get_gene_id_lists.py`](scripts/05_enrichment/get_gene_id_lists.py)
3. Run topGO enrichment for each gene set: [`go_enrichment.R`](scripts/05_enrichment/go_enrichment.R)
4. Summarize enrichment outputs into merged tables: [`summarise_topgo_output.sh`](scripts/05_enrichment/summarise_topgo_output.sh)

Data note: Stage 06 depends on cached HyPhy and odds ratio test result objects for hit-list generation.

## Stage 07: Figures and tables

The scripts and notebooks in [`scripts/07_figures_tables`](scripts/07_figures_tables) generate manuscript figures, significant results intersections, and supplementary data tables.

### Steps:

1. Plot HyPhy omega distributions and selected gene examples: [`Hyphy Omega Plots.ipynb`](scripts/06_figures_tables/Hyphy%20Omega%20Plots.ipynb)
2. Plot odds ratio test results: [`Odds Ratio Test Plots.ipynb`](scripts/06_figures_tables/Odds%20Ratio%20Test%20Plots.ipynb)
3. Generate UpSet plots and intersections of significant results: [`UpSet Plots.ipynb`](scripts/06_figures_tables/UpSet%20Plots.ipynb)
4. Compile supplementary result tables for export: [`Supplementary Data Tables.ipynb`](scripts/06_figures_tables/Supplementary%20Data%20Tables.ipynb)

Data note: Stage 07 expects completed outputs from stages 03-05 and writes figure/table artifacts to `figures/` and `results/`.
