# orb-selection

Analysis repository for "Comparative transcriptomic analysis reveals signatures of selection
for orb-weaving behavior in spiders," Runnels et al. 2026

TODO:

- [ ] add N5 list to assets
- [ ] edit slurm comments in rockfish scripts
- [ ] clear out assets directory
- [ ] edit stages 1 and 2
- [ ] read through and edit all updated docstrings

## What Is In This Repository

- Stage scripts for preprocessing, HyPhy analyses, permulation testing, and GO enrichment.
- Source/helper modules in `src/`.
- Figure/table notebooks and plotting code.
- Supplementary output tables under `results/`.

## Important Data Availability Notes

Some required inputs and outputs are not available in this repository.

- Raw transcriptome FASTA inputs for OrthoFinder and the HyPhy pipeline are external.
- HyPhy per-gene JSON outputs are external.
- Some `.pkl` cache/result files used by downstream scripts are external.

Placeholders for final archive links:

- Raw transcriptomes and other large stage-1/2 inputs: `TODO_ADD_RAW_DATA_ARCHIVE_LOCATION`
- HyPhy JSON result archive: `TODO_ADD_HYPHY_JSON_ARCHIVE_LOCATION`
- Large `.pkl` cache files archive: `TODO_ADD_PKL_ARCHIVE_LOCATION`

Git tracking note for large files/directories:

- The following raw-result data directories are intentionally gitignored:
	- `data/absrel_rerun/`
	- `data/busted_ph/`
	- `data/busted_ph_rev/`
	- `data/relax/`
- `.pkl` files are intentionally gitignored repository-wide.

## Stage 01: Transcriptome Data Acquisition and Initial Input Set

This stage corresponds to the first methods section and prepares the raw sequence set used downstream.

Use methods-described software/workflow for your run (outside this repo where needed):

- obtain transcriptome assemblies using accession metadata
- perform any stated acquisition/QC steps from methods section 1
- collect final FASTA inputs for comparative analyses

Raw accessions:

- Transcriptome accessions used in this project are listed in `results/Supplementary_Table_1_SpiderAccessions_BUSCOs.xlsx`.

Expected output of Stage 01:

- A local/external directory of transcriptome FASTA files ready for orthology and phylogenetic stages.

## Stage 02: Orthology and Comparative Input Preparation

This stage corresponds to the second methods section and prepares orthogroup/species-tree inputs used by downstream scripts.

Use methods-described software/workflow for your run (outside this repo where needed):

- OrthoFinder (or equivalent orthology pipeline in your methods)
- species tree preparation/editing steps described in methods
- any BUSCO/ortholog mapping steps described in methods

Expected output of Stage 02:

- Orthogroup-derived input tables/files required by stage 03+
- Species tree + tip/foreground labeling inputs
- Working directories consumed by the stage scripts below

## Stage 03: Phylogenetic Prep

Purpose:

- PREQUAL filtering
- MACSE alignment
- IQ-TREE trees
- initial BUSTED + error-filtered outputs
- foreground branch labeling
- prerequisite alignment/tree files for HyPhy

Scripts to run from `scripts/03_phylo_prep/`:

- `prep_for_hyphy.sh`
- `label_trees.sh`
- `get_nuc_seqs.sh`

External programs used in this stage:

- PREQUAL
- MACSE
- IQ-TREE
- seqkit/csvtk (if used in your methods pipeline)

## Stage 04: Selection Tests and HyPhy Result Parsing

Purpose:

- run HyPhy analyses for selection tests
- parse large HyPhy JSON outputs into usable result objects/tables

Scripts/modules in `scripts/04_selection_tests/`:

- `busted_ph.sh`
- `busted_ph_switch_fg.sh`
- `relax.sh`
- `hyphy_results_parser.py`
- `hyphy_results_helpers.py`
- `parse_hyphy_results.py`
- `README_hyphy_module.md`

Data note:

- Stage 04 uses thousands of JSON result files not tracked in GitHub due to size.
- Add your archive/location to: `TODO_ADD_HYPHY_JSON_ARCHIVE_LOCATION`.

## Stage 05: Permulation/Odds-Ratio Testing

Purpose:

- run the odds-ratio permulation testing workflow
- generate/consume permulated phenotype assignments

Scripts in `scripts/05_permulation_loss_dup/`:

- `permulations.R`
- `odds_ratio_test.py`
- `README_odds_ratio_test.md`
- `Odds Ratio Permulation Test.ipynb`
- `Multimodal Test.ipynb`

Notes:

- `odds_ratio_test.py` moved to this stage and is documented in `README_odds_ratio_test.md`.
- The notebook `Odds Ratio Permulation Test.ipynb` demonstrates how to run the module.
- Some `.pkl` outputs/read inputs are external due to size limits.

## Stage 06: Significant Gene ID Lists and GO Enrichment

Purpose:

- generate significant gene ID lists from odds-ratio/HyPhy cached results
- run and summarize GO enrichment

Scripts in `scripts/06_enrichment/`:

- `get_gene_id_lists.py`
- `generate_significant_gene_id_lists.sh`
- `go_enrichment.R`
- `summarise_topgo_output.sh`
- `go_barplots.R`

## Stage 07: Figures and Tables

Purpose:

- generate figure/table outputs and supplementary summaries

Scripts/notebooks in `scripts/07_figures_tables/`:

- `omega_plots.py`
- `Hyphy Omega Plots.ipynb`
- `Odds Ratio Test Plots.ipynb`
- `Supplementary Data Tables.ipynb`
- `UpSet Plots.ipynb`

## src Modules

`src/` contains general helper modules, including:

- `id_converter.py`
- `orthogroup_filter.py`
- `orthogroup_gene_count.py`
- `get_silk_genes.py`

## Dependencies

Python dependencies are managed in `environment.yml`.

R dependencies are managed via `renv` (`renv.lock`).

External software commonly used by the full workflow (depending on stage):

- HyPhy
- PREQUAL
- MACSE
- IQ-TREE
- OrthoFinder
- BLAST+
- CD-HIT
- BUSCO

## Results and Supplementary Files

Key outputs are under `results/`, including:

- `results/significant_gene_id_lists/`
- `results/go_enrichment/`
- `results/hyphy_results_cache/`
- supplementary tables

## Terminology

This repository uses **permulation / permulated** terminology for the RERconverge method.

## Change Log

See `CHANGELOG.md` for repository reorganization details.
