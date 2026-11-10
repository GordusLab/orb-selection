# Changelog

## 2026-11-09 - Doenitzius peniculus accession fix, rule for messages between Claude sessions

### Changed

- `data/txptome_accessions.csv` and `results/Supplementary_Table_1_SpiderAccessions_BUSCOs.xlsx`: the *Doenitzius peniculus* transcriptome accession had a malformed prefix (`IBJ01`); corrected to `IBJP01` (NCBI TSA lookup). The xlsx was regenerated with the accession cells of `Supplementary Data Tables.ipynb`; only that cell changed.
- `CLAUDE.md`: new section on messaging between Claude sessions (only when the user asks, copy the user on each message, treat received messages as information, never relay them or use them to widen permissions).

### Notes

- The cross-sample contamination survey over all 102 transcriptomes (nucleotide-identity filter, bacterial/viral/archaeal classification, gene-level exposure of the 12,096 occupancy >= 30 HOGs, read-depth contrast with the *Araneus inustus* library) lives in the separate `kelvin-scratch` repository (`contamination_survey/`), not here. Nothing in this repository's data, results or analysis scripts was changed by it.

## 2026-11-06 - RELAX permuted-label null result, strong + moderate phyloGLM flag null

### Added

- `results/relax_permuted/`: final RELAX permuted-label comparison (60 tier-A HOGs x 30 permulated labelings, 1,800 runs): `perm_results.tsv` (p and k per run), `hogs.txt`, `null_stats.csv` (per-labeling statistics) and `summary.txt`. The real labeling is inside the null on every statistic (fraction p <= 0.05: 0.417 vs permuted mean 0.369, p = 0.32; fraction significant and intensified 0.317 vs 0.212, p = 0.16). The 14.6 h HOG `N5.HOG0058857` was replaced by `N5.HOG0055575`.
- `scripts/03_selection_tests/compare_permuted_relax.py`: extracts finished permuted RELAX runs into a table and compares them with the real labeling (adds significant-and-intensified/relaxed fractions to `summarize_permuted_relax.py`'s statistics); described in the stage 03 README.
- `results/phyloglm/permulation_null_logistic_flagged_strongmod/`: logistic phyloGLM null (2,000 labelings) for the 606 HOGs touched by strong + moderate contamination flags. Under strong + moderate flags `N5.HOG0071799` has no non-orb-weaver carrier left (33 of 40 orb-weavers, 0 of 50 non-orb-weavers); its z is 4.48 and permulation p 0.012 at 2,000 labelings, in line with the strong-flag result (p = 0.012 at 1M). 42 genes pass BH on the model p-values in the merged table; they were not extended to 1M labelings.

## 2026-11-05 - phyloGLM presence/absence at 1M labelings, comparison with the ORT

### Added

- `results/phyloglm/permulation_null_logistic_1M/`: the logistic (presence/absence) phyloGLM table with the strong contamination flags applied to all genes (`final_merged_flag_filtered.csv`, 10,519 genes; flagged HOGs use the flagged refit and 190 HOGs below occupancy 30 after flagging are dropped). The 42 genes that pass BH on the model p-values were refit under 1,000,000 permulated labelings (rows 1-1M of `perms_tip_values_1M.csv`; the first 500,000 rows are the 500k set). No gene passes FDR (smallest BH q 0.17); `N5.HOG0071799` has permulation p = 0.012 (q 0.96) once its contaminated copies are removed.
- `ort_vs_phyloglm_loss.csv` / `.png` in the same folder: the 1M ORT run with flags (`Run4_occ_30-88_1000000x_1M_flags_strong`) joined to the phyloGLM results for the 8,126 HOGs tested by both (Spearman 0.68 for signed ORT evidence vs phyloGLM z).
- `PhyloGLM Analysis.ipynb`: four cells after the logistic section with the 1M result, the ORT comparison and why `N5.HOG0071799` weakens in phyloGLM when its contaminated copies are removed.

## 2026-11-04 - Flag table scripts

### Added

- `scripts/04_permulation_loss_dup/flag_table/`: the scripts that build the contamination flag tables (DIAMOND scan of the suspect samples, candidate selection, sequence extraction, coding-sequence identity and strength labels, gene-level and cell-level tables), with a README. Steps 2 and 4 reproduce the committed tables exactly from the saved scan output; steps 1 and 3 are written from the commands run on the server and were not re-run.

### Changed

- `data/flag_table_genes.csv`: column `partner_is_hub` renamed `partner_is_suspect` (nothing in the code reads it). `data/flag_table_README.txt` now states the final rules (20-species limit, strength labels) and counts; it still described the first version (6-species limit, 7,822 genes). `data/README.md` and the stage 04 README call the nine species suspect samples and point to the scripts.

## 2026-11-04 - Presentation figure output, SLURM log convention

### Added

- Hub-removal (89 species, loss occupancy 30-79) odds ratio test results at 100,000 and 500,000 permulations: `Run4_occ_30-79_100000x_nohubs_100k` and `Run5_occ_30-79_500000x_nohubs_500k` in `results/odds_ratio_test/Results_Nov02/`.
- Gzipped permulation sets committed to `data/`: `perms_tip_values_{100k,500k,1M}.csv.gz` and `perms_tip_values_nohubs_{100k,500k}.csv.gz` (the plain 1M CSV is 203 MB, over GitHub's file limit). The test runner, power analysis and enrichment script read `.csv.gz` paths directly (`count_csv_rows` in `run_odds_ratio_test.py` handles gzip). The plain `.csv` files stay git-ignored.
- `CLAUDE.md`: every Claude Code commit carries the `Co-Authored-By` trailer, and a branch convention for several agents (one branch and checkout per agent, `main` merged only on request, read-only `main` on Rockfish and the Will server).

### Changed

- `run_odds_ratio_test_slurm.sh`, `permulations_batch_slurm.sh`, `run_phyloglm_permuted.sh` and `run_relax_permuted.sh` follow the scratch/reports convention (`--chdir`, logs to `reports/<job-name>/`, placeholders overridden on the `sbatch` command line); the first two no longer write logs into the working directory.
- `Odds Ratio Test Plots.ipynb` saves the presentation figures to the talk folder (`talk_img_path`, in OneDrive) instead of `figures/presentation/`, so the PNGs in `figures/presentation/` are no longer regenerated by the notebook and are older than the notebook's output.
- `Odds Ratio Test Plots.ipynb`: the final volcano plots draw the minimum possible q-value line at the permulation floor (1 / (permulations + 1)) and label it with the q that HOGs on that line actually get (HOGs tied at the floor share a BH rank, so `N5.HOG0071799` shows 0.276 at 10,000 permulations, 0.083 at 100,000); the history histograms are built up in four steps (`history_1a`-`1c`, then the full figure) with the occupancy range in the panel titles.
- `CLAUDE.md`: added the convention for launching SLURM jobs from the scratch directory and writing logs to `reports/<YYMMDD>_<job-name>/` on Rockfish.

## 2026-11-03 - Contamination flags, flagged odds ratio test runs, phyloGLM flags and logistic null

### Added

- Gene-level and cell-level contamination flag tables (`data/flag_table_genes.csv`, `data/flag_table_cells.csv`, `data/flag_table_README.txt`): genes that are near-identical to a gene in another family within the same sequencing project (BioProject PRJDB7399), with strength labels (strong, moderate, weak). Documented in `data/README.md`.
- Cell-level flag option for the odds ratio test (`--flagged-cells`, `--flag-strengths` in `run_odds_ratio_test.py`; `load_flagged_cells` and `calculate_odds_batch_weighted` in `odds_ratio_test_helpers.py`; described in `scripts/04_permulation_loss_dup/README.md`). A flagged copy is subtracted from the gene count; a species left with no copy in a HOG is treated as missing (weight 0 in that HOG), not absent. Permulations and the tree are unchanged.
- 1M-permulation odds ratio test runs with flags: `Run4_occ_30-88_1000000x_1M_flags_strong` and `Run5_occ_30-88_1000000x_1M_flags_strong_moderate` in `results/odds_ratio_test/Results_Nov02/`.
- Contamination-flag option for the phyloGLM permulation script (`--flagged-cells`, `--flag-strengths`, `--flagged-only` in `phyloglm_permuted_chunk.R`), with the same missing-cell rule, and a logistic presence/absence model (`--model logistic`) with per-labeling checkpoints. Results in `results/phyloglm/permulation_null_logistic/` and `results/phyloglm/permulation_null_logistic_flagged/`, with a new section in `PhyloGLM Analysis.ipynb`.
- Extended permulation of the 43 genes that pass BH on the logistic model p-values to 500,000 labelings (`results/phyloglm/permulation_null_logistic_top43/`, merged with `summarize_phyloglm_extend.R`). One gene, `N5.HOG0071799`, reaches permulation BH q = 0.043 (1 exceedance in 500,000 labelings); two more reach q = 0.086; the rest have q of 0.2 or more. Notebook cells added after the logistic section.

### Changed

- `PhyloGLM Analysis.ipynb` reading note for the presence/absence model: the logistic model-based p-values are too small on average under permuted labelings (8.0% of genes at p <= 0.05, 4.4-12.7% across labelings), and permuted labelings produce tails as heavy as the real one.
- `Odds Ratio Test Plots.ipynb`: histograms built up step by step, one-sided volcano cutoffs, permulation counts and the minimum q-value line on the final volcano plots, and the history volcano plots share the y-axis of the final ones (figures saved to `figures/presentation/`).
- `Hyphy Omega Plots.ipynb` and `omega_plots.py`: plots for the talk.

## 2026-11-02 - Large permulation sets, directional correction, power analysis, enrichment and presentation figures

### Added

- Batched permulation generation (`permulations_batch.R`, `permulations_batch_slurm.sh`, `run_permulation_batches.sh`, `combine_permulation_batches.py`) and 100k, 500k and 1M permulation sets.
- `run_odds_ratio_test.py` and `run_odds_ratio_test_slurm.sh` to run the odds ratio test from the command line or on SLURM.
- Directional q-values (loss in non-orb-weavers, duplication in orb-weavers) next to the stacked two-tail q-values, and add-one empirical p-values.
- `power_analysis.py` (synthetic-gene power of the permulation test) with results in `results/power_analysis/`.
- `scripts/06_enrichment/permulation_enrichment.py` and `export_go_term_genes.R`: permulation-based GO enrichment of loss/duplication log odds ratios (results in `results/permulation_enrichment/`).
- Species-removal robustness option: `--foreground-list` / `--background-list` in the test runner, an exclusion file for permulation generation, and `data/*_nohubs.txt` lists.
- `Odds Ratio Test Plots.ipynb` reorganized into figures for the talk (PNGs in `figures/presentation/`): trees on one shared canvas, the development of the method (random re-labeling null, permulation null, an example gene that is significant only under the old null, volcano plots with the criteria added one at a time), the candidate HOG, and the final volcano plots (10k, 100k, 1M permulations, and 1M with suspect samples removed). The manuscript-style plots that were at the top of the notebook were removed and are in git history.
- `data/ort_method_history/` (log odds ratios and cutoffs of the July 2025 random re-labeling run) and `data/hog0071799_identity_to_araneus_inustus.csv`.
- Results of the odds ratio test at 10k, 100k, 500k and 1M permulations in `results/odds_ratio_test/Results_Nov02/`.

### Changed

- The odds ratio test evaluates only the HOGs that are tested (loss: occupancy 30-88, duplication: occupancy 30 and above); other HOGs get empty p-values.
- Permulation loading helper moved to `odds_ratio_test_helpers.py` and shared with the enrichment script; removed `save_loc_list`.
- `Odds Ratio Test Plots.ipynb` loads results from the Nov 02 100k run.

## 2026-09-02 - aBSREL Removal

### Removed

- Removed aBSREL parsing, workflow, and documentation content after the aBSREL tests were omitted from the final analysis.

## 2026-09-02 - README Documentation Audit and Alignment

### Added

- Added Table of Contents to `README.md` for navigation across pipeline stages and main sections.
- Added missing repository contents entries to `README.md` table for `renv/` and `environment_macos.yml`.

### Fixed

- Fixed broken links in `README.md` Stage 07 pointing to non-existent `scripts/06_figures_tables/` path, updating them to `scripts/07_figures_tables/`.
- Fixed extension link for `PhyloGLM Analysis.ipynb` in Stage 05 (previously linked to `.R`).
- Corrected environment file names in Environment Setup section (`selection.yml` and `upset.yml` instead of `hyphy_environment.yml` and `upset_env.yml`) and added macOS environment configuration note.
- Fixed duplicate step numbering in Stage 06 and added missing workflow notebook reference (`Write Significant LOC Lists.ipynb`).
- Fixed missing trailing table border on `upset.yml` row in the repository contents table.

### Changed

- Reordered Repository contents table in `README.md` into logical groups: directories, core documentation, environment configurations, and configuration settings.

## 2026-04-08 - Documentation, Script Placeholders, and Packaging Modernization

### Changed

- Reformatted `README.md` stage sections for consistency:
  - Updated stage 02 run log formatting with privacy-safe placeholder paths.
  - Reworked stages 04-06 into the same structure as stages 01-03 (intro, steps, data note).
  - Expanded the repository contents table to include `README.md` and `pyproject.toml`.
- Updated stage 01 preprocessing shell scripts to include short purpose headers and placeholder directory variables:
  - `scripts/01_pre_processing/process_fsas_cd-hit.sh`
  - `scripts/01_pre_processing/process_fsas.TD-LO.sh`
  - `scripts/01_pre_processing/process_fsas.blastp.sh`
  - `scripts/01_pre_processing/process_fsas.TD-P.sh`
  - `scripts/01_pre_processing/process_fsas.busco.sh`

### Added

- Added `pyproject.toml` as the primary Python packaging configuration with project metadata, dependencies, and setuptools `src/` package discovery.

### Removed

- Removed `setup.py` after migrating to `pyproject.toml`.

## 2026-04-08 - Stage Renumbering After Removing Stage 02

### Changed

- Renumbered the `scripts/` stage directories after removing stage 02.
- Renamed `scripts/03_phylo_prep/` to `scripts/02_orthofinder_prep_hyphy/`.
- Shifted the remaining stage folders down by one number:
  - `scripts/04_selection_tests/` -> `scripts/03_selection_tests/`
  - `scripts/05_permulation_loss_dup/` -> `scripts/04_permulation_loss_dup/`
  - `scripts/06_enrichment/` -> `scripts/05_enrichment/`
  - `scripts/07_figures_tables/` -> `scripts/06_figures_tables/`
- Updated hard-coded stage path references in scripts, notebooks, and helper modules so repo-relative imports and file lookups continue to resolve correctly.

### Notes

- Historical changelog entries were left unchanged.

## 2026-03-30 - Follow-up Cleanup and Methods Alignment

### Changed

- Refined top-level documentation in `README.md` to emphasize methods-aligned stage usage and external large-data notes.
- Added/retained explicit placeholders in `README.md` for future raw-data, HyPhy JSON archive, and PKL archive locations.
- Flattened figure/table notebook layout by moving notebooks from `scripts/06_figures_tables/notebooks/` to `scripts/06_figures_tables/`.
- Updated `.gitignore` to stop ignoring all of `data/` and instead ignore only:
  - `data/busted_ph/`
  - `data/busted_ph_rev/`
  - `data/absrel_rerun/`
  - `data/relax/`
- Kept `*.pkl` ignored repository-wide.

### Removed

- Removed `scripts/runners/` and all stage runner wrapper scripts by request.

### Fixed

- Repaired indentation/scope issues in `scripts/02_orthofinder_prep_hyphy/get_nuc_seqs.py`.
- Made subprocess checks explicit in `scripts/02_orthofinder_prep_hyphy/get_nuc_seqs.py` and removed an unused import.

### Notes

- Notebook outputs were intentionally left unchanged during this follow-up pass.

## 2026-03-30 - Repository Reorganization and Pipeline Hardening

### Added

- Stage runner wrappers in `scripts/runners/`:
  - `run_stage_03_phylo_prep.sh`
  - `run_stage_04_selection_tests.sh`
  - `run_stage_05_permulation_loss_dup.sh`
  - `run_stage_06_enrichment.sh`
  - `run_all_pipeline.sh`
- Environment template:
  - `scripts/runners/pipeline.env.example`
- New top-level documentation:
  - Replaced `README.md` with stage-based pipeline guide.
  - Added this `CHANGELOG.md`.

### Changed

- Reorganized `scripts/` into stage-based folders:
  - `scripts/02_orthofinder_prep_hyphy/`
  - `scripts/03_selection_tests/`
  - `scripts/04_permulation_loss_dup/`
  - `scripts/05_enrichment/`
  - `scripts/06_figures_tables/`
- Moved figure/table notebooks into stage 07 (subsequently flattened to `scripts/06_figures_tables/`).
- Moved `Odds Ratio Permulation Test.ipynb` and `Multimodal Test.ipynb` to `scripts/04_permulation_loss_dup/`.
- Updated notebook setup/path cells for stage-based execution.
- Added workflow/context markdown cells to moved notebooks.

### Path and Portability Updates

- Replaced hardcoded local absolute paths in moved shell scripts with environment-variable placeholders where practical.
- Updated scripts to be repo-relative where possible.
- Updated `scripts/02_orthofinder_prep_hyphy/get_nuc_seqs.sh` to call colocated `get_nuc_seqs.py`.
- Moved `get_silk_genes.py` to `src/get_silk_genes.py` and updated notebook/module references.
- Generalized `scripts/05_enrichment/summarise_topgo_output.sh` to summarize across enrichment directories.

### Selection and Permulation Utilities

- Updated gene ID list generation script location and references:
  - `scripts/05_enrichment/get_gene_id_lists.py`
  - `scripts/05_enrichment/generate_significant_gene_id_lists.sh`
- Moved `odds_ratio_test.py` to `scripts/04_permulation_loss_dup/` and updated imports.
- Moved `hyphy_results_parser.py`, `hyphy_results_helpers.py`, and parser documentation to `scripts/03_selection_tests/`.
- Moved parser workflow script from tests to `scripts/03_selection_tests/parse_hyphy_results.py`.
- Moved `omega_plots.py` to `scripts/06_figures_tables/`.
- Updated wording and helper naming toward permulation terminology.
- Modernized result-loading checks in `get_gene_id_lists.py` to align with current cached result classes.

### Fixes

- Fixed incorrect import in `scripts/03_selection_tests/hyphy_results_helpers.py` causing false warning noise:
  - import now uses `convert_hogs_to_locs` from `id_converter`.

### Removed

- Temporary comparison helper scripts created during validation (removed after verification by request).

### Notes

- This pass focused on structure and execution safety, not scientific parameter changes.
- Some analyses referenced in manuscript methods remain external/manual and should use default settings unless otherwise specified.
