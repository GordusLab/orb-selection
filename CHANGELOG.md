# Changelog

## 2026-03-30 - Follow-up Cleanup and Methods Alignment

### Changed

- Refined top-level documentation in `README.md` to emphasize methods-aligned stage usage and external large-data notes.
- Added/retained explicit placeholders in `README.md` for future raw-data, HyPhy JSON archive, and PKL archive locations.
- Flattened figure/table notebook layout by moving notebooks from `scripts/07_figures_tables/notebooks/` to `scripts/07_figures_tables/`.
- Updated `.gitignore` to stop ignoring all of `data/` and instead ignore only:
  - `data/absrel_rerun/`
  - `data/busted_ph/`
  - `data/busted_ph_rev/`
  - `data/relax/`
- Kept `*.pkl` ignored repository-wide.

### Removed

- Removed `scripts/runners/` and all stage runner wrapper scripts by request.

### Fixed

- Repaired indentation/scope issues in `scripts/03_phylo_prep/get_nuc_seqs.py`.
- Made subprocess checks explicit in `scripts/03_phylo_prep/get_nuc_seqs.py` and removed an unused import.

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
  - `scripts/03_phylo_prep/`
  - `scripts/04_selection_tests/`
  - `scripts/05_permulation_loss_dup/`
  - `scripts/06_enrichment/`
  - `scripts/07_figures_tables/`
- Moved figure/table notebooks into stage 07 (subsequently flattened to `scripts/07_figures_tables/`).
- Moved `Odds Ratio Permulation Test.ipynb` and `Multimodal Test.ipynb` to `scripts/05_permulation_loss_dup/`.
- Updated notebook setup/path cells for stage-based execution.
- Added workflow/context markdown cells to moved notebooks.

### Path and Portability Updates

- Replaced hardcoded local absolute paths in moved shell scripts with environment-variable placeholders where practical.
- Updated scripts to be repo-relative where possible.
- Updated `scripts/03_phylo_prep/get_nuc_seqs.sh` to call colocated `get_nuc_seqs.py`.
- Moved `get_silk_genes.py` to `src/get_silk_genes.py` and updated notebook/module references.
- Generalized `scripts/06_enrichment/summarise_topgo_output.sh` to summarize across enrichment directories.

### Selection and Permulation Utilities

- Updated gene ID list generation script location and references:
  - `scripts/06_enrichment/get_gene_id_lists.py`
  - `scripts/06_enrichment/generate_significant_gene_id_lists.sh`
- Moved `odds_ratio_test.py` to `scripts/05_permulation_loss_dup/` and updated imports.
- Moved `hyphy_results_parser.py`, `hyphy_results_helpers.py`, and parser documentation to `scripts/04_selection_tests/`.
- Moved parser workflow script from tests to `scripts/04_selection_tests/parse_hyphy_results.py`.
- Moved `omega_plots.py` to `scripts/07_figures_tables/`.
- Updated wording and helper naming toward permulation terminology.
- Modernized result-loading checks in `get_gene_id_lists.py` to align with current cached result classes.

### Fixes

- Fixed incorrect import in `scripts/04_selection_tests/hyphy_results_helpers.py` causing false warning noise:
  - import now uses `convert_hogs_to_locs` from `id_converter`.

### Removed

- Temporary comparison helper scripts created during validation (removed after verification by request).

### Notes

- This pass focused on structure and execution safety, not scientific parameter changes.
- Some analyses referenced in manuscript methods remain external/manual and should use default settings unless otherwise specified.
