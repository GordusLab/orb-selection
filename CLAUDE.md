# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Analysis repository for the paper "Comparative transcriptomic analysis reveals signatures of
selection for orb-weaving behavior in spiders" (Runnels et al. 2026). It is a bioinformatics
pipeline, not an application: a `scripts/` directory of numbered stages (01-07) that take raw
transcriptome data through orthology inference, molecular-evolution selection tests, gene
loss/duplication analysis, phylogenetic regression, GO enrichment, and figure/table generation.
Most heavy stages (02, 03) are designed to run on a SLURM HPC cluster (JHU "agordus1" account);
later stages (04-07) are mostly local R/Python scripts and Jupyter notebooks.

Read [README.md](README.md) first — it documents each pipeline stage, its steps, and its inputs/outputs
in detail. [CHANGELOG.md](CHANGELOG.md) records repo reorganizations. Several `scripts/*/README.md` files
document stage-specific module APIs (see below).

## Environments

There is no single environment — different stages use different, separately-managed environments:

- **R** (`renv.lock`, `.Rprofile`): used for stage 04 permulations, stage 05 phyloGLM, stage 06
  topGO enrichment. Requires R 4.6.1 / Bioconductor 3.23.
  ```bash
  R --vanilla -s -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
  Rscript -e 'renv::restore(prompt = FALSE)'
  ```
  `RERconverge` must be installed manually (not via renv) — see its wiki. On macOS, `gdtools`
  needs XQuartz/X11.
- **Python** (`environment.yml` / `environment_macos.yml` conda envs): general pipeline Python
  (biopython, pandas, numpy, matplotlib, scipy). Use the macOS variant on Apple hardware — the
  main `environment.yml` is a Linux snapshot and won't solve on macOS.
- **`selection.yml`**: conda env for running HyPhy selection tests (stage 03) on the cluster.
- **`upset.yml`**: conda env just for the UpSet Plots notebook (stage 07).
- `.venv/` at the repo root is an ad hoc local Python 3.13 venv (not one of the above); don't
  assume it has pipeline dependencies installed.

There is no linter or test suite configured for this repo — don't look for `pytest`/lint commands.

## Pipeline architecture

Each stage directory under `scripts/` corresponds to a methods-section stage and is meant to be
read in the order described in the README. The high-level data flow:

1. **01_pre_processing** — download transcriptomes, cluster (CD-HIT), ORF-call (TransDecoder),
   homology filter (BLASTP), quality-check (BUSCO). Pure shell scripts with placeholder paths.
2. **02_orthofinder_prep_hyphy** — OrthoFinder orthology inference (run manually/externally; the
   exact invocations across several re-runs are logged in
   [scripts/02_orthofinder_prep_hyphy/README.md](scripts/02_orthofinder_prep_hyphy/README.md)),
   then a SLURM array job (`prep_for_hyphy2.sh`, one task per orthogroup) runs
   PREQUAL → MACSE → sequence filtering (`filter_msa_sequences.py`) → ClipKIT trimming →
   remove-duplicates → IQ-TREE → foreground branch labeling (`label-tree.bf`), producing one
   alignment + two labeled trees (orb-weaver fg, non-orb-weaver fg) per HOG.
3. **03_selection_tests** — a second SLURM array job (`run_selection_tests.sh`, 3 tasks per HOG)
   runs BUSTED-PH (orb fg), BUSTED-PH (non-orb fg), and RELAX per orthogroup; `submit_tiers.sh`
   splits HOGs into runtime tiers (A/B/C) from prior run-time records so job time/CPU/memory
   requests are sized appropriately. `hyphy_results_parser.py` / `hyphy_results_helpers.py`
   (documented in [scripts/03_selection_tests/README.md](scripts/03_selection_tests/README.md))
   load the resulting JSONs into `HyphyResultsManager`/`RelaxResult`/`BustedPhResult` objects and
   cache them as `.pkl` under `results/`.
4. **04_permulation_loss_dup** — `permulations.R` generates permulated phenotype tip assignments;
   `odds_ratio_test.py` (API documented in
   [scripts/04_permulation_loss_dup/README.md](scripts/04_permulation_loss_dup/README.md)) runs the
   gene loss/duplication odds-ratio permulation test using those assignments plus BUSCO-based
   correction.
5. **05_phyloglm** — `phyloglm.R` fits a phylogenetic GLM per gene (~12,000 genes, run in
   parallel); results inspected in the paired notebook.
6. **06_enrichment** — BLASTs orthogroups against *P. tepidariorum* / *U. diversus* reference
   genomes, builds significant-gene-ID lists from stage 03/04 outputs, and runs topGO enrichment
   (`go_enrichment.R`), summarized across result sets by `summarise_topgo_output.sh`.
7. **07_figures_tables** — notebooks/scripts turn stage 03-06 outputs into manuscript figures,
   UpSet intersections, and the supplementary data tables; expects completed upstream stages.

`src/` holds cross-stage helper modules, notably:
- `id_converter.py` — converts between *U. diversus* transcript IDs, LOC gene IDs, and HOGs, and
  attaches gene descriptions/*D. melanogaster* orthologs; used throughout stages 03-07 whenever
  HOG- or transcript-level results need to be annotated with gene identity.
- `orthogroup_filter.py` — occupancy / single-copy / taxon-presence filtering of OrthoFinder
  gene-count tables (used to derive the stage-02 input HOG lists).

## Key data-flow objects (`data/`)

See [data/README.md](data/README.md) for the full file-by-file description. The ones most scripts key off:
- `data/N5.tsv` / `data/N5.GeneCount.tsv` — the OrthoFinder N5-node hierarchical orthogroups and
  per-species gene counts; `data/N5.udiv.o75_list.txt` is the resulting list of 4,756 HOGs
  (occupancy ≥ 75, *U. diversus* present) actually run through the stage 02/03 HyPhy pipeline.
- `data/orbweavers-list.txt` / `data/non-orb-weavers-list.txt` — the phenotype grouping (foreground
  labeling for selection tests, foreground/background for the odds-ratio test) used throughout
  stages 02-06.
- `data/species_list_N5.txt` (98 spp. analysis set) vs `data/species_list_all.txt` (102 spp.,
  OrthoFinder input incl. outgroup) — don't conflate these two species sets.

Large/derived inputs (raw transcriptome FASTAs, per-gene HyPhy JSON outputs, some cache `.pkl`
files) are external and not in this repo — see the README's "Data Availability" section and
`.gitignore` (which excludes `*.pkl`, several `data/<test>/` result subdirs, and `/inspect/`).

## Conventions in the SLURM/shell scripts

The stage 02/03 array scripts (`prep_for_hyphy2.sh`, `run_selection_tests.sh`) follow a consistent
resumability pattern worth preserving when editing them:
- Each pipeline step checks whether its output file already exists/is non-empty (sometimes
  grepping the tool's own log for a completion marker, since some tools can exit 0 without
  writing output) before re-running, so a failed/requeued array task resumes rather than restarts.
- `RERUN_FROM_TRIM=1` / `STOP_AFTER_TRIM=1` env vars control re-running from a given step; reruns
  archive (move aside, never delete) prior outputs under a timestamped `old_trim_*` directory
  rather than overwriting them.
- Paths to conda envs, work dirs, JAR files, etc. are left as placeholders (`PATH_TO_...`) or
  shell positional args in the versions meant to be portable/example — don't hardcode the
  original author's cluster paths (`/scratch4/agordus1/crunnel2/...`, `$HOME/orb-selection`) into
  new reusable scripts.
- `scripts/03_selection_tests/submit_tiers.sh` encodes HOG runtime tiers (A/B/C, from prior Pass-2
  run-time records) to right-size `--time`/`--cpus-per-task` per job batch; `requeue_oom.sh`
  handles re-submitting OOM-killed tasks at higher memory.
