# Odds Ratio Permulation Module

Module path: `scripts/04_permulation_loss_dup/odds_ratio_test.py`

This module runs the odds-ratio permulation workflow used for loss/duplication analyses.

## Main Entry Point

```python
odds_ratio_test(
    foreground_list_filename="data/orbweavers-list.txt",
    hog_node_genes_tsv="data/N5.tsv",
    genecount_csv="data/N5.GeneCount.tsv",
    min_occ=0,
    max_occ=None,
    alpha=0.05,
    permulation_reps=10000,
    permulations_tip_values_csv="data/perms_tip_values.csv",
    background_list_filename=None,
    species_of_interest=None,
    results_dir=None,
    fg_name=None,
    bg_name=None,
    buscos_filename="data/buscos.csv",
    correct_for_buscos=True,
    save_pickle=True,
    dir_suffix=None,
)
```

## Key Options

- `min_occ` and `max_occ`: occupancy filter bounds.
- `permulation_reps`: maximum number of permulation replicates to use from the CSV.
- `permulations_tip_values_csv`: CSV of precomputed tip-value assignments.
- `species_of_interest`: focal species for annotation/filtering.
- `results_dir`: output root directory for test results.
- `fg_name` and `bg_name`: labels used in output naming.
- `correct_for_buscos`: apply BUSCO correction when available.
- `save_pickle`: write result objects to `.pkl` for later reuse.
- `dir_suffix`: optional suffix appended to the generated results directory name.

The function returns a `PermulationTestResults` object. When `results_dir` is set,
`fg_name` is required and output files are written to a uniquely named subdirectory.

## Inputs and Large Data Notes

- The module expects local input assets and may read/write large result files.
- Some generated `.pkl` outputs are intentionally not tracked in GitHub due to size.
- Use repository-relative paths where possible to keep runs reproducible.

## Recommended Usage

The notebook `scripts/04_permulation_loss_dup/Odds Ratio Permulation Test.ipynb` demonstrates practical usage and parameter settings for this module.

## Batched log odds ratios

`calculate_odds_batch` in `odds_ratio_test_helpers.py` computes log odds ratios for many permulations at once (a species-by-permulation matrix in,
a HOG-by-permulation matrix out) and matches `calculate_odds` to about 1e-14. `permulated_foreground_matrices` builds those matrices from the permulation tip values.

## Generating more permulations

`permulations_batch.R` generates one batch of permulations (batch id, count and output directory as arguments).
`run_permulation_batches.sh` runs several batches as local parallel processes, and `permulations_batch_slurm.sh` runs them as a SLURM array.
`combine_permulation_batches.py` concatenates the batch CSVs and checks that the foreground count is constant.

## How the p-values are computed

For each HOG and test, the log odds ratio under the true foreground labeling is compared with the log odds ratio under each permulated labeling.
The four per-HOG p-values are `(1 + count) / (permulations + 1)`, where `count` is the number of permulations whose log odds ratio is at least as
extreme as the true one (inclusive of ties, within a 1e-9 tolerance), in each direction: larger (event more likely in the foreground) or smaller (more likely in the background).
The same add-one estimator is used for the global mean tests. A p-value therefore never reaches 0 and cannot go below `1 / (permulations + 1)`.
The two directions of one test are the two tails of the same comparison, so `alpha` (default 0.05) is treated as a two-sided level and each direction is
called at `alpha / 2` (0.025). `results_df` also carries Benjamini-Hochberg q-values (`Q-value ...`) computed over both directions of each test stacked together.

## Pre-specified directions

The hypotheses are directional: loss is more likely in non-orb-weavers (the background) and duplication is more likely in orb-weavers (the foreground).
`directional_q_values` in `odds_ratio_test_helpers.py` computes Benjamini-Hochberg q-values for just those two one-sided tests, each corrected within its own tail over the HOGs
tested for it (loss: occupancy `min_occ`-`max_occ`; duplication: `min_occ` and above). `PermulationTestResults` adds them to `results_df` as
`Q-value loss more likely in bg (directional)` and `Q-value duplication more likely in fg (directional)`, and `results_summary.txt` reports how many HOGs reach q < 0.05.
The two-sided `Q-value ...` columns, which stack both tails of each test, are still produced. For results files written before this existed,
call `directional_q_values(results_df, min_occ, max_occ)` on the `results_all.csv` table.

## Which HOGs get p-values

The permulation loop only evaluates the HOGs that are tested: loss for occupancy `min_occ`-`max_occ` and duplication for occupancy `min_occ` and above.
In `results_all.csv`, a HOG outside a test's range has an empty p-value (and q-value) for that test, for example the loss p-values of HOGs present in more than `max_occ` species.
Log odds ratios are still reported for every HOG.

## Running from the command line

`run_odds_ratio_test.py` wraps `odds_ratio_test()` with the notebook's defaults (use `run_odds_ratio_test_slurm.sh` on a cluster). It uses every row of the permulation CSV unless `--reps` is given.
`--foreground-list` / `--background-list` restrict the analysis to a subset of species (both lists must be given), `--max-occ` sets the upper occupancy bound for loss, and `--no-pickle` skips the large `results.pkl`.

```bash
python scripts/04_permulation_loss_dup/run_odds_ratio_test.py \
    --permulations-csv data/perms_tip_values_1M.csv --dir-suffix 1M_addone --no-pickle
```

## Power analysis

`power_analysis.py` asks what the test could detect if genes really tracked orb-weaving. It generates synthetic genes whose loss (or duplication) follows the true labeling to a chosen degree
("1.0 vs 0.0" follows it perfectly, "0.5 vs 0.5" is a negative control), tests them against every permulation like a real HOG, and mixes them into the real p-values to see whether Benjamini-Hochberg
would call 1, 2, 5 or 20 such genes, with both tails stacked and with the pre-specified tail only. The synthetic genes have independent species-level noise, so the rates are an upper bound on real power.
Results are in `results/power_analysis/`.

## Permulation generation options

`permulations_batch.R` takes `<batch_id> [ntrees] [outdir] [exclude_file] [seed_base]`. `exclude_file` lists extra species to drop from the tree (foreground species that are dropped no longer count),
and `seed_base` (seed = `seed_base + batch_id`) should differ between analyses. `combine_permulation_batches.py` has `--max-batch` and orders batches numerically.
The species-removal analysis uses `data/hub_species_to_exclude.txt`, `data/orbweavers-list_nohubs.txt` and `data/non-orb-weavers-list_nohubs.txt`.

## Presentation figures

The last section of `scripts/07_figures_tables/Odds Ratio Test Plots.ipynb` draws figures explaining what a permulation is, why it differs from a random shuffle, and where the candidate HOG falls in its null distribution. They are saved to `figures/presentation/`.

## Cell-level flags (suspected contamination)

A species-by-HOG cell can be flagged when that species' member of the HOG looks like a contaminant. `--flagged-cells FILE` (`flagged_cells_file` in `odds_ratio_test()`) treats flagged cells as missing data in that HOG only; the species stays in every other HOG and the permulations are unchanged (same tree, same foreground count).

- The table is a CSV/TSV with `HOG` and `species` columns plus either `strength` (gene-level, one row per flagged gene; `--flag-strengths strong` or `--flag-strengths strong moderate` chooses which rows count) or `n_flagged` (cell-level). The column names are provisional until the final flag table exists (`load_flagged_cells` in `odds_ratio_test_helpers.py`).
- Flagged copies are subtracted from the gene counts, which can change a duplication call (count > 1). A species left with no copy is missing, not absent: it counts as neither a loss nor a duplication, and it is left out of that HOG's occupancy.
- Missing cells get weight 0 in that HOG's test. `calculate_odds_batch_weighted` takes a HOG-by-species weight matrix (BUSCO weight times 0 or 1) in place of the single per-species BUSCO vector, and is used only for HOGs with at least one missing cell, for the true values and in every permulation. All other HOGs use the unchanged code, so their results are identical.
- The log shows how many HOGs were touched, how many changed occupancy, and how many HOGs are tested for loss and duplication before and after. `results_all.csv` gets `Missing flagged cells` and `Occupancy (before flagging)` columns, and `results_summary.txt` reports the flagged cells.
- Checks: with an empty table the log odds ratios equal the existing 100k `results_all.csv` to 1e-15; `calculate_odds_batch_weighted` equals `calculate_odds_batch` with BUSCO weights; a hand-computed toy example matches; masking the 8 non-orb-weaver carriers of N5.HOG0071799 gives a loss log odds ratio of -6.03 (-3.26 as analysed) and duplication 3.49 (3.64).

