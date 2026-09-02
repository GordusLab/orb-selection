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
