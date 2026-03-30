# Odds Ratio Permulation Module

Module path: `scripts/05_permulation_loss_dup/odds_ratio_test.py`

This module runs the odds-ratio permulation workflow used for loss/duplication analyses.

## Main Entry Point

```python
odds_ratio_test(
    test,
    foreground_list_filename="data/orbweavers-list.txt",
    hog_node_genes_tsv="data/N5.tsv",
    genecount_csv="data/N5.GeneCount.tsv",
    occupancy_threshold=0,
    max_occ=None,
    alternative="less",
    alpha=0.05,
    permutation_reps=10000,
    test_triple_gaussian_params=True,
    permulation_tip_values=None,
    permulations_tip_values_csv="data/perms_tip_values.csv",
    background_list_filename=None,
    species_of_interest=None,
    results_dir=None,
    fg_name=None,
    bg_name=None,
    buscos_filename="data/buscos.csv",
    correct_for_buscos=True,
    save_pickle=True,
    save_two_tailed_hits=False,
)
```

## Key Options

- `test`: analysis type, typically `"loss"` or `"duplication"`.
- `alternative`: tail direction, typically `"less"` (loss) or `"greater"` (duplication).
- `occupancy_threshold` and `max_occ`: occupancy filter bounds.
- `permutation_reps`: number of permulation/permutation replicates.
- `permulation_tip_values`: in-memory tip-value assignments (optional).
- `permulations_tip_values_csv`: CSV of precomputed tip-value assignments.
- `species_of_interest`: focal species for annotation/filtering.
- `results_dir`: output root directory for test results.
- `fg_name` and `bg_name`: labels used in output naming.
- `correct_for_buscos`: apply BUSCO correction when available.
- `save_pickle`: write result objects to `.pkl` for later reuse.
- `save_two_tailed_hits`: optionally save two-tailed hit lists.

## Inputs and Large Data Notes

- The module expects local input assets and may read/write large result files.
- Some generated `.pkl` outputs are intentionally not tracked in GitHub due to size.
- Use repository-relative paths where possible to keep runs reproducible.

## Recommended Usage

The notebook `scripts/05_permulation_loss_dup/Odds Ratio Permulation Test.ipynb` demonstrates practical usage and parameter settings for this module.
