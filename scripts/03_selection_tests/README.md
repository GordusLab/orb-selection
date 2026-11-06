# HyPhy Results Analysis Module

This module provides classes to store and manage results from HyPhy selection analyses (RELAX and BUSTED-PH).

## Files

- `scripts/03_selection_tests/hyphy_results_parser.py` - Main module with result classes
- `scripts/03_selection_tests/hyphy_results_helpers.py` - Helper functions
- `scripts/03_selection_tests/parse_hyphy_results.py` - Project parser workflow script
- `scripts/03_selection_tests/README.md` - This documentation

## Quick Start

Run the following examples from `scripts/03_selection_tests/`.

```python
from hyphy_results_parser import HyphyResultsManager
from hyphy_results_helpers import convert_hyphy_results_to_locs

# Initialize manager
manager = HyphyResultsManager()

# Load results from JSON directories
relax_result = manager.load_relax_from_json('/path/to/relax/jsons/')

# Convert HOGs to LOCs using orb-selection's id_converter
significant_results = relax_result.get_significant_results()
results_with_locs = convert_hyphy_results_to_locs(significant_results)

# Save in orb-selection results directory
manager.save_all_results('../../results/hyphy_results_cache/')
```

## Integration Features

### 1. HOG to LOC Conversion
Uses `src/id_converter.py`

```python
from hyphy_results_helpers import convert_hyphy_results_to_locs

# Convert any HyPhy results DataFrame
results_with_locs = convert_hyphy_results_to_locs(hyphy_results)
```

### 2. Data Access
Automatically accesses data files from `orb-selection/data/`:
- `data/N5.tsv` - Default hierarchical orthogroup file
- `data/id_converter.tsv` - Gene ID conversion data
- `data/Uloborus_diversus__v__Drosophila_melanogaster.tsv` - Ortholog mapping

### 3. Results Storage
Results are saved in the `orb-selection/results/` directory:
```python
# Cache results for future use
manager.save_all_results('../../results/hyphy_results_cache/')

# Results are saved as:
# results/hyphy_results_cache/relax_results.pkl
# results/hyphy_results_cache/busted_ph_orb_results.pkl
# results/hyphy_results_cache/busted_ph_non_orb_results.pkl
```

## Classes

### HyphyResultsManager
Central manager for handling multiple analysis results.

**Key Methods:**
- `load_relax_from_json(directory)` - Load RELAX results
- `load_busted_ph_from_json(directory)` - Load BUSTED-PH results  
- `compare_significant_genes(analyses)` - Compare gene sets
- `get_overlap_stats(analyses)` - Get overlap statistics
- `save_all_results(directory)` - Save all results to cache

### RelaxResult, BustedPhResult
Individual result classes with analysis-specific methods:

```python
# RELAX-specific
relax_result = RelaxResult(relax_df)
relaxed_genes = relax_result.get_relaxed_results()
intensified_genes = relax_result.get_intensified_results()
selection_counts = relax_result.count_selection_types()

# BUSTED-PH-specific  
busted_result = BustedPhResult(busted_df)
hits = busted_result.get_hits()

```

## Common Operations

### Loading and Converting to LOCs
```python
from hyphy_results_parser import HyphyResultsManager
from hyphy_results_helpers import convert_hyphy_results_to_locs

# Load results
manager = HyphyResultsManager()
relax_result = manager.load_relax_from_json('/path/to/relax/jsons/')

# Get significant results and convert to LOCs
significant = relax_result.get_significant_results()
results_with_locs = convert_hyphy_results_to_locs(significant)

# Now you have LOC IDs and gene descriptions
print(results_with_locs[['LOC', 'Description', 'result', 'p_value']])
```

### Filtering and Analysis
```python
from hyphy_results_helpers import filter_omega, omega_filter_summary

# Filter by omega values (removes outliers)
filtered = relax_result.filter_omega(10000)

# Show filtering summary
omega_filter_summary(relax_result.results_df, [10, 100, 1000, 10000], 'result', 'relaxed')

# Get specific result types
relaxed = relax_result.get_relaxed_results()
intensified = relax_result.get_intensified_results()
```

### Cross-Analysis Comparisons
```python
# Load multiple analyses
manager.load_relax_from_json('/path/to/relax/')
manager.load_busted_ph_from_json('/path/to/busted_ph_orb/', name='busted_ph_orb')
manager.load_busted_ph_from_json('/path/to/busted_ph_non_orb/', name='busted_ph_non_orb')

# Compare significant genes
overlap_stats = manager.get_overlap_stats(['relax', 'busted_ph_orb', 'busted_ph_non_orb'])
gene_sets = manager.compare_significant_genes(['relax', 'busted_ph_orb', 'busted_ph_non_orb'])

print(f"RELAX hits: {len(gene_sets['relax'])}")
print(f"BUSTED-PH orb hits: {len(gene_sets['busted_ph_orb'])}")
print(f"BUSTED-PH non-orb hits: {len(gene_sets['busted_ph_non_orb'])}")
```

## Example Workflow

```python
#!/usr/bin/env python3

from hyphy_results_parser import HyphyResultsManager
from hyphy_results_helpers import convert_hyphy_results_to_locs, filter_omega

# 1. Initial setup (run once)
manager = HyphyResultsManager()
relax_result = manager.load_relax_from_json('/path/to/relax/jsons/')
busted_ph_orb_result = manager.load_busted_ph_from_json(
	'/path/to/busted_ph_orb/jsons/', name='busted_ph_orb'
)
busted_ph_non_orb_result = manager.load_busted_ph_from_json(
	'/path/to/busted_ph_non_orb/jsons/', name='busted_ph_non_orb'
)

# 2. Cache results for future use
manager.save_all_results('../../results/hyphy_results_cache/')

# 3. Analysis sessions (load from cache)
manager = HyphyResultsManager()
manager.load_all_results_from_directory('../../results/hyphy_results_cache/')
relax_result = manager.get_result('relax')

# 4. Filter and analyze
significant = relax_result.get_significant_results()
filtered = relax_result.filter_omega(10000)

# 5. Convert to LOCs using orb-selection id_converter
results_with_locs = convert_hyphy_results_to_locs(filtered.results_df)

# 6. Save final results
results_with_locs.to_csv('../../results/relax_significant_with_locs.csv', index=False)
```

## Permuted-label RELAX null

Calibrates the real RELAX results against permulated phenotype labelings (the same
`data/perms_tip_values_100k.csv` used for the odds ratio test). A fixed random subset of
HOGs is re-run under many permuted foreground labelings, so each labeling is applied to the
same genes; genome-wide statistics (fraction significant, fraction with k > 1, median log2 k)
are then compared with the real labeling on the same HOGs. RELAX foreground is the
non-orb-weavers, i.e. the tips with value 0 in a permulation row.

1. `python make_permuted_relax_manifest.py --outdir DIR` -- samples HOGs (default: runtime tier A,
   60 HOGs) and labelings (default: 30) and writes `DIR/manifest.tsv` and `DIR/hogs.txt`
2. `sbatch --array=1-<n_tasks> run_relax_permuted.sh DIR/manifest.tsv WD OUT_ROOT HYPHY_ANALYSES` --
   one RELAX run per (HOG, labeling), reusing the real run's alignment and IQ-TREE tree
3. `python summarize_permuted_relax.py --perm-root OUT_ROOT --hogs DIR/hogs.txt`

Before submitting the full array, check that one permuted tree under `OUT_ROOT/<perm_id>/` has the
expected `{Foreground}` labels.

`summarize_permuted_relax.py` needs every sampled HOG in every labeling. An alternative is
`compare_permuted_relax.py`: `extract --perm-root OUT_ROOT --out perm_results.tsv` (standard library only, can
run on the cluster) writes one row per finished run, and `compare --results perm_results.tsv --hogs hogs.txt`
reports the same statistics plus the fractions that are significant and intensified or relaxed.

Result (60 HOGs x 30 labelings, 4 Nov 2026; `results/relax_permuted/`): the real labeling lies inside the permuted null on
every statistic (fraction p <= 0.05: real 0.417, permuted mean 0.369, p = 0.32; significant and intensified 0.317 vs 0.212,
p = 0.16). The original sample included `N5.HOG0058857` (a 14.6 h outlier); it was replaced by `N5.HOG0055575` using the
same 30 labelings (`results/relax_permuted/hogs.txt` is the final list). With 60 HOGs and 30 labelings the test has
low power.

## Dependencies

The module requires:
- pandas
- numpy  
- Access to orb-selection/data/ files
- The orb-selection src/id_converter.py module

All standard dependencies should already be available in the orb-selection environment.
