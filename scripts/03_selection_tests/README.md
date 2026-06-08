# HyPhy Results Analysis Module

This module provides classes to store and manage results from HyPhy selection analyses (RELAX, aBSREL, BUSTED-PH).

## Files

- `scripts/03_selection_tests/hyphy_results_parser.py` - Main module with result classes
- `scripts/03_selection_tests/hyphy_results_helpers.py` - Helper functions
- `scripts/03_selection_tests/parse_hyphy_results.py` - Project parser workflow script
- `scripts/03_selection_tests/README.md` - This documentation

## Quick Start

```python
# Add repo src + stage directory to sys.path, then import modules
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
manager.save_all_results('results/hyphy_cache/')
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
manager.save_all_results('results/hyphy_cache/')

# Results are saved as:
# results/hyphy_cache/relax_results.pkl
# results/hyphy_cache/busted_ph_results.pkl
# results/hyphy_cache/absrel_results.pkl
```

## Classes

### HyphyResultsManager
Central manager for handling multiple analysis results.

**Key Methods:**
- `load_relax_from_json(directory)` - Load RELAX results
- `load_busted_ph_from_json(directory)` - Load BUSTED-PH results  
- `load_absrel_from_json(directory)` - Load aBSREL results
- `compare_significant_genes(analyses)` - Compare gene sets
- `get_overlap_stats(analyses)` - Get overlap statistics
- `save_all_results(directory)` - Save all results to cache

### RelaxResult, BustedPhResult, AbsrelResult
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

# aBSREL-specific
absrel_result = AbsrelResult(absrel_df)
gene_specific = absrel_result.get_gene_specific_results()
species_results = absrel_result.get_results_by_species('Uloborus_diversus')
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
omega_filter_summary(relax_result.results_df, [10, 100, 1000, 10000])

# Get specific result types
relaxed = relax_result.get_relaxed_results()
intensified = relax_result.get_intensified_results()
```

### Cross-Analysis Comparisons
```python
# Load multiple analyses
manager.load_relax_from_json('/path/to/relax/')
manager.load_busted_ph_from_json('/path/to/busted_ph/')

# Compare significant genes
overlap_stats = manager.get_overlap_stats(['relax', 'busted_ph'])
gene_sets = manager.compare_significant_genes(['relax', 'busted_ph'])

print(f"RELAX hits: {len(gene_sets['relax'])}")
print(f"BUSTED-PH hits: {len(gene_sets['busted_ph'])}")
print(f"Overlap: {len(gene_sets['relax'] & gene_sets['busted_ph'])}")
```

## Example Workflow

```python
#!/usr/bin/env python3

from hyphy_results_parser import HyphyResultsManager
from hyphy_results_helpers import convert_hyphy_results_to_locs, filter_omega

# 1. Initial setup (run once)
manager = HyphyResultsManager()
relax_result = manager.load_relax_from_json('/path/to/relax/jsons/')
busted_ph_result = manager.load_busted_ph_from_json('/path/to/busted_ph/jsons/')

# 2. Cache results for future use
manager.save_all_results('results/hyphy_cache/')

# 3. Analysis sessions (load from cache)
manager = HyphyResultsManager()
manager.load_all_results_from_directory('results/hyphy_cache/')
relax_result = manager.get_result('relax')

# 4. Filter and analyze
significant = relax_result.get_significant_results()
filtered = relax_result.filter_omega(10000)

# 5. Convert to LOCs using orb-selection id_converter
results_with_locs = convert_hyphy_results_to_locs(filtered.results_df)

# 6. Save final results
results_with_locs.to_csv('results/relax_significant_with_locs.csv', index=False)
```

## Dependencies

The module requires:
- pandas
- numpy  
- Access to orb-selection/data/ files
- The orb-selection src/id_converter.py module

All standard dependencies should already be available in the orb-selection environment.
