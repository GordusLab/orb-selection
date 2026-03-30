#!/usr/bin/env python3

"""Parse and summarize HyPhy JSON outputs for downstream analyses.

This script was used to:
1. Load HyPhy JSON result directories into parser classes
2. Filter and summarize result sets
3. Convert HOG hits to LOCs for enrichment workflows
4. Cache parsed outputs as pickle files for reproducibility
5. Compare hit sets across analyses
"""

import os
import sys
import importlib.util
from pathlib import Path

# Add repo src plus this stage directory to the path.
repo_root = Path(__file__).parent.parent.parent
src_dir = repo_root / "src"
stage_dir = Path(__file__).parent
for path in (src_dir, stage_dir):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


def _load_module(module_name: str, module_path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(module_path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module {module_name} from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


hyphy_results_parser = _load_module(
    "hyphy_results_parser", stage_dir / "hyphy_results_parser.py"
)
hyphy_results_helpers = _load_module(
    "hyphy_results_helpers", stage_dir / "hyphy_results_helpers.py"
)

HyphyResultsManager = hyphy_results_parser.HyphyResultsManager
omega_filter_summary = hyphy_results_helpers.omega_filter_summary
convert_hyphy_results_to_locs = hyphy_results_helpers.convert_hyphy_results_to_locs

def main():
    """Main parser workflow used in this project."""
    
    print("=== HyPhy Analysis Results Parsing (orb-selection) ===\n")
    
    # Initialize the results manager
    manager = HyphyResultsManager()
    
    # Data paths in the orb-selection repository
    # Point to the actual data directory in orb-selection
    data_dir = repo_root / "data"
    relax_path = str(data_dir / "relax")
    busted_ph_path = str(data_dir / "busted_ph") 
    busted_ph_rev_path = str(data_dir / "busted_ph_rev")
    absrel_path = str(data_dir / "absrel_rerun")
    
    # Check if directories exist
    paths_exist = {
        'relax': os.path.exists(relax_path),
        'busted_ph': os.path.exists(busted_ph_path),
        'busted_ph_rev': os.path.exists(busted_ph_rev_path),
        'absrel': os.path.exists(absrel_path)
    }
    
    print("Available data paths:")
    for analysis, exists in paths_exist.items():
        status = "✓" if exists else "✗"
        print(f"  {status} {analysis}: {locals()[f'{analysis}_path']}")
    print()
    
    # Load results from available paths
    results_loaded = []
    
    if paths_exist['relax']:
        try:
            print("Loading RELAX results...")
            relax_result = manager.load_relax_from_json(relax_path)
            results_loaded.append('relax')
            print(f"  ✓ Loaded {len(relax_result)} RELAX results")
            
            # Demonstrate RELAX-specific functionality
            print("  RELAX Analysis:")
            significant = relax_result.get_significant_results()
            print(f"    - Significant results: {len(significant)}")
            
            relaxed = relax_result.get_relaxed_results()
            intensified = relax_result.get_intensified_results()
            print(f"    - Relaxed: {len(relaxed)}, Intensified: {len(intensified)}")
            
            # Convert to LOCs using orb-selection id_converter
            print("    - Converting HOGs to LOCs using orb-selection id_converter...")
            try:
                relax_with_locs = convert_hyphy_results_to_locs(significant)
                if 'LOC' in relax_with_locs.columns:
                    unique_locs = relax_with_locs['LOC'].dropna().nunique()
                    print(f"      ✓ Converted to {unique_locs} unique LOCs")
                else:
                    print("      ✗ LOC conversion failed")
            except Exception as e:
                print(f"      ✗ Error converting to LOCs: {e}")
            
            # Filter by omega values
            filtered = relax_result.filter_omega(10000)
            print(f"    - After omega filtering (< 10000): {len(filtered)}")
            
            # Show omega filtering summary
            print("    - Omega filtering summary:")
            omega_filter_summary(relax_result.results_df, [10, 100, 1000, 10000], 'result', 'relaxed')
            
            # Selection type classification
            selection_counts = relax_result.count_selection_types()
            print("    - Selection type counts:")
            for sel_type, count in selection_counts.items():
                print(f"      {sel_type}: {count}")
            
        except Exception as e:
            print(f"  ✗ Error loading RELAX results: {e}")
        print()
    
    if paths_exist['busted_ph']:
        try:
            print("Loading BUSTED-PH results...")
            busted_ph_result = manager.load_busted_ph_from_json(busted_ph_path, name='busted_ph')
            results_loaded.append('busted_ph')
            print(f"  ✓ Loaded {len(busted_ph_result)} BUSTED-PH results")
            
            # Demonstrate BUSTED-PH-specific functionality
            hits = busted_ph_result.get_hits()
            non_sig = busted_ph_result.get_non_significant()
            print(f"    - Hits: {len(hits)}, Non-significant: {len(non_sig)}")
            
            # Convert hits to LOCs
            if len(hits) > 0:
                print("    - Converting BUSTED-PH hits to LOCs...")
                try:
                    hits_with_locs = convert_hyphy_results_to_locs(hits)
                    if 'LOC' in hits_with_locs.columns:
                        unique_locs = hits_with_locs['LOC'].dropna().nunique()
                        print(f"      ✓ Converted {len(hits)} hits to {unique_locs} unique LOCs")
                    else:
                        print("      ✗ LOC conversion failed")
                except Exception as e:
                    print(f"      ✗ Error converting hits to LOCs: {e}")
            
            # Filter by omega values
            filtered = busted_ph_result.filter_omega(10000)
            print(f"    - After omega filtering (< 10000): {len(filtered)}")
            
        except Exception as e:
            print(f"  ✗ Error loading BUSTED-PH results: {e}")
        print()
    
    if paths_exist['busted_ph_rev']:
        try:
            print("Loading BUSTED-PH-REV results...")
            busted_ph_rev_result = manager.load_busted_ph_from_json(busted_ph_rev_path, name='busted_ph_rev')
            results_loaded.append('busted_ph_rev')
            print(f"  ✓ Loaded {len(busted_ph_rev_result)} BUSTED-PH-REV results")
            
            # Demonstrate BUSTED-PH-REV-specific functionality
            hits_rev = busted_ph_rev_result.get_hits()
            non_sig_rev = busted_ph_rev_result.get_non_significant()
            print(f"    - Hits: {len(hits_rev)}, Non-significant: {len(non_sig_rev)}")
            
            # Convert hits to LOCs
            if len(hits_rev) > 0:
                print("    - Converting BUSTED-PH-REV hits to LOCs...")
                try:
                    hits_rev_with_locs = convert_hyphy_results_to_locs(hits_rev)
                    if 'LOC' in hits_rev_with_locs.columns:
                        unique_locs_rev = hits_rev_with_locs['LOC'].dropna().nunique()
                        print(f"      ✓ Converted {len(hits_rev)} hits to {unique_locs_rev} unique LOCs")
                    else:
                        print("      ✗ LOC conversion failed")
                except Exception as e:
                    print(f"      ✗ Error converting hits to LOCs: {e}")
            
            # Filter by omega values
            filtered_rev = busted_ph_rev_result.filter_omega(10000)
            print(f"    - After omega filtering (< 10000): {len(filtered_rev)}")
            
        except Exception as e:
            print(f"  ✗ Error loading BUSTED-PH-REV results: {e}")
        print()
    
    if paths_exist['absrel']:
        try:
            print("Loading aBSREL results...")
            absrel_result = manager.load_absrel_from_json(absrel_path)
            results_loaded.append('absrel')
            print(f"  ✓ Loaded {len(absrel_result)} aBSREL results")
            
            # Demonstrate aBSREL-specific functionality
            significant = absrel_result.get_significant_results()
            print(f"    - Significant results: {len(significant)}")
            
            gene_specific = absrel_result.get_gene_specific_results()
            node_specific = absrel_result.get_node_specific_results()
            print(f"    - Gene-specific: {len(gene_specific)}, Node-specific: {len(node_specific)}")
            
            # Convert significant results to LOCs if any exist
            if len(significant) > 0:
                print("    - Converting significant aBSREL results to LOCs...")
                try:
                    # For aBSREL, we need to create a DataFrame with HOG index from the HOG column
                    hog_df = significant.set_index('HOG') if 'HOG' in significant.columns else significant
                    absrel_with_locs = convert_hyphy_results_to_locs(hog_df)
                    if 'LOC' in absrel_with_locs.columns:
                        unique_locs = absrel_with_locs['LOC'].dropna().nunique()
                        print(f"      ✓ Converted to {unique_locs} unique LOCs")
                except Exception as e:
                    print(f"      ✗ Error converting to LOCs: {e}")
            
        except Exception as e:
            print(f"  ✗ Error loading aBSREL results: {e}")
        print()
    
    # Demonstrate cross-analysis comparisons if multiple results loaded
    if len(results_loaded) > 1:
        print("=== Cross-Analysis Comparisons ===")
        try:
            overlap_stats = manager.get_overlap_stats(results_loaded)
            print("Overlap statistics between analyses:")
            if not overlap_stats.empty:
                print(overlap_stats.to_string(index=False))
            else:
                print("No overlap data available")
            print()
            
            gene_sets = manager.compare_significant_genes(results_loaded)
            print("Gene set sizes:")
            for analysis, genes in gene_sets.items():
                print(f"  {analysis}: {len(genes)} genes")
            
        except Exception as e:
            print(f"Error in cross-analysis comparison: {e}")
        print()
    
    # Demonstrate saving/loading functionality
    if results_loaded:
        # Save in the orb-selection results directory
        save_dir = str(Path(__file__).parent.parent / "results" / "hyphy_results_cache")
        print(f"=== Saving Results to {save_dir} ===")
        try:
            manager.save_all_results(save_dir)
            print("✓ Results saved successfully")
            
            # Test loading
            new_manager = HyphyResultsManager()
            new_manager.load_all_results_from_directory(save_dir)
            loaded_results = new_manager.list_results()
            print(f"✓ Loaded results: {loaded_results}")
            
        except Exception as e:
            print(f"✗ Error saving/loading results: {e}")
        print()
    
    # Show example of accessing data
    print("=== Data Access ===")
    for result_name in results_loaded:
        result = manager.get_result(result_name)
        print(f"{result_name.upper()} Result:")
        print(f"  Type: {type(result).__name__}")
        print(f"  Shape: {result.results_df.shape}")
        print(f"  Columns: {list(result.results_df.columns)[:5]}{'...' if len(result.results_df.columns) > 5 else ''}")
        
        # Show summary stats for numeric columns
        if len(result.results_df.select_dtypes(include=['number']).columns) > 0:
            print("  Numeric column summary:")
            summary = result.get_summary_stats()
            print(f"    Mean values: {summary.loc['mean'].head(3).to_dict()}")
        print()
    
    
    print("=== Parsing Complete ===")
    print("\nReusable parser workflow summary:")
    print("1. Import classes: from hyphy_results_parser import HyphyResultsManager")
    print("2. Create a manager: manager = HyphyResultsManager()")
    print("3. Load your data: result = manager.load_relax_from_json('/path/to/json/files')")
    print("4. Convert to LOCs: from hyphy_results_helpers import convert_hyphy_results_to_locs")
    print("5. Analyze: significant = result.get_significant_results()")
    print("6. Save for reuse: result.save_to_pickle('my_results.pkl')")

if __name__ == "__main__":
    main()
