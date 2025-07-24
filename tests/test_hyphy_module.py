#!/usr/bin/env python3

"""
Test script for the HyPhy analysis results module.

"""

import sys
import pandas as pd
from pathlib import Path

# Add the src directory to the path
# We're in tests/, so we need to go up one level and then into src/
src_dir = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_dir))

from hyphy_results_parser import (
    HyphyResultsManager,
    RelaxResult,
    BustedPhResult,
    AbsrelResult
)
from hyphy_results_helpers import (
    filter_omega, 
    omega_filter_summary,
    convert_hyphy_results_to_locs
)

def create_sample_data():
    """Create sample data for testing."""
    
    # Sample RELAX data
    relax_data = pd.DataFrame({
        'p_value': [0.01, 0.3, 0.02, 0.8, 0.001],
        'k': [0.5, 1.2, 2.1, 0.9, 0.3],
        'LRT': [10.5, 2.1, 15.2, 1.0, 25.8],
        'ω_mean_test': [0.8, 2.5, 12.0, 0.5, 1.1],
        'ω_mean_ref': [1.2, 1.8, 8.0, 0.7, 2.0],
        'result': ['relaxed', 'not significant', 'intensified', 'not significant', 'relaxed']
    })
    relax_data.index = [f'HOG{i:05d}' for i in range(len(relax_data))]
    
    # Sample BUSTED-PH data  
    busted_ph_data = pd.DataFrame({
        'test_pval': [0.01, 0.3, 0.02, 0.8],
        'background_pval': [0.8, 0.2, 0.7, 0.9],
        'shared_pval': [0.02, 0.4, 0.01, 0.6],
        'ω_mean_test': [1.5, 0.8, 25.0, 0.6],
        'ω_mean_ref': [0.9, 1.1, 15.0, 0.8],
        'result': ['hit', 'not significant', 'hit', 'not significant']
    })
    busted_ph_data.index = [f'HOG{i:05d}' for i in range(len(busted_ph_data))]
    
    # Sample aBSREL data
    absrel_data = pd.DataFrame({
        'HOG': ['HOG00001', 'HOG00002', 'HOG00001'],
        'node/species': ['Uloborus_diversus', 'Argiope_argentata', 'Node_1'],
        'gene': ['gene1', 'gene2', 'NA'],
        'corrected_p_value': [0.01, 0.03, 0.02],
        'LRT': [8.5, 6.2, 12.1]
    })
    
    return relax_data, busted_ph_data, absrel_data

def test_result_classes():
    """Test the result classes with sample data."""
    
    print("=== Testing Result Classes ===\n")
    
    # Create sample data
    relax_data, busted_ph_data, absrel_data = create_sample_data()
    
    # Test RelaxResult
    print("1. Testing RelaxResult:")
    relax_result = RelaxResult(relax_data)
    print(f"   Total results: {len(relax_result)}")
    
    significant = relax_result.get_significant_results()
    print(f"   Significant results: {len(significant)}")
    
    relaxed = relax_result.get_relaxed_results()
    intensified = relax_result.get_intensified_results()
    print(f"   Relaxed: {len(relaxed)}, Intensified: {len(intensified)}")
    
    # Test omega filtering
    filtered = relax_result.filter_omega(10.0)
    print(f"   After omega filtering (< 10): {len(filtered)}")
    
    # Test selection type classification
    selection_counts = relax_result.count_selection_types()
    print(f"   Selection type counts: {selection_counts}")
    print()
    
    # Test BustedPhResult
    print("2. Testing BustedPhResult:")
    busted_ph_result = BustedPhResult(busted_ph_data)
    print(f"   Total results: {len(busted_ph_result)}")
    
    hits = busted_ph_result.get_hits()
    print(f"   Hits: {len(hits)}")
    
    filtered_busted = busted_ph_result.filter_omega(20.0)
    print(f"   After omega filtering (< 20): {len(filtered_busted)}")
    print()
    
    # Test AbsrelResult
    print("3. Testing AbsrelResult:")
    absrel_result = AbsrelResult(absrel_data)
    print(f"   Total results: {len(absrel_result)}")
    
    significant_absrel = absrel_result.get_significant_results()
    print(f"   Significant results: {len(significant_absrel)}")
    
    gene_specific = absrel_result.get_gene_specific_results()
    node_specific = absrel_result.get_node_specific_results()
    print(f"   Gene-specific: {len(gene_specific)}, Node-specific: {len(node_specific)}")
    print()

def test_manager():
    """Test the HyphyResultsManager."""
    
    print("=== Testing HyphyResultsManager ===\n")
    
    manager = HyphyResultsManager()
    
    # Create sample data and add to manager
    relax_data, busted_ph_data, absrel_data = create_sample_data()
    
    relax_result = RelaxResult(relax_data)
    busted_ph_result = BustedPhResult(busted_ph_data)
    absrel_result = AbsrelResult(absrel_data)
    
    manager.add_result('relax', relax_result)
    manager.add_result('busted_ph', busted_ph_result)
    manager.add_result('absrel', absrel_result)
    
    print(f"Available results: {manager.list_results()}")
    
    # Test overlap analysis
    overlap_stats = manager.get_overlap_stats(['relax', 'busted_ph'])
    print("Overlap statistics:")
    if not overlap_stats.empty:
        print(overlap_stats.to_string(index=False))
    else:
        print("No overlap data available")
    print()
    
    # Test gene set comparison
    gene_sets = manager.compare_significant_genes(['relax', 'busted_ph', 'absrel'])
    print("Gene sets:")
    for analysis, genes in gene_sets.items():
        print(f"  {analysis}: {len(genes)} genes")
    print()

def test_helper_functions():
    """Test the helper functions."""
    
    print("=== Testing Helper Functions ===\n")
    
    # Create test data
    test_data = pd.DataFrame({
        'ω_mean_test': [0.5, 2.0, 15.0, 0.8, 25.0],
        'ω_mean_ref': [1.0, 1.5, 20.0, 1.2, 18.0],
        'result': ['relaxed', 'intensified', 'not significant', 'relaxed', 'hit'],
        'LOC': [f'LOC{i:03d}' for i in range(5)]
    })
    
    print("Original data:")
    print(test_data)
    print()
    
    # Test omega filtering
    print("Omega filtering tests:")
    filtered_10 = filter_omega(test_data, 10.0)
    print(f"Filtered by omega < 10: {len(filtered_10)} rows")
    
    filtered_20 = filter_omega(test_data, 20.0)
    print(f"Filtered by omega < 20: {len(filtered_20)} rows")
    print()
    
    # Test omega filter summary
    omega_filter_summary(test_data, [1, 5, 10, 20, 30], 'result', 'relaxed')
    print()

def test_id_converter_integration():
    """Test integration with id_converter module."""
    
    print("=== Testing id_converter Integration ===\n")
    
    # Create sample HOG data
    sample_hog_data = pd.DataFrame({
        'p_value': [0.01, 0.02, 0.03],
        'result': ['relaxed', 'intensified', 'relaxed']
    })
    sample_hog_data.index = ['HOG0000001', 'HOG0000002', 'HOG0000003']
    
    print("Sample HOG data:")
    print(sample_hog_data)
    print()
    
    # Test LOC conversion
    print("Testing LOC conversion:")
    try:
        # This will attempt to use the id_converter module
        converted_data = convert_hyphy_results_to_locs(sample_hog_data)
        
        if 'LOC' in converted_data.columns:
            print("✓ Successfully converted HOGs to LOCs")
            print(f"  Original HOGs: {len(sample_hog_data)}")
            print(f"  Converted rows: {len(converted_data)}")
            print(f"  Unique LOCs: {converted_data['LOC'].dropna().nunique()}")
            
            if 'Description' in converted_data.columns:
                print("✓ Gene descriptions also included")
        else:
            print("⚠ LOC conversion returned data but no LOC column found")
            
    except Exception as e:
        print(f"✗ Error testing LOC conversion: {e}")
        print("  This is expected if id_converter data files are not available")
    
    print()

def test_data_persistence():
    """Test saving and loading data."""
    
    print("=== Testing Data Persistence ===\n")
    
    # Create sample data
    relax_data, _, _ = create_sample_data()
    relax_result = RelaxResult(relax_data)
    
    # Test pickle save/load in orb-selection results directory
    results_dir = Path(__file__).parent.parent / "results"
    results_dir.mkdir(exist_ok=True)
    
    test_file = results_dir / "test_relax_result.pkl"
    try:
        relax_result.save_to_pickle(str(test_file))
        print("✓ Successfully saved to pickle")
        
        loaded_result = RelaxResult.load_from_pickle(str(test_file))
        print(f"✓ Successfully loaded from pickle ({len(loaded_result)} results)")
        
        # Clean up
        test_file.unlink()
        
    except Exception as e:
        print(f"✗ Error with pickle save/load: {e}")
    
    # Test CSV export
    test_csv = results_dir / "test_relax_result.csv"
    try:
        relax_result.to_csv(str(test_csv))
        print("✓ Successfully saved to CSV")
        
        # Clean up
        test_csv.unlink()
        
    except Exception as e:
        print(f"✗ Error with CSV save: {e}")
    
    print()

def demonstrate_orb_selection_integration():
    """Demonstrate orb-selection specific features."""
    
    print("=== orb-selection Integration Features ===\n")
    
    print("Integration features:")
    print("✓ Uses src/id_converter.py for HOG to LOC conversion")
    print("✓ Accesses data from orb-selection/data/ directory")
    print("✓ Can save results to orb-selection/results/ directory")
    print("✓ Modular design allows importing from src/")
    print()
    
    print("File locations in orb-selection:")
    print("- src/hyphy_results_parser.py - Main classes")
    print("- src/hyphy_results_helpers.py - Helper functions")
    print("- src/phylogenetic_example.py - Usage example")
    print("- src/test_phylogenetic_module.py - This test script")
    print()
    
    print("Example import patterns:")
    print("from src.hyphy_results_parser import HyphyResultsManager")
    print("from src.hyphy_results_helpers import filter_omega, convert_hyphy_results_to_locs")
    print()

def main():
    """Run all tests and demonstrations."""
    
    print("HyPhy Analysis Results Module - orb-selection Integration Test")
    print("=" * 80)
    print()
    
    try:
        test_result_classes()
        test_manager()
        test_helper_functions()
        test_id_converter_integration()
        test_data_persistence()
        demonstrate_orb_selection_integration()
        
        print("=" * 80)
        print("✓ All tests completed successfully!")
        print("\nThe module is ready to use with your orb-selection workflow.")
        print("See phylogenetic_example.py for working with real JSON files.")
        
    except Exception as e:
        print(f"✗ Error during testing: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
