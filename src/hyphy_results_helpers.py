#!/usr/bin/env python3

"""
Helper functions for HyPhy analysis results, particularly
for filtering and extracting LOC IDs.
"""

import os
import sys
from typing import Tuple, List, Union
import pandas as pd
import numpy as np

# Add the src directory to the path for imports
src_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, src_dir)

try:
    from id_converter import main as convert_hogs_to_locs
except ImportError:
    print("Warning: id_converter module not found. LOC conversion functions will not work.")
    convert_hogs_to_locs = None


def filter_omega(df: pd.DataFrame, omega_threshold: float) -> pd.DataFrame:
    """
    Filter dataframe by omega values for both test and reference branches.
    
    Args:
        df: DataFrame with ω_mean_test and ω_mean_ref columns
        omega_threshold: Maximum omega value threshold
        
    Returns:
        Filtered DataFrame
    """
    return df[(df['ω_mean_test'] < omega_threshold) & (df['ω_mean_ref'] < omega_threshold)]


def filter_omega_rev(df: pd.DataFrame, omega_threshold: float) -> pd.DataFrame:
    """
    Reverse omega filter - keep rows with omega values >= threshold.
    
    Args:
        df: DataFrame with ω_mean_test and ω_mean_ref columns
        omega_threshold: Minimum omega value threshold
        
    Returns:
        Filtered DataFrame
    """
    return df[(df['ω_mean_test'] >= omega_threshold) | (df['ω_mean_ref'] >= omega_threshold)]


def get_fltrd_LOCs(LOCs_df: pd.DataFrame, omega: float, test: str) -> Union[Tuple, None]:
    """
    Get filtered LOC lists based on omega filtering and test type.
    
    Args:
        LOCs_df: DataFrame with LOC column and results
        omega: Omega threshold for filtering
        test: Type of test ('relax' or 'busted-ph')
        
    Returns:
        Tuple of LOC arrays based on test type
    """
    all_LOCs = LOCs_df['LOC'].dropna().unique()
    df_fltrd = filter_omega(LOCs_df, omega)
    all_LOCs_fltrd = df_fltrd['LOC'].dropna().unique()
    df_ns = df_fltrd[df_fltrd['result'] == 'not significant']
    ns_LOCs = df_ns['LOC'].dropna().unique()

    if test == 'relax':
        df_rel = df_fltrd[df_fltrd['result'] == 'relaxed']
        df_int = df_fltrd[df_fltrd['result'] == 'intensified']
        df_hits = df_fltrd[(df_fltrd['result'] == 'intensified') | (df_fltrd['result'] == 'relaxed')]
        rel_LOCs = df_rel['LOC'].dropna().unique()
        int_LOCs = df_int['LOC'].dropna().unique()
        hit_LOCs = df_hits['LOC'].dropna().unique()

        return all_LOCs, all_LOCs_fltrd, hit_LOCs, ns_LOCs, rel_LOCs, int_LOCs
    
    elif test == 'busted-ph':
        df_hits = df_fltrd[df_fltrd['result'] == 'hit']
        hit_LOCs = df_hits['LOC'].dropna().unique()

        return all_LOCs, all_LOCs_fltrd, hit_LOCs, ns_LOCs
    
    else:
        raise ValueError(f"Test type '{test}' not supported. Use 'relax' or 'busted-ph'.")


def count_rel_int(df: pd.DataFrame) -> pd.DataFrame:
    """
    Count and classify relaxation/intensification types and add selection type column.
    
    Args:
        df: DataFrame with RELAX results
        
    Returns:
        DataFrame with added 'type of selection' column
    """
    df = df.copy()
    
    # Count different types
    df_relax_of_pos = df[(df['result'] == 'relaxed') & (df['ω_mean_test'] < df['ω_mean_ref'])]
    df_relax_of_neg = df[(df['result'] == 'relaxed') & (df['ω_mean_test'] > df['ω_mean_ref'])]
    i = len(df_relax_of_pos)
    j = len(df_relax_of_neg)
    print(f'{i} relaxation of positive selection, {j} relaxation of negative selection')

    df_int_of_pos = df[(df['result'] == 'intensified') & (df['ω_mean_test'] > df['ω_mean_ref'])]
    df_int_of_neg = df[(df['result'] == 'intensified') & (df['ω_mean_test'] < df['ω_mean_ref'])]
    i = len(df_int_of_pos)
    j = len(df_int_of_neg)
    print(f'{i} intensification of positive selection, {j} intensification of negative selection')

    # Add classification column
    conditions = [
        df['result'].eq('relaxed') & (df['ω_mean_test'] < df['ω_mean_ref']),
        df['result'].eq('relaxed') & (df['ω_mean_test'] > df['ω_mean_ref']),
        df['result'].eq('intensified') & (df['ω_mean_test'] < df['ω_mean_ref']),
        df['result'].eq('intensified') & (df['ω_mean_test'] > df['ω_mean_ref'])
    ]

    choices = ['positive', 'negative', 'negative', 'positive']
    df['type of selection'] = np.select(conditions, choices, default='')

    return df


def get_universe_LOCs(NX_genecount_tsv: str, thresholds: List[int], output_path: str) -> None:
    """
    Generate universe LOC files for topGO analysis at different occupancy thresholds.
    Uses the id_converter module from orb-selection.
    
    Args:
        NX_genecount_tsv: Path to gene count TSV file
        thresholds: List of occupancy thresholds
        output_path: Output directory path
    """
    if convert_hogs_to_locs is None:
        print("Error: id_converter module not available. Cannot generate universe LOCs.")
        return
    
    try:
        genecount_df = pd.read_csv(NX_genecount_tsv, sep='\t', index_col='HOG')
        
        # Use the id_converter to get LOCs - assumes N5.tsv is in ../data/
        data_dir = os.path.join(src_dir, '..', 'data')
        n5_tsv_path = os.path.join(data_dir, 'N5.tsv')
        
        LOCs_df = convert_hogs_to_locs(genecount_df, n5_tsv_path)
        LOCs_df = LOCs_df.set_index('HOG')
        LOCs_df['occupancy'] = LOCs_df.select_dtypes(include='number').astype('bool').sum(axis=1)

        for threshold in thresholds:
            LOCs_fltrd = LOCs_df[LOCs_df['occupancy'] >= threshold]
            LOC_list = LOCs_fltrd['LOC'].dropna().unique()
            
            os.makedirs(output_path, exist_ok=True)
            with open(f'{output_path}/occ_{threshold}_universe.txt', 'w+') as f:
                for line in LOC_list:
                    f.write(f'{line}\n')
            
            print(f"Generated universe LOC file for occupancy >= {threshold}: {len(LOC_list)} LOCs")
            
    except Exception as e:
        print(f"Error generating universe LOCs: {e}")


def calculate_mean_omega(df: pd.DataFrame, branch_type: str = 'both') -> pd.DataFrame:
    """
    Calculate mean omega values from rate distributions.
    
    Args:
        df: DataFrame with omega and proportion columns
        branch_type: 'test', 'ref', or 'both' branches to calculate
        
    Returns:
        DataFrame with added mean omega columns
    """
    df = df.copy()
    
    if branch_type in ['test', 'both']:
        omega_cols_test = ['ω1_test', 'ω2_test', 'ω3_test']
        prop_cols_test = ['ω1_test_P', 'ω2_test_P', 'ω3_test_P']
        
        if all(col in df.columns for col in omega_cols_test + prop_cols_test):
            df['ω_mean_test'] = (
                df['ω1_test'] * df['ω1_test_P'] +
                df['ω2_test'] * df['ω2_test_P'] +
                df['ω3_test'] * df['ω3_test_P']
            )
    
    if branch_type in ['ref', 'both']:
        omega_cols_ref = ['ω1_ref', 'ω2_ref', 'ω3_ref']
        prop_cols_ref = ['ω1_ref_P', 'ω2_ref_P', 'ω3_ref_P']
        
        if all(col in df.columns for col in omega_cols_ref + prop_cols_ref):
            df['ω_mean_ref'] = (
                df['ω1_ref'] * df['ω1_ref_P'] +
                df['ω2_ref'] * df['ω2_ref_P'] +
                df['ω3_ref'] * df['ω3_ref_P']
            )
    
    return df


def omega_filter_summary(df: pd.DataFrame, omega_thresholds: List[float], 
                        result_column: str = 'result', hit_value: str = 'hit') -> None:
    """
    Print summary of how many results remain after different omega filtering thresholds.
    
    Args:
        df: DataFrame to filter
        omega_thresholds: List of omega thresholds to test
        result_column: Column name containing results
        hit_value: Value in result_column considered a "hit"
    """
    print("Omega filtering summary:")
    print("=" * 40)
    
    for omega in omega_thresholds:
        filtered_df = filter_omega(df, omega)
        total_remaining = len(filtered_df)
        
        if hit_value in filtered_df[result_column].values:
            hits_remaining = len(filtered_df[filtered_df[result_column] == hit_value])
            print(f"Max omega = {omega:>8}: {total_remaining:>5} total, {hits_remaining:>5} hits")
        else:
            print(f"Max omega = {omega:>8}: {total_remaining:>5} total")
    
    print("=" * 40)


def convert_hyphy_results_to_locs(results_df: pd.DataFrame, 
                                hog_node_genes_tsv: str = None) -> pd.DataFrame:
    """
    Convert HyPhy analysis results to include LOC information using id_converter.
    
    Args:
        results_df: DataFrame with HOG indices from HyPhy analysis
        hog_node_genes_tsv: Path to hierarchical orthogroup file (optional, defaults to N5.tsv)
        
    Returns:
        DataFrame with added LOC and gene description columns
    """
    if convert_hogs_to_locs is None:
        print("Error: id_converter module not available. Cannot convert to LOCs.")
        return results_df
    
    try:
        # Use default N5.tsv path if not provided
        if hog_node_genes_tsv is None:
            data_dir = os.path.join(src_dir, '..', 'data')
            hog_node_genes_tsv = os.path.join(data_dir, 'N5.tsv')
        
        # Check if the required files exist
        if not os.path.exists(hog_node_genes_tsv):
            print(f"Warning: HOG node genes file not found: {hog_node_genes_tsv}")
            return results_df
        
        # Convert using the id_converter module
        results_with_locs = convert_hogs_to_locs(results_df, hog_node_genes_tsv)
        
        print(f"Successfully converted {len(results_df)} HOGs to LOCs using id_converter")
        return results_with_locs
        
    except Exception as e:
        print(f"Error converting to LOCs: {e}")
        print("This may be due to missing data files or HOGs not found in the reference.")
        return results_df


# Example usage
if __name__ == "__main__":
    # Example of using helper functions
    print("Helper functions for HyPhy analysis results")
    
    # Create example data
    example_data = pd.DataFrame({
        'HOG': ['HOG001', 'HOG002', 'HOG003'],
        'ω_mean_test': [0.5, 2.0, 15.0],
        'ω_mean_ref': [1.0, 1.5, 20.0],
        'result': ['relaxed', 'intensified', 'not significant'],
        'LOC': ['LOC001', 'LOC002', 'LOC003']
    })
    
    print("\nExample DataFrame:")
    print(example_data)
    
    print("\nFiltered by omega < 10:")
    filtered = filter_omega(example_data, 10)
    print(filtered)
    
    print("\nOmega filtering summary:")
    omega_filter_summary(example_data, [1, 5, 10, 100], 'result', 'relaxed')
