#!/usr/bin/env python3

"""
Phylogenetic Analysis Results Module

This module provides classes to store and manage results from phylogenetic selection analyses
(RELAX, aBSREL, BUSTED-PH) to avoid re-parsing JSON files repeatedly.

Classes:
    - HyphyResult: Base class for all analysis results
    - RelaxResult: Storage and analysis of RELAX results
    - AbsrelResult: Storage and analysis of aBSREL results  
    - BustedPhResult: Storage and analysis of BUSTED-PH results
    - HyphyResultsManager: Manager class to handle multiple analysis types
"""

import json
import os
from typing import List, Dict, Optional
from abc import ABC, abstractmethod
import pickle
from pathlib import Path
import pandas as pd
import numpy as np

class HyphyResult(ABC):
    """
    Abstract base class for phylogenetic analysis results.
    
    This class provides common functionality for storing and manipulating
    results from different phylogenetic selection analysis tools.
    """
    
    def __init__(self, analysis_type: str, results_df: pd.DataFrame):
        """
        Initialize the analysis result.
        
        Args:
            analysis_type: Type of analysis ('relax', 'absrel', 'busted-ph')
            results_df: DataFrame containing the parsed results
        """
        self.analysis_type = analysis_type
        self.results_df = results_df.copy()
        self._original_df = results_df.copy()  # Keep original for reset functionality
        
    def __len__(self) -> int:
        """Return the number of results."""
        return len(self.results_df)
    
    def __repr__(self) -> str:
        """String representation of the results."""
        return f"{self.__class__.__name__}(n_results={len(self)}, analysis_type='{self.analysis_type}')"
    
    def get_summary_stats(self) -> pd.Series:
        """Get summary statistics for numeric columns."""
        return self.results_df.select_dtypes(include=[np.number]).describe()
    
    def filter_omega(self, max_omega: float) -> 'HyphyResult':
        """
        Filter results by maximum omega values for both test and reference branches.
        
        Args:
            max_omega: Maximum omega value threshold
            
        Returns:
            New instance with filtered results
        """
        if hasattr(self, '_has_omega_columns') and self._has_omega_columns:
            filtered_df = self.results_df[
                (self.results_df['ω_mean_test'] < max_omega) & 
                (self.results_df['ω_mean_ref'] < max_omega)
            ]
            return self.__class__(filtered_df)
        else:
            return self
    
    def filter_omega_reverse(self, min_omega: float) -> 'HyphyResult':
        """
        Filter results by minimum omega values (exclude low omega values).
        
        Args:
            min_omega: Minimum omega value threshold
            
        Returns:
            New instance with filtered results
        """
        if hasattr(self, '_has_omega_columns') and self._has_omega_columns:
            filtered_df = self.results_df[
                (self.results_df['ω_mean_test'] >= min_omega) | 
                (self.results_df['ω_mean_ref'] >= min_omega)
            ]
            return self.__class__(filtered_df)
        else:
            return self
    
    def reset_filters(self) -> None:
        """Reset DataFrame to original unfiltered state."""
        self.results_df = self._original_df.copy()
    
    def save_to_pickle(self, filepath: str) -> None:
        """Save the results object to a pickle file."""
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
    
    @classmethod
    def load_from_pickle(cls, filepath: str) -> 'HyphyResult':
        """Load results object from a pickle file."""
        with open(filepath, 'rb') as f:
            return pickle.load(f)
    
    def to_csv(self, filepath: str, **kwargs) -> None:
        """Save results DataFrame to CSV."""
        self.results_df.to_csv(filepath, **kwargs)
    
    @abstractmethod
    def get_significant_results(self, alpha: float = 0.05) -> pd.DataFrame:
        """Get statistically significant results."""
        pass


class RelaxResult(HyphyResult):
    """
    Storage and analysis class for RELAX results.
    
    RELAX tests for relaxation or intensification of selection pressure
    between test and reference branches.
    """
    
    def __init__(self, results_df: pd.DataFrame):
        """
        Initialize RELAX results.
        
        Args:
            results_df: DataFrame with RELAX results
        """
        super().__init__('relax', results_df)
        self._has_omega_columns = True
        
        # Calculate mean omega values if not present
        if 'ω_mean_test' not in self.results_df.columns:
            self._calculate_mean_omega()
    
    def _calculate_mean_omega(self) -> None:
        """Calculate mean omega values from rate distributions."""
        omega_cols_test = ['ω1_test', 'ω2_test', 'ω3_test']
        prop_cols_test = ['ω1_test_P', 'ω2_test_P', 'ω3_test_P']
        omega_cols_ref = ['ω1_ref', 'ω2_ref', 'ω3_ref']
        prop_cols_ref = ['ω1_ref_P', 'ω2_ref_P', 'ω3_ref_P']
        
        if all(col in self.results_df.columns for col in omega_cols_test + prop_cols_test):
            self.results_df['ω_mean_test'] = (
                self.results_df['ω1_test'] * self.results_df['ω1_test_P'] +
                self.results_df['ω2_test'] * self.results_df['ω2_test_P'] +
                self.results_df['ω3_test'] * self.results_df['ω3_test_P']
            )
            
        if all(col in self.results_df.columns for col in omega_cols_ref + prop_cols_ref):
            self.results_df['ω_mean_ref'] = (
                self.results_df['ω1_ref'] * self.results_df['ω1_ref_P'] +
                self.results_df['ω2_ref'] * self.results_df['ω2_ref_P'] +
                self.results_df['ω3_ref'] * self.results_df['ω3_ref_P']
            )
    
    def get_significant_results(self, alpha: float = 0.05) -> pd.DataFrame:
        """
        Get statistically significant RELAX results.
        
        Args:
            alpha: Significance threshold (default: 0.05)
            
        Returns:
            DataFrame with significant results (relaxed or intensified)
        """
        return self.results_df[
            (self.results_df['p_value'] <= alpha) &
            (self.results_df['result'].isin(['relaxed', 'intensified']))
        ]
    
    def get_relaxed_results(self, alpha: float = 0.05) -> pd.DataFrame:
        """Get genes with relaxed selection."""
        return self.results_df[
            (self.results_df['p_value'] <= alpha) &
            (self.results_df['result'] == 'relaxed')
        ]
    
    def get_intensified_results(self, alpha: float = 0.05) -> pd.DataFrame:
        """Get genes with intensified selection."""
        return self.results_df[
            (self.results_df['p_value'] <= alpha) &
            (self.results_df['result'] == 'intensified')
        ]
    
    def classify_selection_type(self) -> pd.DataFrame:
        """
        Classify whether relaxation/intensification affects positive or negative selection.
        
        Returns:
            DataFrame with additional 'selection_type' column
        """
        df = self.results_df.copy()
        
        conditions = [
            df['result'].eq('relaxed') & (df['ω_mean_test'] < df['ω_mean_ref']),
            df['result'].eq('relaxed') & (df['ω_mean_test'] > df['ω_mean_ref']),
            df['result'].eq('intensified') & (df['ω_mean_test'] < df['ω_mean_ref']),
            df['result'].eq('intensified') & (df['ω_mean_test'] > df['ω_mean_ref'])
        ]
        
        choices = ['positive', 'negative', 'negative', 'positive']
        
        df['selection_type'] = np.select(conditions, choices, default='')
        
        return df
    
    def count_selection_types(self) -> Dict[str, int]:
        """Count different types of selection changes."""
        df = self.classify_selection_type()
        
        counts = {
            'relaxation_of_positive': len(df[
                (df['result'] == 'relaxed') & 
                (df['ω_mean_test'] < df['ω_mean_ref'])
            ]),
            'relaxation_of_negative': len(df[
                (df['result'] == 'relaxed') & 
                (df['ω_mean_test'] > df['ω_mean_ref'])
            ]),
            'intensification_of_positive': len(df[
                (df['result'] == 'intensified') & 
                (df['ω_mean_test'] > df['ω_mean_ref'])
            ]),
            'intensification_of_negative': len(df[
                (df['result'] == 'intensified') & 
                (df['ω_mean_test'] < df['ω_mean_ref'])
            ])
        }
        
        return counts


class AbsrelResult(HyphyResult):
    """
    Storage and analysis class for aBSREL results.
    
    aBSREL detects lineage-specific episodic diversifying selection.
    """
    
    def __init__(self, results_df: pd.DataFrame):
        """
        Initialize aBSREL results.
        
        Args:
            results_df: DataFrame with aBSREL results
        """
        super().__init__('absrel', results_df)
        self._has_omega_columns = False
    
    def get_significant_results(self, alpha: float = 0.05) -> pd.DataFrame:
        """
        Get statistically significant aBSREL results.
        
        Args:
            alpha: Significance threshold (default: 0.05)
            
        Returns:
            DataFrame with significant results
        """
        return self.results_df[self.results_df['corrected_p_value'] <= alpha]
    
    def get_results_by_species(self, species_name: str) -> pd.DataFrame:
        """Get results for a specific species."""
        return self.results_df[
            self.results_df['node/species'].str.contains(species_name, na=False)
        ]
    
    def get_gene_specific_results(self) -> pd.DataFrame:
        """Get results that are gene-specific (not node-specific)."""
        return self.results_df[self.results_df['gene'] != 'NA']
    
    def get_node_specific_results(self) -> pd.DataFrame:
        """Get results that are node-specific."""
        return self.results_df[self.results_df['gene'] == 'NA']


class BustedPhResult(HyphyResult):
    """
    Storage and analysis class for BUSTED-PH results.
    
    BUSTED-PH tests for gene-wide episodic diversifying selection
    associated with a phenotype.
    """
    
    def __init__(self, results_df: pd.DataFrame):
        """
        Initialize BUSTED-PH results.
        
        Args:
            results_df: DataFrame with BUSTED-PH results
        """
        super().__init__('busted-ph', results_df)
        self._has_omega_columns = True
        
        # Calculate mean omega values if not present
        if 'ω_mean_test' not in self.results_df.columns:
            self._calculate_mean_omega()
    
    def _calculate_mean_omega(self) -> None:
        """Calculate mean omega values from rate distributions."""
        omega_cols_test = ['ω1_test', 'ω2_test', 'ω3_test']
        prop_cols_test = ['ω1_test_P', 'ω2_test_P', 'ω3_test_P']
        omega_cols_ref = ['ω1_ref', 'ω2_ref', 'ω3_ref']
        prop_cols_ref = ['ω1_ref_P', 'ω2_ref_P', 'ω3_ref_P']
        
        if all(col in self.results_df.columns for col in omega_cols_test + prop_cols_test):
            self.results_df['ω_mean_test'] = (
                self.results_df['ω1_test'] * self.results_df['ω1_test_P'] +
                self.results_df['ω2_test'] * self.results_df['ω2_test_P'] +
                self.results_df['ω3_test'] * self.results_df['ω3_test_P']
            )
            
        if all(col in self.results_df.columns for col in omega_cols_ref + prop_cols_ref):
            self.results_df['ω_mean_ref'] = (
                self.results_df['ω1_ref'] * self.results_df['ω1_ref_P'] +
                self.results_df['ω2_ref'] * self.results_df['ω2_ref_P'] +
                self.results_df['ω3_ref'] * self.results_df['ω3_ref_P']
            )
    
    def get_significant_results(self, alpha: float = 0.05) -> pd.DataFrame:
        """
        Get statistically significant BUSTED-PH results.
        
        Args:
            alpha: Significance threshold (default: 0.05)
            
        Returns:
            DataFrame with significant hits
        """
        return self.results_df[self.results_df['result'] == 'hit']
    
    def get_hits(self) -> pd.DataFrame:
        """Get BUSTED-PH hits (significant results)."""
        return self.get_significant_results()
    
    def get_non_significant(self) -> pd.DataFrame:
        """Get non-significant results."""
        return self.results_df[self.results_df['result'] == 'not significant']


class HyphyResultsManager:
    """
    Manager class to handle multiple phylogenetic analysis results.
    
    This class provides functionality to load, store, and compare results
    from different phylogenetic analysis tools.
    """
    
    def __init__(self):
        """Initialize the results manager."""
        self.results = {}
    
    def add_result(self, name: str, result: HyphyResult) -> None:
        """
        Add a phylogenetic analysis result.
        
        Args:
            name: Name identifier for the result
            result: HyphyResult instance
        """
        self.results[name] = result
    
    def get_result(self, name: str) -> Optional[HyphyResult]:
        """Get a specific result by name."""
        return self.results.get(name)
    
    def list_results(self) -> List[str]:
        """List all available result names."""
        return list(self.results.keys())
    
    def load_relax_from_json(self, json_directory: str, name: str = 'relax') -> RelaxResult:
        """
        Load RELAX results from JSON files.
        
        Args:
            json_directory: Directory containing RELAX JSON files
            name: Name to store the result under
            
        Returns:
            RelaxResult instance
        """
        json_files = [f for f in os.listdir(json_directory) if f.endswith('.json')]
        relax_results = []
        
        for json_file in json_files:
            with open(os.path.join(json_directory, json_file), 'r') as f:
                data = json.load(f)
                pval = data['test results']['p-value']
                k = data['test results']['relaxation or intensification parameter']
                
                # Determine result classification
                if (float(pval) <= 0.05) & (float(k) < 1):
                    result = 'relaxed'
                elif (float(pval) <= 0.05) & (float(k) > 1):
                    result = 'intensified'
                else:
                    result = 'not significant'
                
                relax_results.append({
                    'HOG': json_file.split('_')[0],
                    'p_value': pval,
                    'k': k,
                    'LRT': data['test results']['LRT'],
                    'MG94xREV_ω_reference': data["fits"]["MG94xREV with separate rates for branch sets"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *Reference*"][0][0],
                    'MG94xREV_ω_test': data["fits"]["MG94xREV with separate rates for branch sets"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *Test*"][0][0],
                    'ω1_test': data["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["0"]["omega"],
                    'ω1_test_P': data["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["0"]["proportion"],
                    'ω2_test': data["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["1"]["omega"],
                    'ω2_test_P': data["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["1"]["proportion"],
                    'ω3_test': data["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["2"]["omega"],
                    'ω3_test_P': data["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["2"]["proportion"],
                    'ω1_ref': data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["0"]["omega"],
                    'ω1_ref_P': data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["0"]["proportion"],
                    'ω2_ref': data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["1"]["omega"],
                    'ω2_ref_P': data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["1"]["proportion"],
                    'ω3_ref': data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["2"]["omega"],
                    'ω3_ref_P': data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["2"]["proportion"],
                    'result': result
                })
        
        relax_df = pd.DataFrame(relax_results).set_index('HOG')
        result_obj = RelaxResult(relax_df)
        self.add_result(name, result_obj)
        return result_obj
    
    def load_busted_ph_from_json(self, json_directory: str, name: str = 'busted_ph') -> BustedPhResult:
        """
        Load BUSTED-PH results from JSON files.
        
        Args:
            json_directory: Directory containing BUSTED-PH JSON files
            name: Name to store the result under
            
        Returns:
            BustedPhResult instance
        """
        json_files = [f for f in os.listdir(json_directory) if f.endswith('.json')]
        busted_ph_results = []
        
        for json_file in json_files:
            with open(os.path.join(json_directory, json_file), 'r') as f:
                data = json.load(f)
                test_pval = data['test results']['p-value']
                background_pval = data['test results background']['p-value']
                shared_pval = data['test results shared distributions']['p-value']
                
                # Determine result classification
                if (float(test_pval) <= 0.05) & (float(background_pval) > 0.05) & (float(shared_pval) <= 0.05):
                    result = 'hit'
                else:
                    result = 'not significant'
                
                busted_ph_results.append({
                    'HOG': json_file.split('_')[0],
                    'test_pval': test_pval,
                    'test_LRT': data['test results']['LRT'],
                    'background_pval': background_pval,
                    'background_LRT': data['test results background']['LRT'],
                    'shared_pval': shared_pval,
                    'shared_LRT': data['test results shared distributions']['LRT'],
                    'MG94xREV_ω_test': data["fits"]["MG94xREV with separate rates for branch sets"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *test*"][0][0],
                    'MG94xREV_ω_ref': data["fits"]["MG94xREV with separate rates for branch sets"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *background*"][0][0],
                    'ω1_test': data["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]["0"]["omega"],
                    'ω1_test_P': data["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]["0"]["proportion"],
                    'ω2_test': data["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]["1"]["omega"],
                    'ω2_test_P': data["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]["1"]["proportion"],
                    'ω3_test': data["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]["2"]["omega"],
                    'ω3_test_P': data["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]["2"]["proportion"],
                    'ω1_ref': data["fits"]["Unconstrained model"]["Rate Distributions"]["Background"]["0"]["omega"],
                    'ω1_ref_P': data["fits"]["Unconstrained model"]["Rate Distributions"]["Background"]["0"]["proportion"],
                    'ω2_ref': data["fits"]["Unconstrained model"]["Rate Distributions"]["Background"]["1"]["omega"],
                    'ω2_ref_P': data["fits"]["Unconstrained model"]["Rate Distributions"]["Background"]["1"]["proportion"],
                    'ω3_ref': data["fits"]["Unconstrained model"]["Rate Distributions"]["Background"]["2"]["omega"],
                    'ω3_ref_P': data["fits"]["Unconstrained model"]["Rate Distributions"]["Background"]["2"]["proportion"],
                    'result': result
                })
        
        busted_ph_df = pd.DataFrame(busted_ph_results).set_index('HOG')
        result_obj = BustedPhResult(busted_ph_df)
        self.add_result(name, result_obj)
        return result_obj
    
    def load_absrel_from_json(self, json_directory: str, name: str = 'absrel') -> AbsrelResult:
        """
        Load aBSREL results from JSON files.
        
        Args:
            json_directory: Directory containing aBSREL JSON files
            name: Name to store the result under
            
        Returns:
            AbsrelResult instance
        """
        json_files = [f for f in os.listdir(json_directory) if f.endswith('.json')]
        absrel_hits = []
        
        for json_file in json_files:
            with open(os.path.join(json_directory, json_file), 'r') as f:
                data = json.load(f)
                
                for gene in data['branch attributes']['0']:
                    if str(data['branch attributes']['0'][gene]['Corrected P-value']) != 'None':
                        pval = data['branch attributes']['0'][gene]['Corrected P-value']
                        rate_classes = int(data['branch attributes']['0'][gene]['Rate classes'])
                        
                        if float(pval) <= 0.05:
                            if "Node" not in gene:
                                absrel_hits.append({
                                    'HOG': json_file.split('_')[0],
                                    'node/species': gene.rsplit('_', 3)[0],
                                    'gene': gene.rsplit('_', 3)[1],
                                    'corrected_p_value': pval,
                                    'non-synonymous_subs/site': data['branch attributes']['0'][gene]['Full adaptive model (non-synonymous subs/site)'],
                                    'synonymous subs/site': data['branch attributes']['0'][gene]['Full adaptive model (synonymous subs/site)'],
                                    'LRT': data['branch attributes']['0'][gene]['LRT'],
                                    'rate_classes': rate_classes,
                                    'rate_distributions': data['branch attributes']['0'][gene]['Rate Distributions']
                                })
                            else:
                                absrel_hits.append({
                                    'HOG': json_file.split('_')[0],
                                    'node/species': gene,
                                    'gene': 'NA',
                                    'corrected_p_value': pval,
                                    'LRT': data['branch attributes']['0'][gene]['LRT'],
                                    'rate_classes': rate_classes,
                                    'rate_distributions': data['branch attributes']['0'][gene]['Rate Distributions']
                                })
        
        absrel_df = pd.DataFrame(absrel_hits)
        result_obj = AbsrelResult(absrel_df)
        self.add_result(name, result_obj)
        return result_obj
    
    def save_all_results(self, directory: str) -> None:
        """Save all results to pickle files in the specified directory."""
        Path(directory).mkdir(parents=True, exist_ok=True)
        
        for name, result in self.results.items():
            filepath = os.path.join(directory, f"{name}_results.pkl")
            result.save_to_pickle(filepath)
    
    def load_all_results_from_directory(self, directory: str) -> None:
        """Load all pickle files from a directory."""
        pickle_files = [f for f in os.listdir(directory) if f.endswith('_results.pkl')]
        
        for pickle_file in pickle_files:
            name = pickle_file.replace('_results.pkl', '')
            filepath = os.path.join(directory, pickle_file)
            result = HyphyResult.load_from_pickle(filepath)
            self.add_result(name, result)
    
    def compare_significant_genes(self, result_names: List[str], alpha: float = 0.05) -> Dict[str, set]:
        """
        Compare significant genes across different analysis results.
        
        Args:
            result_names: Names of results to compare
            alpha: Significance threshold
            
        Returns:
            Dictionary with gene sets for each analysis
        """
        gene_sets = {}
        
        for name in result_names:
            if name in self.results:
                result = self.results[name]
                sig_results = result.get_significant_results(alpha)
                
                if hasattr(sig_results, 'index'):
                    gene_sets[name] = set(sig_results.index)
                elif 'HOG' in sig_results.columns:
                    gene_sets[name] = set(sig_results['HOG'])
                else:
                    gene_sets[name] = set()
        
        return gene_sets
    
    def get_overlap_stats(self, result_names: List[str], alpha: float = 0.05) -> pd.DataFrame:
        """
        Get overlap statistics between different analysis results.
        
        Args:
            result_names: Names of results to compare
            alpha: Significance threshold
            
        Returns:
            DataFrame with overlap statistics
        """
        gene_sets = self.compare_significant_genes(result_names, alpha)
        
        overlap_data = []
        for i, name1 in enumerate(result_names):
            for name2 in result_names[i+1:]:
                if name1 in gene_sets and name2 in gene_sets:
                    set1, set2 = gene_sets[name1], gene_sets[name2]
                    overlap = len(set1 & set2)
                    union = len(set1 | set2)
                    jaccard = overlap / union if union > 0 else 0
                    
                    overlap_data.append({
                        'analysis_1': name1,
                        'analysis_2': name2,
                        'n_genes_1': len(set1),
                        'n_genes_2': len(set2),
                        'overlap': overlap,
                        'union': union,
                        'jaccard_index': jaccard
                    })
        
        return pd.DataFrame(overlap_data)


# Example usage and helper functions
def example_usage():
    """Example of how to use the phylogenetic results classes."""
    
    # Initialize the manager
    # manager = HyphyResultsManager()
    
    # Load results from JSON directories
    # relax_result = manager.load_relax_from_json('/path/to/relax/jsons/')
    # busted_ph_result = manager.load_busted_ph_from_json('/path/to/busted_ph/jsons/')
    # absrel_result = manager.load_absrel_from_json('/path/to/absrel/jsons/')
    
    # Filter results by omega values
    # filtered_relax = relax_result.filter_omega(10000)
    # significant_relax = relax_result.get_significant_results()
    
    # Get selection type classifications
    # selection_types = relax_result.classify_selection_type()
    # counts = relax_result.count_selection_types()
    
    # Save results for future use
    # manager.save_all_results('/path/to/save/directory/')
    
    # Compare results across analyses
    # overlap_stats = manager.get_overlap_stats(['relax', 'busted_ph'])
    
    print("Example usage completed. See function comments for actual usage.")

if __name__ == "__main__":
    example_usage()
