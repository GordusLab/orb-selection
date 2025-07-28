import os
import pickle
import warnings
warnings.filterwarnings('ignore', category=SyntaxWarning)
import odds_ratio_test
import id_converter
from hyphy_results_helpers import get_fltrd_LOCs, filter_omega
from hyphy_results_parser import HyphyResult
import pandas as pd

scripts = os.path.dirname(__file__)
data = os.path.join(scripts, "..", "data")


def get_universe_bootstrap(genecount_df: pd.DataFrame, min_occupancy: int = 0, max_occupancy: int = None) -> pd.Series:
    """
    Generate universe LOC files for topGO analysis at different occupancy thresholds.

    Args:
        NX_genecount_tsv: Path to gene count TSV file
        min_occupancy: Minimum occupancy threshold
        max_occupancy: Maximum occupancy threshold (optional, defaults to None)
    """
    genes_tsv = f"{data}/N5.tsv"

    LOCs_df = id_converter.main(genecount_df, genes_tsv)
    LOCs_df = LOCs_df.set_index('HOG')
    LOCs_df['occupancy'] = LOCs_df.select_dtypes(include='number').astype('bool').sum(axis=1)

    LOCs_fltrd = LOCs_df[LOCs_df['occupancy'] >= min_occupancy and (LOCs_df['occupancy'] <= max_occupancy if max_occupancy is not None else True)]
    universe_locs = LOCs_fltrd['LOC'].dropna().unique()

    return universe_locs

def get_universe_hyphy(hyphy_results: HyphyResult, omega: int):
    """
    Generate universe LOCs for HyPhy results based on omega threshold.

    Args:
        hyphy_results: HyPhyResult object containing results DataFrame
        omega: Maximum omega threshold
    """
    genes_tsv = f"{data}/N5.tsv"
    
    LOCs_df = id_converter.main(hyphy_results.results_df, genes_tsv)
    LOCs_df = LOCs_df.set_index('HOG')

    LOCs_fltrd = filter_omega(LOCs_df, omega)
    universe_locs = LOCs_fltrd['LOC'].dropna().unique()

    return universe_locs

def load_pickle_file(fname):
    """
    Loads a previously saved the bootstrap results pickle file.
    """

    with open(fname, 'rb') as file:
        results = pickle.load(file)

    return results

def main_bs(bs_results, tail=None):
    """
    Main function to load the bootstrap results and extract LOCs.
    """

    print("Extracting LOCs...")
    df = bs_results.results_fltrd_df
    genes_tsv = bs_results.true_odds.hog_node_genes_tsv

    df = id_converter.main(df, genes_tsv)

    if tail == 'left':
        df = df[df['Log odds ratio'] < 0]
    elif tail == 'right':
        df = df[df['Log odds ratio'] > 0]

    locs_list = df["LOC"].dropna().unique()

    return locs_list

def main_hyphy(hyphy_results, omega=10000, relax_result=None):
    """
    Main function to load the HyPhy results and extract LOCs.
    """

    hyphy_df = hyphy_results.results_df
    test = hyphy_results.analysis_type

    genes_tsv = f"{data}/N5.tsv"
    print("Extracting LOCs...")
    merged_df = id_converter.main(hyphy_df, genes_tsv)

    loc_list = get_fltrd_LOCs(merged_df, omega=omega, test=test, relax_result=relax_result)

    return loc_list

if __name__ == "__main__":
     
    import sys
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Extract LOCs from bootstrap results pickle file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python script.py results.pkl
  python script.py results.pkl --tail left
  python script.py results.pkl --tail right --output locs.txt
  python script.py results.pkl --save-to locs.txt
        """
    )
    
    parser.add_argument('pickle_file', help='Path to the pickle file (.pkl)')
    parser.add_argument('--tail', choices=['left', 'right'], 
                       help='Filter bootstrap results by tail direction (left or right)')
    parser.add_argument('--omega', type=float, default=10000, 
                        help='Omega value for filtering HyPhy results (default: 10000)')
    parser.add_argument('--relax-result', choices=['relaxed', 'intensified'], 
                        help='Relaxation result type for filtering (default: None)')
    parser.add_argument('--output', '--save-to', dest='hits_file', 
                    help='Save results to file instead of printing to console')
    parser.add_argument('--universe-file', dest='universe_file', default=None,
                        help='Save universe LOCs to file')

    args = parser.parse_args()
    
    fname = args.pickle_file
    tail = args.tail
    omega = args.omega
    relax_result = args.relax_result
    hits_file = args.hits_file
    universe_file = args.universe_file

    if not fname.endswith('.pkl'):
        print("Error: The file must be a pickle file with .pkl extension.")
        sys.exit(1)
    if not os.path.exists(fname):
        print(f"Error: The file {fname} does not exist.")
        sys.exit(1)

    print(f"Loading results from {fname}...")
    results = load_pickle_file(fname)

    if isinstance(results, odds_ratio_test.BootstrapTestResults):
        locs_list = main_bs(results, tail=tail)
        universe_file = get_universe_bootstrap(results.true_odds.genecount_df, results.occupancy_threshold, results.maximum)
    elif isinstance(results, HyphyResult):
        locs_list = main_hyphy(results, omega=omega, relax_result=relax_result)
        universe_file = get_universe_hyphy(results, omega=omega)
    else:
        print("Error: Unsupported results type. Expected BootstrapTestResults or HyphyResult.")
        sys.exit(1) 

    print(f"Extracted {len(locs_list)} unique LOCs.")

    if hits_file:
        # Save to file
        try:
            with open(hits_file, 'w') as f:
                for loc in locs_list:
                    f.write(f"{loc}\n")
            print(f"Results saved to {hits_file}")
        except IOError as e:
            print(f"Error writing to file {hits_file}: {e}")
            sys.exit(1)

    if universe_file:
        try:
            with open(universe_file, 'w') as f:
                for loc in universe_file:
                    f.write(f"{loc}\n")
            print(f"Universe LOCs saved to {universe_file}")
        except IOError as e:
            print(f"Error writing to file {universe_file}: {e}")
            sys.exit(1)
