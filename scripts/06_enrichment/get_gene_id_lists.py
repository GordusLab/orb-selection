"""Extract gene LOCs from permulation-test results or HyPhy results."""

import os
import pickle
import importlib.util
import warnings
import pandas as pd

scripts = os.path.dirname(__file__)
repo_root = os.path.abspath(os.path.join(scripts, "..", ".."))
src_path = os.path.join(repo_root, "src")
stage04_path = os.path.join(repo_root, "scripts", "04_selection_tests")
stage05_path = os.path.join(repo_root, "scripts", "05_permulation_loss_dup")

import sys
for path in (src_path, stage04_path, stage05_path):
    if path not in sys.path:
        sys.path.insert(0, path)

from id_converter import convert_hogs_to_locs


def _load_module(module_name, module_path):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module {module_name} from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


odds_ratio_test = _load_module(
    "odds_ratio_test", os.path.join(stage05_path, "odds_ratio_test.py")
)
hyphy_results_helpers = _load_module(
    "hyphy_results_helpers", os.path.join(stage04_path, "hyphy_results_helpers.py")
)
hyphy_results_parser = _load_module(
    "hyphy_results_parser", os.path.join(stage04_path, "hyphy_results_parser.py")
)

get_fltrd_LOCs = hyphy_results_helpers.get_fltrd_LOCs
filter_omega = hyphy_results_helpers.filter_omega
HyphyResult = hyphy_results_parser.HyphyResult

warnings.filterwarnings("ignore", category=SyntaxWarning)

data_dir = os.path.join(repo_root, "data")


def load_pickle_file(fname):
    """
    Load a previously saved results pickle file.
    """

    with open(fname, "rb") as file:
        results = pickle.load(file)

    return results


def get_universe_permulation(perm_results) -> pd.Series:
    """
    Generate universe LOC files for topGO analysis at different occupancy thresholds.
    """
    min_occupancy = perm_results.occupancy_threshold
    max_occupancy = perm_results.maximum
    genecount_df = perm_results.true_odds.genecount_df
    genes_tsv = os.path.join(data_dir, "N5.tsv")

    print("Generating universe LOCs...")

    genecount_df["occupancy"] = (
        genecount_df.select_dtypes(include="number").astype("bool").sum(axis=1)
    )

    # Drop all columns except 'HOG' and 'occupancy'
    genecount_df = genecount_df[["occupancy"]]

    if max_occupancy is not None:
        genecount_df = genecount_df[
            (genecount_df["occupancy"] >= min_occupancy)
            & (genecount_df["occupancy"] <= max_occupancy)
        ]
    else:
        genecount_df = genecount_df[genecount_df["occupancy"] >= min_occupancy]

    LOCs_df = convert_hogs_to_locs(genecount_df, genes_tsv)

    universe_loc_ids = LOCs_df["LOC"].dropna().unique()

    return universe_loc_ids


def main_permulation(perm_results, tail=None):
    """
    Main function to load permulation-test results and extract LOCs.
    """

    print("Extracting LOCs for hits...")
    df = perm_results.results_fltrd_df
    genes_tsv = perm_results.true_odds.hog_node_genes_tsv

    df = convert_hogs_to_locs(df, genes_tsv)

    if tail == "left":
        df = df[df["Log odds ratio"] < 0]
    elif tail == "right":
        df = df[df["Log odds ratio"] > 0]

    hit_loc_ids = df["LOC"].dropna().unique()

    return hit_loc_ids


def main_hyphy(hyphy_results, omega=10000, relax_result=None):
    """
    Main function to load the HyPhy results and extract LOCs.
    """

    hyphy_df = hyphy_results.results_df
    test = hyphy_results.analysis_type

    genes_tsv = os.path.join(data_dir, "N5.tsv")
    print("Extracting LOCs...")
    merged_df = convert_hogs_to_locs(hyphy_df, genes_tsv)

    hit_loc_ids = get_fltrd_LOCs(
        merged_df, omega=omega, test=test, relax_result=relax_result
    )

    universe_fltrd = filter_omega(merged_df, omega)
    universe_loc_ids = universe_fltrd["LOC"].dropna().unique()

    return hit_loc_ids, universe_loc_ids


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Extract LOCs from permulation-test or HyPhy results pickle file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python scripts/06_enrichment/get_gene_id_lists.py results.pkl
    python scripts/06_enrichment/get_gene_id_lists.py results.pkl --tail left
    python scripts/06_enrichment/get_gene_id_lists.py results.pkl --tail right --hits-file locs.txt
        """,
    )

    parser.add_argument("pickle_file", help="Path to the pickle file (.pkl)")
    parser.add_argument(
        "--tail",
        choices=["left", "right"],
        help="Filter permulation-test results by tail direction (left or right)",
    )
    parser.add_argument(
        "--omega",
        type=float,
        default=10000,
        help="Omega value for filtering HyPhy results (default: 10000)",
    )
    parser.add_argument(
        "--relax-result",
        choices=["relaxed", "intensified"],
        help="Relaxation result type for filtering (default: None)",
    )
    parser.add_argument(
        "--hits-file",
        dest="hits_file",
        help="Save results to file instead of printing to console",
    )
    parser.add_argument(
        "--universe-file", dest="universe_file", help="Save universe LOCs to file"
    )

    args = parser.parse_args()

    pickle_file_path = args.pickle_file
    tail_arg = args.tail
    omega_arg = args.omega
    relax_result_arg = args.relax_result
    hits_file = args.hits_file
    universe_file = args.universe_file

    universe_locs = []

    if not pickle_file_path.endswith(".pkl"):
        print("Error: The file must be a pickle file with .pkl extension.")
        sys.exit(1)
    if not os.path.exists(pickle_file_path):
        print(f"Error: The file {pickle_file_path} does not exist.")
        sys.exit(1)

    print(f"Loading results from {pickle_file_path}...")
    loaded_results = load_pickle_file(pickle_file_path)

    if isinstance(loaded_results, odds_ratio_test.PermutationTestResults):
        if universe_file:
            universe_locs = get_universe_permulation(loaded_results)
        hit_locs = main_permulation(loaded_results, tail=tail_arg)
    elif isinstance(loaded_results, HyphyResult):
        hit_locs, universe_locs = main_hyphy(
            loaded_results, omega=omega_arg, relax_result=relax_result_arg
        )
    else:
        print(
            "Error: Unsupported results type. Expected PermutationTestResults or HyphyResult."
        )
        sys.exit(1)

    print(f"Extracted {len(hit_locs)} unique significant LOCs.")

    if hits_file:
        # Save to file
        try:
            with open(hits_file, "w", encoding="utf-8") as f:
                for loc in hit_locs:
                    f.write(f"{loc}\n")
            print(f"Results saved to {hits_file}")
        except IOError as e:
            print(f"Error writing to file {hits_file}: {e}")
            sys.exit(1)

    if universe_file:
        try:
            with open(universe_file, "w", encoding="utf-8") as f:
                for loc in universe_locs:
                    f.write(f"{loc}\n")
            print(f"Universe LOCs saved to {universe_file}")
        except IOError as e:
            print(f"Error writing to file {universe_file}: {e}")
            sys.exit(1)
