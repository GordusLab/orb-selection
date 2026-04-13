"""
Odds ratio analysis of gene loss and gain: testing whether
the odds that genes in an orthogroup are lost or duplicated
differs significantly between a test group and a reference
group in a phylogeny

Follow the instructions at https://github.com/nclark-lab/RERconverge/wiki/Install
to set up RERconverge and run the src/permulations.R script before proceeding.
"""

import sys
import os
import random
import pickle
from datetime import datetime
from typing import Dict, List, Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from tqdm.auto import tqdm
import seaborn as sns
import orthogroup_gene_count
import id_converter

# Set the random seed for reproducibility
random.seed(42)

plt.rcParams["font.family"] = "Verdana"

_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_SRC_DIR, os.pardir, os.pardir))
CONSOLE_PRINT_WIDTH = 96


def _fixed_width_lines(text: str, width: int = CONSOLE_PRINT_WIDTH) -> List[str]:
    """Return fixed-width wrapped lines for aligned console printing."""

    def _wrap_no_word_breaks(line: str) -> List[str]:
        if len(line) <= width:
            return [line]

        wrapped: List[str] = []
        remaining = line
        break_chars = (" ", "\t", "/", "_", "-", ",", ";")

        while len(remaining) > width:
            window = remaining[: width + 1]
            break_at = max(window.rfind(ch) for ch in break_chars)

            if break_at <= 0:
                # No useful breakpoints; hard-wrap at fixed width.
                wrapped.append(remaining[:width])
                remaining = remaining[width:]
                continue

            split_idx = break_at + 1
            wrapped.append(remaining[:split_idx].rstrip())

            # Drop leading whitespace only when splitting on whitespace.
            if remaining[break_at] in {" ", "\t"}:
                remaining = remaining[split_idx:].lstrip()
            else:
                remaining = remaining[split_idx:]

        wrapped.append(remaining)
        return wrapped

    wrapped_lines: List[str] = []
    for raw_line in str(text).splitlines() or [""]:
        segments = _wrap_no_word_breaks(raw_line)
        if not segments:
            wrapped_lines.append("")
            continue
        wrapped_lines.extend(segments)
    return wrapped_lines


def _cprint(text: str = "", width: int = CONSOLE_PRINT_WIDTH) -> None:
    """Print fixed-width lines to stdout."""
    for line in _fixed_width_lines(text, width=width):
        print(line)


def _emit(text: str, file_obj, width: int = CONSOLE_PRINT_WIDTH) -> None:
    """Write output to file object; fixed-width only when writing to stdout."""
    if file_obj is sys.stdout:
        _cprint(text, width=width)
    else:
        print(text, file=file_obj)


def _resolve_repo_path(path: Optional[str]) -> Optional[str]:
    """Resolve user paths with repo-relative defaults and external path support.

    Resolution order:
    1. Absolute paths and '~' expanded paths are used as-is.
    2. Explicit './' or '../' paths are resolved from current working directory.
    3. Bare relative paths are resolved from repository root.
    """
    if path is None:
        return None

    expanded = os.path.expanduser(path)
    if os.path.isabs(expanded):
        return expanded

    if (
        expanded.startswith(f".{os.sep}")
        or expanded == "."
        or expanded.startswith(f"..{os.sep}")
        or expanded == ".."
    ):
        return os.path.abspath(expanded)

    return os.path.abspath(os.path.join(_REPO_ROOT, expanded))


def _unique_results_dir(
    parent_dir: str,
    time_obj: datetime,
    min_occ: int,
    max_occ: Optional[int],
    permulation_reps: int,
) -> str:
    """
    Create a dated parent folder, then return a unique run subdirectory inside it
    named with sequential run number and test parameters.

    Example: If parent_results_dir is "results/", this creates:
      results/Results_Jan29/Run1_occ_0-100_10000x/
      results/Results_Jan29/Run2_occ_50-75_10000x/ (different params, Run2)
      results/Results_Jan29/Run3_occ_0-100_10000x/ (same params as Run1, Run3)
    """
    # Create the dated parent folder (e.g., Results_Jan29)
    date_short = time_obj.strftime("%b%d")
    dated_parent = f"{parent_dir}/Results_{date_short}"
    os.makedirs(dated_parent, exist_ok=True)

    # Build the directory name with test parameters
    # Format occupancy range
    occ_range = f"{min_occ}-{max_occ}"

    # Build base directory name
    base_name = f"occ_{occ_range}_{permulation_reps}x"

    # Find the highest run number across all existing directories
    existing_run_nums = []
    if os.path.exists(dated_parent):
        for entry in os.listdir(dated_parent):
            if entry.startswith("Run"):
                # Extract run number from "Run1_...", "Run2_...", etc.
                try:
                    run_num = int(entry.split("_")[0][3:])
                    existing_run_nums.append(run_num)
                except (ValueError, IndexError):
                    pass

    # Determine next run number
    next_run_num = max(existing_run_nums) + 1 if existing_run_nums else 1

    subdir_name = f"Run{next_run_num}_{base_name}"
    subdir_path = f"{dated_parent}/{subdir_name}"

    return subdir_path


def drop_empty_cols(df, print_txt=True):
    """
    Drops columns where all entries (ignoring headers) are 0,
    to get rid of species not included in the current node/hierarchy.
    """

    # Print the number of columns before cleaning
    num_columns_before = df.shape[1]
    if print_txt:
        _cprint(
            f"Number of columns before dropping empty columns: {num_columns_before}"
        )

    # Drop columns where all entries (ignoring headers) are 0
    df_cleaned = df.loc[:, (df.ne(0)).any(axis=0)]

    # Print the number of columns after cleaning
    num_columns_after = df_cleaned.shape[1]

    if print_txt:
        _cprint(f"Number of columns after dropping empty columns: {num_columns_after}")
        _cprint("Species with no sequences in any orthogroup have been dropped.")

    return df_cleaned


def occupancy_filter(arr, min_occ, max_occ, total_occ_arr):
    """
    Function to filter an array of values according to whether the
    orthogroup meets a certain occupancy threshold.

    Args:
        arr, min_occ, max_occ, total_occ_arr

    Returns:
        fltrd_arr (numpy array): occupancy array filtered according
        to thresholds provided
    """

    if max_occ is None:
        max_occ = total_occ_arr.max()

    idx = np.asarray((total_occ_arr >= min_occ) & (total_occ_arr <= max_occ)).nonzero()[
        0
    ]
    fltrd_arr = arr[idx]

    return fltrd_arr

def filter_for_sp_of_interest(df, genecount_df, species_name):
    """Filter the DataFrame for HOGs that include a species of interest"""

    sp_of_interest_present = genecount_df[genecount_df[species_name] != 0]

    sp_of_interest_present_hogs = sp_of_interest_present.index.values
    sp_of_interest_present_hogs = set(sp_of_interest_present_hogs)
    df_fltrd_sp_of_int = df[df.index.isin(sp_of_interest_present_hogs)]

    return len(df_fltrd_sp_of_int)


def calculate_odds(foreground_bool_arr, background_bool_arr, test_bool_mat, busco_arr):
    """Function to calculate the odds ratio and log odds ratio"""

    # I don't know why this is necessary but my kernel crashes without it
    foreground_bool_arr = foreground_bool_arr.reshape(foreground_bool_arr.size, 1)
    background_bool_arr = background_bool_arr.reshape(background_bool_arr.size, 1)

    # If a busco array is provided, weight the foreground and background arrays
    if busco_arr is not None:
        foreground_bool_arr = foreground_bool_arr.flatten() * busco_arr
        background_bool_arr = background_bool_arr.flatten() * busco_arr

        foreground_bool_arr = foreground_bool_arr.reshape(foreground_bool_arr.size, 1)
        background_bool_arr = background_bool_arr.reshape(background_bool_arr.size, 1)

    # Calculate the number of foreground and background that are
    # [missing, duplicated]

    ## foreground yes
    foreground_yes_arr = np.matmul(test_bool_mat, foreground_bool_arr)

    ## background yes
    background_yes_arr = np.matmul(test_bool_mat, background_bool_arr)

    # Calculate the number of foreground and background that are
    # [not missing, not duplicated]
    test_inv_bool_mat = 1 - test_bool_mat

    ## foreground no
    foreground_no_arr = np.matmul(test_inv_bool_mat, foreground_bool_arr)

    ## background no
    background_no_arr = np.matmul(test_inv_bool_mat, background_bool_arr)

    # Use the Haldane-Anscombe correction to account for any 0 entries
    # when calculating the odds ratios
    foreground_yes_arr += 0.5
    background_yes_arr += 0.5
    foreground_no_arr += 0.5
    background_no_arr += 0.5

    # Calculate odds missing for foreground & background, odds ratio, and log odds ratio

    ## odds
    odds_foreground_arr = foreground_yes_arr / foreground_no_arr
    odds_background_arr = background_yes_arr / background_no_arr

    ## odds ratio
    odds_ratio_arr = odds_foreground_arr / odds_background_arr

    ## log odds ratio
    log_odds_ratio_arr = np.log(odds_ratio_arr)

    return log_odds_ratio_arr


class OddsRatioResults:
    """Class to hold the results of the odds ratio calculations"""

    def __init__(
        self,
        genecount_csv,
        hog_node_genes_tsv,
        time,
        foreground_list_filename,
        background_list_filename=None,
        buscos_filename=None,
    ):
        """Initialize the inputs for the odds ratio calculations"""
        self.genecount_csv = genecount_csv
        self.hog_node_genes_tsv = hog_node_genes_tsv
        self.time = time
        self.foreground_list_filename = foreground_list_filename
        self.background_list_filename = background_list_filename

        self.foreground_list_arr = np.loadtxt(foreground_list_filename, dtype=str)

        if background_list_filename is not None:
            self.background_list_arr = np.loadtxt(background_list_filename, dtype=str)
        else:
            self.background_list_arr = None

        if buscos_filename is not None:
            self.buscos = pd.read_csv(buscos_filename, header=None)
            self.buscos.columns = [
                "Species",
                "Single_copy_buscos",
                "Duplicated_buscos",
                "Total_buscos",
                "Fraction_sc",
                "Fraction_total",
            ]

        # Get the genecount arrays and matrices
        self._get_genecount_arrays()

        # Define the foreground and background arrays
        self._define_foreground()

        # Calculate the log odds ratios of loss
        self.loss_lor_arr = calculate_odds(
            self.foreground_bool_arr,
            self.background_bool_arr,
            self.loss_bool_mat,
            busco_arr=getattr(self, "loss_busco_arr", None),
        )

        # Calculate the log odds ratios of duplication
        self.dup_lor_arr = calculate_odds(
            self.foreground_bool_arr,
            self.background_bool_arr,
            self.dup_bool_mat,
            busco_arr=getattr(self, "dup_busco_arr", None),
        )

        # Create the results DataFrame
        self.results_df = self.results_to_df()

    def __repr__(self):
        """String representation of the OddsRatioResults object"""
        return (
            f"OddsRatioResults(genecount_csv={self.genecount_csv}, "
            f"hog_node_genes_tsv={self.hog_node_genes_tsv}, "
            f"foreground_list_filename={self.foreground_list_filename}, "
            f"background_list_filename={self.background_list_filename})"
        )

    def print_attributes(self):
        """Print all instance attribute names and types."""
        for attr_name in sorted(vars(self)):
            attr_value = getattr(self, attr_name)
            _cprint(f"{attr_name}: {type(attr_value).__name__}")

    def _get_genecount_arrays(self):
        """Function to create table of counts of genes per species in each orthogroup
        at the specified node of my family tree, and matrices of 1s and 0s
        corresponding to whether a gene is missing, duplicated, or single copy
        for each species."""

        # Counts of genes per species in each HOG
        genecount_df = pd.read_csv(self.genecount_csv, index_col="HOG", sep="\t")
        genecount_df = genecount_df.drop(
            columns=["OG", "Gene Tree Parent Clade"], errors="ignore"
        )

        # If there is a 'Total' column, drop it
        if "Total" in genecount_df.columns:
            genecount_df = genecount_df.drop(columns=["Total"])

        # If there is an 'Occupancy' column, drop it
        if "Occupancy" in genecount_df.columns:
            genecount_df = genecount_df.drop(columns=["Occupancy"])

        # Remove any empty columns (species with no genes in any HOGs)
        genecount_df = drop_empty_cols(genecount_df, print_txt=False)

        # Remove any species not in the foreground or background lists
        if self.background_list_arr is not None:
            species_to_keep = self.foreground_list_arr

            species_to_keep = np.concatenate(
                (species_to_keep, self.background_list_arr), axis=0
            )
            genecount_df = genecount_df[species_to_keep]

        # List of HOG IDs
        hog_list = genecount_df.index.values

        # List of all species
        all_species_arr = genecount_df.columns.to_numpy()

        # Get busco scores for each species
        if hasattr(self, "buscos"):
            sc_busco_arr = (
                self.buscos.set_index("Species")
                .loc[all_species_arr]["Fraction_sc"]
                .values
            )
            total_busco_arr = (
                self.buscos.set_index("Species")
                .loc[all_species_arr]["Fraction_total"]
                .values
            )

        # Convert counts df to numpy
        genecount_mat = genecount_df.to_numpy()

        # Matrix where 1=species is present in the HOG, 0=species has no gene in the HOG
        occ_bool_mat = genecount_mat.astype(bool).astype(float)

        # Calculate the total occupancy of each HOG (number of species present)
        total_occ_arr = occ_bool_mat.sum(axis=1).reshape(genecount_mat.shape[0], 1)

        # Matrix where 1=species is absent from the HOG, 0=species is present
        loss_bool_mat = np.asarray(genecount_mat == 0).astype(float)

        # Matrix where 1=gene with more than 1 copy, 0=gene with 1 or no copies
        dup_bool_mat = np.asarray(genecount_mat > 1).astype(float)

        self.genecount_df = genecount_df
        self.hog_list = hog_list
        self.all_species_arr = all_species_arr
        self.total_occ_arr = total_occ_arr
        self.loss_bool_mat = loss_bool_mat
        self.dup_bool_mat = dup_bool_mat

        # When testing for loss, the total fraction of complete buscos, whether single-copy or
        # duplicated, should be a good proxy for how reliable loss calls will be
        if hasattr(self, "buscos"):
            self.loss_busco_arr = total_busco_arr

        # When testing for duplication, the fraction of recovered buscos that are single-copy
        # should be a good proxy for how reliable duplication calls will be
        if hasattr(self, "buscos"):
            self.dup_busco_arr = sc_busco_arr

    def _define_foreground(self):
        """Function to make a numpy vector corresponding to whether
        each species is foreground or not: 1 = foreground, 0 = background.
        If the background are not specified, any animal that is
        not specified to be OW is considered non-OW."""

        # assign True to index of foreground in the species array
        foreground_bool_arr = np.isin(self.all_species_arr, self.foreground_list_arr)

        total_species_count = len(foreground_bool_arr)
        foreground_count = foreground_bool_arr.sum()

        # If no background is specified, all species not identified as
        # foreground are considered background
        if self.background_list_arr is None:
            background_bool_arr = np.isin(
                self.all_species_arr, self.foreground_list_arr, invert=True
            )

        # If the background are specified
        else:
            background_bool_arr = np.isin(
                self.all_species_arr, self.background_list_arr
            )

        background_count = background_bool_arr.sum()

        background_bool_arr = background_bool_arr.reshape(
            background_bool_arr.size, 1
        ).astype(float)
        foreground_bool_arr = foreground_bool_arr.reshape(
            foreground_bool_arr.size, 1
        ).astype(float)

        _cprint(
            f"{total_species_count} species total, {foreground_count} foreground, "
            f"{background_count} background"
        )

        self.total_species_count = total_species_count
        self.foreground_count = foreground_count
        self.background_count = background_count
        self.foreground_bool_arr = foreground_bool_arr
        self.background_bool_arr = background_bool_arr

    def results_to_df(self):
        """Function to convert the results of the odds ratio calculations
        into a DataFrame for easier analysis and plotting"""

        # Create a DataFrame with the results
        results_df = pd.DataFrame(
            {
                "HOG": self.hog_list,
                "Occupancy": self.total_occ_arr.flatten().astype(int),
                "Log odds ratio of loss": self.loss_lor_arr.flatten(),
                "Log odds ratio of duplication": self.dup_lor_arr.flatten(),
            }
        )
        results_df = results_df.set_index("HOG")

        return results_df


def load_permulation_tip_values_from_csv(csv_path: str) -> List[Dict[str, float]]:
    """Load permulated tip values from a CSV file.

    Expected format:
    - One row per permulation.
    - Species names as columns.
    - Optional first column named `perm_id`.
    - Cell values are numeric (typically 0/1) tip assignments.
    """

    csv_path = _resolve_repo_path(csv_path)
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Permulation tip-values CSV not found: {csv_path}")

    tip_df = pd.read_csv(csv_path)
    if tip_df.empty:
        raise ValueError(f"Permulation tip-values CSV is empty: {csv_path}")

    species_cols = [c for c in tip_df.columns if c != "perm_id"]
    if not species_cols:
        raise ValueError(
            f"Permulation tip-values CSV has no species columns: {csv_path}"
        )

    tip_numeric = tip_df[species_cols].apply(pd.to_numeric, errors="raise")
    if tip_numeric.isnull().any().any():
        raise ValueError(f"Permulation tip-values CSV contains NaN values: {csv_path}")

    return [
        {species: float(value) for species, value in row.items()}
        for row in tip_numeric.to_dict(orient="records")
    ]

def save_loc_list(df, NX_path, output_path):
    locs_df = id_converter.convert_hogs_to_locs(df, NX_path, show_progress=False)
    list = locs_df['LOC'].dropna().unique()
    with open(output_path, 'w', encoding="utf-8") as f:
        for loc in list:
            f.write(f"{loc}\n")
    print(f"Wrote {len(list)} items to {output_path}")
    return list

class PermulationTestResults:
    """Class to perform the permulation test for the odds ratio"""

    def __init__(
        self,
        true_odds,
        permulation_tip_values,
        min_occ=0,
        max_occ=None,
        alpha=0.05,
        permulation_reps=10000,
        permulation_rows_available: Optional[int] = None,
        sampled_from_available: bool = False,
        species_of_interest=None,
    ):
        """Initialize the permulation test for a given odds ratio results object"""

        self.true_odds = true_odds

        self.min_occ = min_occ

        self.permulation_tip_values = permulation_tip_values
        self.permulation_rows_available = permulation_rows_available
        self.sampled_from_available = sampled_from_available

        if self.permulation_tip_values is not None:
            if len(self.permulation_tip_values) == 0:
                raise ValueError("permulation_tip_values cannot be empty.")
            self.permulation_reps = len(self.permulation_tip_values)
        else:
            self.permulation_reps = permulation_reps

        self.species_of_interest = species_of_interest

        if max_occ is not None:
            self.max_occ = max_occ
        else:
            self.max_occ = self.true_odds.total_species_count

        self.alpha = alpha
        self.loss_ci_av = None
        self.dup_ci_av = None

        # Filter the true log odds ratios based on the occupancy threshold
        self.true_fltrd_loss_lors = occupancy_filter(
            self.true_odds.loss_lor_arr,
            self.min_occ,
            self.max_occ,
            self.true_odds.total_occ_arr,
        )

        self.true_fltrd_dup_lors = occupancy_filter(
            self.true_odds.dup_lor_arr,
            self.min_occ,
            None,
            self.true_odds.total_occ_arr,
        )

        # Calculate the mean and standard deviation of the true log odds ratios
        self.true_mean_loss = np.mean(self.true_fltrd_loss_lors)
        self.true_stddev_loss = np.std(self.true_fltrd_loss_lors)
        self.true_mean_dup = np.mean(self.true_fltrd_dup_lors)
        self.true_stddev_dup = np.std(self.true_fltrd_dup_lors)

        _cprint(
            f"Log odds ratios of LOSS: \n Mean: {self._fmt_stat(self.true_mean_loss)}, Stddev: {self._fmt_stat(self.true_stddev_loss)}, Count of HOGs: {len(self.true_fltrd_loss_lors)}"
        )
        _cprint(
            f"Log odds ratios of DUPLICATION: \n Mean: {self._fmt_stat(self.true_mean_dup)}, Stddev: {self._fmt_stat(self.true_stddev_dup)}, Count of HOGs: {len(self.true_fltrd_dup_lors)}"
        )

        # Run the permulation test
        self._run_permulation()

        # Get the list of HOGs which are significantly different between
        # the test groups according to the permulation thresholds
        self._get_hits_dfs()

        # Print the results of the permulation test
        self.print_permulation_results()

    def __repr__(self):
        """String representation of the PermulationTestResults object"""
        return (
            f"PermulationTestResults(true_odds={self.true_odds}, "
            #            f"stat={self.stat}, "
            f"min_occ={self.min_occ}, "
            f"permulation_reps={self.permulation_reps})"
        )

    def print_attributes(self):
        """Print all instance attribute names and types."""
        for attr_name in sorted(vars(self)):
            attr_value = getattr(self, attr_name)
            _cprint(f"{attr_name}: {type(attr_value).__name__}")

    def get_pval_thresholds(self, alpha, alternative):
        """Calculate the permulation p-value thresholds for the log odds ratio"""

        if alternative == "two-tailed":
            z_crit = norm.ppf(1 - alpha / 2)  # two-tailed z critical value
        else:
            z_crit = norm.ppf(1 - alpha)  # one-tailed z critical value

        def compute_ci(means, stddevs):
            lowers = means - z_crit * stddevs
            uppers = means + z_crit * stddevs

            # Shape: (n_permulations, 2) where each row is [lower, upper].
            cis = np.column_stack((lowers, uppers))

            # Average CI as [mean_lower, mean_upper].
            ci_av = np.array(
                [
                    np.mean(cis[:, 0]),
                    np.mean(cis[:, 1]),
                ]
            )
            return ci_av

        self.loss_ci_av = compute_ci(self.means_loss, self.stddevs_loss)
        self.dup_ci_av = compute_ci(self.means_dup, self.stddevs_dup)

        return self.loss_ci_av, self.dup_ci_av
    
    def get_universe_permulation(self) -> tuple:
        """
        Generate universe LOC files for topGO analysis at different occupancy thresholds.
        """
        min_occupancy = self.min_occ
        max_occupancy = self.max_occ
        genecount_df = self.true_odds.genecount_df

        print("Generating universe LOCs...")

        genecount_df["occupancy"] = (
            genecount_df.select_dtypes(include="number").astype("bool").sum(axis=1)
        )

        # Drop all columns except 'HOG' and 'occupancy'
        genecount_df = genecount_df[["occupancy"]]

        # Filter for minimum occupancy    
        genecount_df = genecount_df[genecount_df["occupancy"] >= min_occupancy]

        loss_uni_df = genecount_df[
            (genecount_df["occupancy"] >= min_occupancy)
            & (genecount_df["occupancy"] <= max_occupancy)
        ]
        dup_uni_df = genecount_df.copy()

        return loss_uni_df, dup_uni_df

    def save_go_lists(self, results_dir: str, use_perm_pvals=False):
        
        dfs = self.results_fltrd_dfs
        if use_perm_pvals:
            dfs["loss_fg_perm_pval"] = dfs["loss_fg"][dfs["loss_fg"]["Significant in permulation test"] == True]
            dfs["loss_bg_perm_pval"] = dfs["loss_bg"][dfs["loss_bg"]["Significant in permulation test"] == True]
            dfs["dup_fg_perm_pval"] = dfs["dup_fg"][dfs["dup_fg"]["Significant in permulation test"] == True]
            dfs["dup_bg_perm_pval"] = dfs["dup_bg"][dfs["dup_bg"]["Significant in permulation test"] == True]

        # Convert HOG hit list to LOCs + descriptions and save as a companion file.
        for key in dfs.keys():
            save_loc_list(self.results_fltrd_dfs[key], self.true_odds.hog_node_genes_tsv, f"{results_dir}/{key}_sig_locs.txt")
        
        loss_uni, dup_uni = self.get_universe_permulation()
        
        save_loc_list(loss_uni, self.true_odds.hog_node_genes_tsv, f"{results_dir}/loss_universe_locs.txt")
        save_loc_list(dup_uni, self.true_odds.hog_node_genes_tsv, f"{results_dir}/dup_universe_locs.txt")


    @staticmethod
    def _fmt_stat(value: float, ndigits: int = 2) -> str:
        """Format a float while normalizing negative zero to positive zero."""
        rounded = round(float(value), ndigits)
        if rounded == 0.0:
            rounded = 0.0
        return f"{rounded:.{ndigits}f}"

    @staticmethod
    def _fmt_ci(ci_values, ndigits: int = 3) -> str:
        """Format CI vectors like [lower, upper] with fixed precision."""
        ci_arr = np.asarray(ci_values, dtype=float).flatten()
        if ci_arr.size == 0:
            return "[]"
        formatted = ", ".join(f"{x:.{ndigits}f}" for x in ci_arr)
        return f"[{formatted}]"

    def _foreground_background_from_permulation(
        self,
        perm_tips: Dict[str, float],
        species_incl_idx: np.ndarray,
    ):
        """Map a permulated named tip-value vector onto analysis species order."""

        new_foregrounds = np.zeros(self.true_odds.total_species_count)

        mapped = 0
        for idx in species_incl_idx:
            species = self.true_odds.all_species_arr[idx]
            if species in perm_tips:
                new_foregrounds[idx] = 1 if float(perm_tips[species]) == 1.0 else 0
                mapped += 1

        if mapped == 0:
            raise ValueError(
                "No species names from permulation tips matched analysis species names."
            )

        new_backgrounds = np.zeros(self.true_odds.total_species_count)
        new_backgrounds[species_incl_idx] = 1.0 - new_foregrounds[species_incl_idx]

        return new_foregrounds, new_backgrounds

    def _permulation_iter(
        self,
        i,
        counter_mean_loss,
        counter_mean_dup,
        counters_hogs,
        species_incl_idx,
        perm_tips=None,
    ):
        """Function to perform a single iteration of the permulation test"""

        # Create a new set of foreground and background arrays for permulation.
        new_foregrounds, new_backgrounds = self._foreground_background_from_permulation(
            perm_tips,
            species_incl_idx,
        )

        # Recalculate the odds ratios and log odds ratios for the new foreground/background arrays
        # for loss
        new_log_odds_ratio_loss = calculate_odds(
            new_foregrounds,
            new_backgrounds,
            self.true_odds.loss_bool_mat,
            busco_arr=getattr(self.true_odds, "loss_busco_arr", None),
        )

        # for duplication
        new_log_odds_ratio_dup = calculate_odds(
            new_foregrounds,
            new_backgrounds,
            self.true_odds.dup_bool_mat,
            busco_arr=getattr(self.true_odds, "dup_busco_arr", None),
        )

        # Compare the log odds ratio of each gene in the new permulated distribution
        # to the log odds ratio of the same gene in the true distribution, and count
        # how many times the permulated log odds ratio is more extreme than the true
        # log odds ratio according to the alternative hypothesis

        # for loss, greater
        counters_hogs[:, 0] += (
            new_log_odds_ratio_loss > self.true_odds.loss_lor_arr
        ).flatten()
        # for loss, less
        counters_hogs[:, 1] += (
            new_log_odds_ratio_loss < self.true_odds.loss_lor_arr
        ).flatten()

        # for duplication, greater
        counters_hogs[:, 2] += (
            new_log_odds_ratio_dup > self.true_odds.dup_lor_arr
        ).flatten()
        # for duplication, less
        counters_hogs[:, 3] += (
            new_log_odds_ratio_dup < self.true_odds.dup_lor_arr
        ).flatten()

        # Filter the new log odds ratio arrays based on the occupancy threshold
        new_loss_lor_fltrd = occupancy_filter(
            new_log_odds_ratio_loss,
            self.min_occ,
            self.max_occ,
            self.true_odds.total_occ_arr,
        )
        new_dup_lor_fltrd = occupancy_filter(
            new_log_odds_ratio_dup,
            self.min_occ,
            None,
            self.true_odds.total_occ_arr,
        )
        # Calculate the statistics for the new filtered, permulated distribution
        new_mean_loss = np.mean(new_loss_lor_fltrd)
        new_stddev_loss = np.std(new_loss_lor_fltrd)
        new_mean_dup = np.mean(new_dup_lor_fltrd)
        new_stddev_dup = np.std(new_dup_lor_fltrd)

        # Compare permulated mean to true mean according to the alternative hypothesis and update counter
        if new_mean_loss < self.true_mean_loss:
            counter_mean_loss += 1
        if new_mean_dup > self.true_mean_dup:
            counter_mean_dup += 1

        self.means_loss[i] = new_mean_loss
        self.stddevs_loss[i] = new_stddev_loss
        self.means_dup[i] = new_mean_dup
        self.stddevs_dup[i] = new_stddev_dup

        return counter_mean_loss, counter_mean_dup, counters_hogs

    def _run_permulation(self):
        """Creating 10,000 test log odds ratio distributions with the set
        of species defined as "foreground" chosen via permulations (but still the
        same number of foreground  as the true number)"""

        _cprint("\nLAUNCHING PERMULATION TEST\n\n")

        if (
            self.permulation_rows_available is not None
            and self.permulation_reps < self.permulation_rows_available
        ):
            if self.sampled_from_available:
                _cprint(
                    f"Using {self.permulation_reps} randomly sampled permulation rows "
                    f"(out of {self.permulation_rows_available} available)."
                )
            else:
                _cprint(
                    f"Using {self.permulation_reps} permulation rows "
                    f"(out of {self.permulation_rows_available} available)."
                )

        if self.max_occ != self.true_odds.total_species_count:
            _cprint(f"** Maximum occupancy set to {self.max_occ} for loss test **\n")
        else:
            pass

        if self.min_occ > 0:
            _cprint(f"** Minimum occupancy set to {self.min_occ} **\n\n")

        if self.permulation_tip_values is not None:
            _cprint(
                f"Using {self.permulation_reps} permulation-derived foreground/background assignments.\n"
            )

        # Initialize lists to store the results of the permulation
        self.means_loss = np.zeros(self.permulation_reps)
        self.means_dup = np.zeros(self.permulation_reps)
        self.stddevs_loss = np.zeros(self.permulation_reps)
        self.stddevs_dup = np.zeros(self.permulation_reps)

        # Initialize counter for how often the means of the permulated
        # distributions exceed the true mean
        counter_mean_loss = 0
        counter_mean_dup = 0

        # Initialize an array which counts how many times the gene's log
        # odds ratio in the permulated distribution is more extreme than the
        # gene's log odds ratio in the true distribution, according to the
        # alternative hypothesis, across all genes in the analysis
        counters_hogs = np.zeros((len(self.true_odds.loss_lor_arr), 4))

        # Get a list of the indices of the species being compared,
        # in case not all species are included in the analysis
        species_incl_arr = (
            self.true_odds.foreground_bool_arr + self.true_odds.background_bool_arr
        ).astype(bool)
        species_incl_idx = species_incl_arr.nonzero()[0]

        for i, perm_tips in enumerate(tqdm(self.permulation_tip_values)):
            counter_mean_loss, counter_mean_dup, counters_hogs = self._permulation_iter(
                i,
                counter_mean_loss,
                counter_mean_dup,
                counters_hogs,
                species_incl_idx,
                perm_tips=perm_tips,
            )

        pval_mean_loss = counter_mean_loss / self.permulation_reps
        pval_mean_dup = counter_mean_dup / self.permulation_reps
        pvals_hogs = counters_hogs / self.permulation_reps
        _cprint(
            f"Permulation counter for MEANS:\n Loss: {str(counter_mean_loss)}, Duplication: {str(counter_mean_dup)}\n"
        )

        self.pval_mean_loss = pval_mean_loss
        self.pval_mean_dup = pval_mean_dup
        self.pvals_hogs = pvals_hogs
        self.results_df = self.true_odds.results_df.copy()
        self.results_df["P-value loss more likely in fg"] = pvals_hogs[:, 0]
        self.results_df["P-value loss more likely in bg"] = pvals_hogs[:, 1]
        self.results_df["P-value duplication more likely in fg"] = pvals_hogs[:, 2]
        self.results_df["P-value duplication more likely in bg"] = pvals_hogs[:, 3]

        # Getting average stats from across the 10000 permulated distributions
        self.loss_mean_av = np.mean(self.means_loss)
        self.loss_stddev_av = np.mean(self.stddevs_loss)
        self.dup_mean_av = np.mean(self.means_dup)
        self.dup_stddev_av = np.mean(self.stddevs_dup)

    def filter_for_permulation_hits(
            self, 
            min_occ=None,
            max_occ=None
            ):
        """Filter the DataFrame based on the occupancy and confidence threshold"""
        df = self.results_df

        if min_occ is None:
            min_occ = self.min_occ
        if max_occ is None:
            max_occ = self.max_occ

        loss_fg_df = df[
            (df["Occupancy"] >= min_occ)
            & (df["Occupancy"] <= max_occ)
            & (df["Log odds ratio of loss"] > self.loss_ci_av[1])
        ].copy()

        loss_fg_df["Significant in permulation test"] = (
            loss_fg_df["P-value loss more likely in fg"] <= self.alpha
        )

        loss_bg_df = df[
            (df["Occupancy"] >= min_occ)
            & (df["Occupancy"] <= max_occ)
            & (df["Log odds ratio of loss"] < self.loss_ci_av[0])
        ].copy()

        loss_bg_df["Significant in permulation test"] = (
            loss_bg_df["P-value loss more likely in bg"] <= self.alpha
        )

        dup_fg_df = df[
            (df["Occupancy"] >= min_occ)
            & (df["Log odds ratio of duplication"] > self.dup_ci_av[1])
        ].copy()

        dup_fg_df["Significant in permulation test"] = (
            dup_fg_df["P-value duplication more likely in fg"] <= self.alpha
        )

        dup_bg_df = df[
            (df["Occupancy"] >= min_occ)
            & (df["Log odds ratio of duplication"] < self.dup_ci_av[0])
        ].copy()

        dup_bg_df["Significant in permulation test"] = (
            dup_bg_df["P-value duplication more likely in bg"] <= self.alpha
        )

        dfs = {
            "loss_fg": loss_fg_df,
            "loss_bg": loss_bg_df,
            "dup_fg": dup_fg_df,
            "dup_bg": dup_bg_df
        }

        df_all = pd.concat(list(dfs.values())).drop_duplicates()
        total_count = len(df_all)

        counts = {
            key: len(df) for key, df in dfs.items()
        }

        return dfs, counts, df_all, total_count
    
    def plot_permulation_stats(
        self,
        test,
        fg_name="foreground",
        bg_name="background",
        include_stddev=True,
        title=True,
        subplot_titles=True,
        hist_color="blue",
        hist_alpha=0.3,
        edgecolor=None,
        legend_fontsize=10,
        axis_label_fontsize=12,
        xlim=None,
        ylim=None,
    ):
        """Plotting the permulation means and standard deviations
        and alpha thresholds to ensure the results are relatively
        tightly distributed"""

        ncols = 2 if include_stddev else 1
        fig_width = 12 if include_stddev else 6.5
        fig, axs = plt.subplots(1, ncols, figsize=(fig_width, 5))
        axs = np.atleast_1d(axs)
        if test=="loss":
            test_name = "loss"
            maximum=self.max_occ
            binwidth=0.05

        elif test=="dup":
            test_name = "duplication"
            maximum=self.true_odds.total_species_count
            binwidth=0.01
        
        else:
            raise ValueError(f"Invalid test type: {test}. Must be 'loss' or 'dup'.")

        if title:
            fig.suptitle(
                f"Permulated (null) distribution stats for gene {test_name},\n"
                f"{fg_name} vs. {bg_name}\n"
                f"Maximum occupancy = {maximum}, minimum occupancy = {self.min_occ}",
                fontsize=16,
            )

        sns.histplot(
            data=getattr(self, f"means_{test}"),
            binwidth=binwidth,
            stat="count",
            ax=axs[0],
            legend=False,
            color=hist_color,
            alpha=hist_alpha,
            edgecolor=edgecolor,
        )

        if subplot_titles:
            axs[0].set_title("permulated means")
        axs[0].set(xlabel="Means", ylabel="Count")
        axs[0].xaxis.label.set_fontsize(axis_label_fontsize)
        axs[0].xaxis.label.set_fontweight("bold")
        axs[0].yaxis.label.set_fontsize(axis_label_fontsize)
        axs[0].yaxis.label.set_fontweight("bold")
        axs[0].axvline(
            x=getattr(self, f"{test}_mean_av"),
            linestyle="dotted",
            color="black",
            label="Avg. permulated mean",
        )
        axs[0].axvline(
            x=getattr(self, f"true_mean_{test}"),
            linestyle="--",
            color="salmon",
            label="True mean",
        )
        axs[0].legend(fontsize=legend_fontsize)
        if xlim is not None:
            axs[0].set_xlim(xlim)
        if ylim is not None:
            axs[0].set_ylim(ylim)

        if include_stddev:
            sns.histplot(
                data=getattr(self, f"stddevs_{test}"),
                binwidth=binwidth,
                stat="count",
                ax=axs[1],
                legend=False,
                color=hist_color,
                alpha=hist_alpha,
                edgecolor=edgecolor,
            )

            if subplot_titles:
                axs[1].set_title("Standard deviations")
            axs[1].set(xlabel="Standard deviations", ylabel="Count")

            axs[1].xaxis.label.set_fontsize(axis_label_fontsize)
            axs[1].xaxis.label.set_fontweight("bold")

            axs[1].yaxis.label.set_fontsize(axis_label_fontsize)
            axs[1].yaxis.label.set_fontweight("bold")

            axs[1].axvline(
                x=getattr(self, f"{test}_stddev_av"),
                linestyle="dotted",
                color="black",
                label="Avg. permulated stddev",
            )
            axs[1].axvline(
                x=getattr(self, f"true_stddev_{test}"),
                linestyle="--",
                color="salmon",
                label="True stddev",
            )
            axs[1].legend(fontsize=legend_fontsize)
            if xlim is not None:
                axs[1].set_xlim(xlim)
            if ylim is not None:
                axs[1].set_ylim(ylim)

        plt.tight_layout()

        return fig, axs

    def plot_permulation_results(
        self,
        test,
        fg_name,
        bg_name="background",
        gaussfit_color="blue",
        avpermulation_color="red",
        hist_color="red",
        thresholds_color="darkred",
        hist_alpha=0.3,
        edgecolor=None,
        bins=100,
        title=True,
        legend_fontsize=10,
        textbox_fontsize=10,
        axis_label_fontsize=12,
    ):
        """Function to plot the results of the permulation test"""

        fig, ax = plt.subplots(figsize=(6, 5))

        if title:
            fig.suptitle(
                f"Log odds ratio of gene {test}, {fg_name} vs. {bg_name}\n"
                f"Maximum occupancy = {self.max_occ}, minimum occupancy = {self.min_occ}",
                fontsize=14,
            )

        # Histogram of the true log odds ratios, filtered for occupancy
        ax.hist(
            getattr(self, f"true_fltrd_{test}_lors"),
            bins=bins,
            density=True,
            color=hist_color,
            alpha=hist_alpha,
            edgecolor=edgecolor,
        )

        x = np.linspace(
            getattr(self, f"true_fltrd_{test}_lors").min(),
            getattr(self, f"true_fltrd_{test}_lors").max(),
            100,
        )

        # Gaussian fit to the true log odds ratios
        ax.plot(
            x,
            norm.pdf(
                x,
                np.mean(getattr(self, f"true_fltrd_{test}_lors")),
                np.std(getattr(self, f"true_fltrd_{test}_lors")),
            ),
            color=gaussfit_color,
            linestyle="--",
            label="Gaussian fit to\ntrue distribution",
        )

        # Normal distribution using the average permulated stats
        ax.plot(
            x,
            norm.pdf(
                x, getattr(self, f"{test}_mean_av"), getattr(self, f"{test}_stddev_av")
            ),
            color=avpermulation_color,
            linestyle="--",
            label="Average permulated\ndistribution",
        )

        # permulation-derived confidence intervals
        ax.axvline(
            x=getattr(self, f"{test}_ci_av")[0],
            label=f"Mean permulated\nthresholds for\nalpha={self.alpha}",
            linestyle="dotted",
            color=thresholds_color,
        )

        ax.axvline(
            x=getattr(self, f"{test}_ci_av")[1],
            linestyle="dotted",
            color=thresholds_color,
        )

        ax.text(
            0.03,
            0.95,
            f"Permulated mean = {self._fmt_stat(getattr(self, f'{test}_mean_av'))}\n"
            f"True mean = {self._fmt_stat(getattr(self, f'true_mean_{test}'))}\n\n"
            f"Permulated std. dev. = {self._fmt_stat(getattr(self, f'{test}_stddev_av'))}\n"
            f"True std. dev. = {self._fmt_stat(getattr(self, f'true_stddev_{test}'))}",
            transform=ax.transAxes,
            fontsize=textbox_fontsize,
            ha="left",
            va="top",
            bbox=dict(
                facecolor="white",
                alpha=0.7,
                edgecolor="0.5",
                linewidth=0.6,
                boxstyle="round,pad=0.2",
            ),
        )

        plt.xlabel("Log odds ratio", fontsize=axis_label_fontsize, fontweight="bold")
        plt.ylabel("Density", fontsize=axis_label_fontsize, fontweight="bold")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        plt.legend(fontsize=legend_fontsize, loc="upper right")
        plt.tight_layout()

        return fig, ax

    def plot_permulation_results_layered(
        self,
        test,
        fg_name,
        bg_name="background",
        gaussfit_color="blue",
        avpermulation_color="red",
        hist_color="red",
        thresholds_color="darkred",
        bins=100,
        title=True,
    ):
        """Function to create 4 sequential plots with elements layered on top of each other.

        Returns 4 figures showing progressive buildup:
        1. Average permulated distribution
        2. + permulation-derived thresholds (with permulated stats)
        3. + Histogram of true log odds ratios (with permulated stats)
        4. + Gaussian fit to the histogram (with true and permulated stats)
        """

        # Define x-axis range
        x = np.linspace(
            getattr(self, f"true_fltrd_{test}_lors").min(),
            getattr(self, f"true_fltrd_{test}_lors").max(),
            100,
        )

        # ----------- Compute a common y-limit so all panels share the same scale -----------
        avg_pdf = norm.pdf(
            x, getattr(self, f"{test}_mean_av"), getattr(self, f"{test}_stddev_av")
        )
        true_pdf = norm.pdf(
            x,
            getattr(self, f"true_mean_{test}"),
            getattr(self, f"true_stddev_{test}"),
        )
        hist_vals, _ = np.histogram(
            getattr(self, f"true_fltrd_{test}_lors"), bins=bins, density=True
        )
        hist_max = hist_vals.max() if hist_vals.size else 0.0
        y_max = max(avg_pdf.max(), true_pdf.max(), hist_max) * 1.05
        if y_max == 0:
            y_max = 1  # fallback to avoid zero-height axis
        # -------------------------------------------------------------------------------

        title_str = (
            f"Log odds ratio of gene {test}, {fg_name} vs. {bg_name}\n"
            f"Maximum occupancy = {self.max_occ}, minimum occupancy = {self.min_occ}"
        )

        figs = []
        axes = []

        # ========== PLOT 1: Average permulated distribution ==========
        fig1, ax1 = plt.subplots(figsize=(8, 6))

        if title:
            fig1.suptitle(title_str, fontsize=14)

        ax1.plot(
            x,
            norm.pdf(
                x, getattr(self, f"{test}_mean_av"), getattr(self, f"{test}_stddev_av")
            ),
            color=avpermulation_color,
            linewidth=2.5,
            label="Average permulated\ndistribution",
        )
        ax1.fill_between(
            x,
            norm.pdf(
                x, getattr(self, f"{test}_mean_av"), getattr(self, f"{test}_stddev_av")
            ),
            alpha=0.2,
            color=avpermulation_color,
        )

        ax1.set_xlabel("Log odds ratio", fontsize=14, fontweight="bold")
        ax1.set_ylabel("Density", fontsize=14, fontweight="bold")
        ax1.set_ylim(bottom=0, top=y_max)
        ax1.set_xlim(x.min(), x.max())
        plt.setp(ax1.get_xticklabels(), fontsize=13)
        plt.setp(ax1.get_yticklabels(), fontsize=13)
        ax1.legend(
            fontsize=13, loc="upper right", ncol=1, labelspacing=0.8, handlelength=1.5
        )
        plt.tight_layout()
        figs.append(fig1)
        axes.append(ax1)

        # ========== PLOT 2: + permulation-derived thresholds ==========
        fig2, ax2 = plt.subplots(figsize=(8, 6))

        if title:
            fig2.suptitle(title_str, fontsize=14)

        ax2.plot(
            x,
            norm.pdf(
                x, getattr(self, f"{test}_mean_av"), getattr(self, f"{test}_stddev_av")
            ),
            color=avpermulation_color,
            linewidth=2.5,
            label="Average permulated\ndistribution",
        )
        ax2.fill_between(
            x,
            norm.pdf(
                x, getattr(self, f"{test}_mean_av"), getattr(self, f"{test}_stddev_av")
            ),
            alpha=0.2,
            color=avpermulation_color,
            zorder=0,
        )

        ax2.axvline(
            x=getattr(self, f"{test}_ci_av")[0],
            label=f"Mean permulated\nthresholds for\nalpha={self.alpha}",
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )
        ax2.axvline(
            x=getattr(self, f"{test}_ci_av")[1],
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )

        ax2.text(
            0.03,
            0.95,
            f"permulated mean = {self._fmt_stat(getattr(self, f'{test}_mean_av'))}\n"
            f"permulated std. dev. = {self._fmt_stat(getattr(self, f'{test}_stddev_av'))}",
            transform=ax2.transAxes,
            fontsize=12,
            ha="left",
            va="top",
            bbox=dict(
                facecolor="white",
                alpha=0.7,
                edgecolor="0.5",
                linewidth=0.6,
                boxstyle="round,pad=0.2",
            ),
        )

        ax2.set_xlabel("Log odds ratio", fontsize=14, fontweight="bold")
        ax2.set_ylabel("Density", fontsize=14, fontweight="bold")
        ax2.set_ylim(bottom=0, top=y_max)
        ax2.set_xlim(x.min(), x.max())
        plt.setp(ax2.get_xticklabels(), fontsize=13)
        plt.setp(ax2.get_yticklabels(), fontsize=13)
        ax2.legend(
            fontsize=13, loc="upper right", ncol=1, labelspacing=0.8, handlelength=1.5
        )
        plt.tight_layout()
        figs.append(fig2)
        axes.append(ax2)

        # ========== PLOT 3: + Histogram ==========
        fig3, ax3 = plt.subplots(figsize=(8, 6))

        if title:
            fig3.suptitle(title_str, fontsize=14)

        ax3.plot(
            x,
            norm.pdf(
                x, getattr(self, f"{test}_mean_av"), getattr(self, f"{test}_stddev_av")
            ),
            color=avpermulation_color,
            linewidth=2.5,
            label="Average permulated\ndistribution",
        )
        ax3.fill_between(
            x,
            norm.pdf(
                x, getattr(self, f"{test}_mean_av"), getattr(self, f"{test}_stddev_av")
            ),
            alpha=0.2,
            color=avpermulation_color,
            zorder=0,
        )

        ax3.axvline(
            x=getattr(self, f"{test}_ci_av")[0],
            label=f"Mean permulated\nthresholds for\nalpha={self.alpha}",
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )
        ax3.axvline(
            x=getattr(self, f"{test}_ci_av")[1],
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )

        ax3.hist(
            getattr(self, f"true_fltrd_{test}_lors"),
            bins=bins,
            density=True,
            color=hist_color,
            alpha=0.3,
            edgecolor=hist_color,
            label="True distribution",
            zorder=3,
        )

        ax3.text(
            0.03,
            0.95,
            f"permulated mean = {self._fmt_stat(getattr(self, f'{test}_mean_av'))}\n"
            f"permulated std. dev. = {self._fmt_stat(getattr(self, f'{test}_stddev_av'))}",
            transform=ax3.transAxes,
            fontsize=12,
            ha="left",
            va="top",
            bbox=dict(
                facecolor="white",
                alpha=0.7,
                edgecolor="0.5",
                linewidth=0.6,
                boxstyle="round,pad=0.2",
            ),
        )

        ax3.set_xlabel("Log odds ratio", fontsize=14, fontweight="bold")
        ax3.set_ylabel("Density", fontsize=14, fontweight="bold")
        ax3.set_xlim(x.min(), x.max())
        ax3.set_ylim(bottom=0, top=y_max)
        plt.setp(ax3.get_xticklabels(), fontsize=13)
        plt.setp(ax3.get_yticklabels(), fontsize=13)
        ax3.legend(
            fontsize=13, loc="upper right", ncol=1, labelspacing=0.8, handlelength=1.5
        )
        plt.tight_layout()
        figs.append(fig3)
        axes.append(ax3)

        # ========== PLOT 4: + Gaussian fit ==========
        fig4, ax4 = plt.subplots(figsize=(8, 6))

        if title:
            fig4.suptitle(title_str, fontsize=14)

        ax4.plot(
            x,
            norm.pdf(
                x, getattr(self, f"{test}_mean_av"), getattr(self, f"{test}_stddev_av")
            ),
            color=avpermulation_color,
            linewidth=2.5,
            label="Average permulated\ndistribution",
        )
        ax4.fill_between(
            x,
            norm.pdf(
                x, getattr(self, f"{test}_mean_av"), getattr(self, f"{test}_stddev_av")
            ),
            alpha=0.2,
            color=avpermulation_color,
            zorder=0,
        )

        ax4.axvline(
            x=getattr(self, f"{test}_ci_av")[0],
            label=f"Mean permulated\nthresholds for\nalpha={self.alpha}",
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )
        ax4.axvline(
            x=getattr(self, f"{test}_ci_av")[1],
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )

        ax4.hist(
            getattr(self, f"true_fltrd_{test}_lors"),
            bins=bins,
            density=True,
            color=hist_color,
            alpha=0.3,
            edgecolor=hist_color,
            label="True distribution",
            zorder=3,
        )

        ax4.plot(
            x,
            norm.pdf(
                x,
                getattr(self, f"true_mean_{test}"),
                getattr(self, f"true_stddev_{test}"),
            ),
            color=gaussfit_color,
            linestyle="--",
            linewidth=2,
            label="Gaussian fit to\ntrue distribution",
        )

        ax4.text(
            0.03,
            0.95,
            f"Permulated mean = {self._fmt_stat(getattr(self, f'{test}_mean_av'))}\n"
            f"Permulated std. dev. = {self._fmt_stat(getattr(self, f'{test}_stddev_av'))}\n\n"
            f"True mean = {self._fmt_stat(getattr(self, f'{test}_true_mean'))}\n"
            f"True std. dev. = {self._fmt_stat(getattr(self, f'{test}_true_stddev'))}",
            transform=ax4.transAxes,
            fontsize=12,
            ha="left",
            va="top",
            bbox=dict(
                facecolor="white",
                alpha=0.7,
                edgecolor="0.5",
                linewidth=0.6,
                boxstyle="round,pad=0.2",
            ),
        )

        ax4.set_xlabel("Log odds ratio", fontsize=14, fontweight="bold")
        ax4.set_ylabel("Density", fontsize=14, fontweight="bold")
        ax4.set_xlim(x.min(), x.max())
        ax4.set_ylim(bottom=0, top=y_max)
        plt.setp(ax4.get_xticklabels(), fontsize=13)
        plt.setp(ax4.get_yticklabels(), fontsize=13)
        ax4.legend(
            fontsize=13, loc="upper right", ncol=1, labelspacing=0.8, handlelength=1.5
        )
        plt.tight_layout()
        figs.append(fig4)
        axes.append(ax4)

        return figs, axes

    def _get_hits_dfs(self):
        """Filter the dataframe for HOGs that have LORs exceeding the
        alpha values determined by permulation, and filter for those
        which include a species of interest if specified."""

        self.get_pval_thresholds(self.alpha, "two-tailed")

        dfs, counts, df_all, total_count = self.filter_for_permulation_hits()
        self.results_fltrd_dfs = dfs
        self.results_fltrd_df_all = df_all
        self.all_hits_count = total_count
        self.counts_hits = counts

    def print_permulation_results(self, fname=sys.stdout):
        """Function to print the results of the permulation test"""

        def _write_results(file_obj):

            _emit(
                "*********************** RESULTS ***********************\n\n"
                f"Permulation test with {self.permulation_reps} repetitions\n"
                f"Minimum occupancy: {self.min_occ}; max occupancy (loss): {self.max_occ}\n"
                f"Analysis run on {self.true_odds.time}",
                file_obj,
            )

            _emit(
                f"Foreground list: {self.true_odds.foreground_list_filename}", file_obj
            )
            if self.true_odds.background_list_filename is not None:
                _emit(
                    f"Background list: {self.true_odds.background_list_filename}",
                    file_obj,
                )
            _emit(f"Gene count file: {self.true_odds.genecount_csv}", file_obj)
            _emit(
                f"Hierarchical orthogroup file: {self.true_odds.hog_node_genes_tsv}\n",
                file_obj,
            )

            _emit(
                f"Total number of HOGs in node: {len(self.true_odds.hog_list)}\n"
                f"Total number of HOGs within occupancy threshold: "
                f"{len(self.true_fltrd_loss_lors)} (loss), "
                f"{len(self.true_fltrd_dup_lors)} (duplication)\n"
                f"Total species: {self.true_odds.total_species_count}\n"
                f"Foreground count: {self.true_odds.foreground_count}\n"
                f"Background count: {self.true_odds.background_count}\n"
                f"True mean, loss: {self._fmt_stat(self.true_mean_loss, ndigits=3)}\n"
                f"True standard deviation, loss: {self._fmt_stat(self.true_stddev_loss, ndigits=3)}\n"
                f"True mean, duplication: {self._fmt_stat(self.true_mean_dup, ndigits=3)}\n"
                f"True standard deviation, duplication: {self._fmt_stat(self.true_stddev_dup, ndigits=3)}\n\n",
                file_obj,
            )

            _emit(
                "** Permulation P-VALUES ** \n\n"
                f"Probability that the null is true for MEAN, loss (alt=less): {self._fmt_stat(self.pval_mean_loss, ndigits=3)}\n"
                f"Probability that the null is true for MEAN, duplication (alt=greater): {self._fmt_stat(self.pval_mean_dup, ndigits=3)}\n\n",
                file_obj,
            )

            _emit(
                f"Permulated average mean, loss: {self._fmt_stat(self.loss_mean_av, ndigits=3)}\n"
                f"Permulated average standard deviation, loss: {self._fmt_stat(self.loss_stddev_av, ndigits=3)}\n"
                f"Permulated average mean, duplication: {self._fmt_stat(self.dup_mean_av, ndigits=3)}\n"
                f"Permulated average standard deviation, duplication: {self._fmt_stat(self.dup_stddev_av, ndigits=3)}\n"
                f"User-defined significance threshold: {self.alpha}\n"
                f"Permulation-derived alpha threshold, loss: {self._fmt_ci(self.loss_ci_av, ndigits=3)}\n"
                f"Permulation-derived alpha threshold, duplication: {self._fmt_ci(self.dup_ci_av, ndigits=3)}\n\n"
                f"Significant HOGs (all): {self.all_hits_count}\n",
                file_obj,
            )
            _emit(
                    f"  - loss_fg: {self.counts_hits['loss_fg']}\n"
                    f"  - loss_bg: {self.counts_hits['loss_bg']}\n"
                    f"  - dup_fg: {self.counts_hits['dup_fg']}\n"
                    f"  - dup_bg: {self.counts_hits['dup_bg']}\n",
                    file_obj,
                )

            if self.species_of_interest is not None:
                sp_of_int_counts = {
                    "loss_fg": filter_for_sp_of_interest(self.results_fltrd_dfs['loss_fg'], self.true_odds.genecount_df, self.species_of_interest),
                    "loss_bg": filter_for_sp_of_interest(self.results_fltrd_dfs['loss_bg'], self.true_odds.genecount_df, self.species_of_interest),
                    "dup_fg": filter_for_sp_of_interest(self.results_fltrd_dfs['dup_fg'], self.true_odds.genecount_df, self.species_of_interest),
                    "dup_bg": filter_for_sp_of_interest(self.results_fltrd_dfs['dup_bg'], self.true_odds.genecount_df, self.species_of_interest)
                }
                _emit(
                    f"Significant HOGs, {self.species_of_interest} present:\n",
                    file_obj,
                )
                _emit(
                    f"  - loss_fg: {sp_of_int_counts['loss_fg']}\n"
                    f"  - loss_bg: {sp_of_int_counts['loss_bg']}\n"
                    f"  - dup_fg: {sp_of_int_counts['dup_fg']}\n"
                    f"  - dup_bg: {sp_of_int_counts['dup_bg']}\n",
                    file_obj,
                )

                _emit(
                    f"Significant HOGs, using permulation p-values:\n",
                    file_obj,
                )
                _emit(
                    f"  - loss_fg: {self.results_fltrd_dfs['loss_fg'][self.results_fltrd_dfs['loss_fg']['Significant in permulation test'] == True].shape[0]}\n"
                    f"  - loss_bg: {self.results_fltrd_dfs['loss_bg'][self.results_fltrd_dfs['loss_bg']['Significant in permulation test'] == True].shape[0]}\n"
                    f"  - dup_fg: {self.results_fltrd_dfs['dup_fg'][self.results_fltrd_dfs['dup_fg']['Significant in permulation test'] == True].shape[0]}\n"
                    f"  - dup_bg: {self.results_fltrd_dfs['dup_bg'][self.results_fltrd_dfs['dup_bg']['Significant in permulation test'] == True].shape[0]}\n",
                    file_obj,
                )

        if fname is sys.stdout:
            _write_results(sys.stdout)
        else:
            with open(fname, "w", encoding="utf-8") as file_obj:
                _write_results(file_obj)

    def save_pickle_file(self, fname):
        """
        Saves the permulation results object to a pickle file.
        """

        with open(fname, "wb") as file:
            pickle.dump(self, file)

    @classmethod
    def load_from_pickle(cls, filepath: str) -> "PermulationTestResults":
        """Load results object from a pickle file."""
        with open(filepath, "rb") as f:
            return pickle.load(f)

    def save_results_files(
        self, results_dir, save_pickle, fg_name, bg_name="background"
    ):
        """
        Takes in permulation test results instance and saves all
        relevant plots and tables.
        """

        # Save the permulation results object as a pickle file
        if save_pickle:
            self.save_pickle_file(fname=f"{results_dir}/results.pkl")
        else:
            _cprint(
                "save_pickle is set to False, so the full results object will not be saved as a pickle file.\n"
                "To save the results object for future use, set save_pickle to True. If you stored the test \n"
                "results object in a variable, you can also save it as a pickle file using the `save_pickle_file` method, e.g. \n"
                "my_results.save_pickle_file('path/to/save/my_results.pkl')\n"
            )

        # Save the full table of odds ratios and p-values to a csv
        # Overwrites the results_all.csv from OddsRatioTest
        self.results_df.to_csv(f"{results_dir}/results_all.csv", index=True)

        self.results_fltrd_df_all.to_csv((f"{results_dir}/" + "fltrd_hits.csv"), index=True)

        # Save a text file summarizing results from the analysis
        self.print_permulation_results(fname=f"{results_dir}/" + "results_summary.txt")

        #### Save figures ####

        # permulation statistics
        lossfig, _ = self.plot_permulation_stats("loss", fg_name, bg_name)

        lossfig.savefig(
            os.path.join(results_dir, "loss_stats_dists.png"),
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.3,
        )

        dupfig, _ = self.plot_permulation_stats(
            "dup",
            fg_name,
            bg_name,
        )

        dupfig.savefig(
            os.path.join(results_dir, "dup_stats_dists.png"),
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.3,
        )

        # True distribution overlaid with average permulated distribution
        lossfig, _ = self.plot_permulation_results("loss", fg_name, bg_name)

        lossfig.savefig(
            os.path.join(results_dir, "loss_results.png"),
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.3,
        )

        dupfig, _ = self.plot_permulation_results("dup", fg_name, bg_name)
        dupfig.savefig(
            os.path.join(results_dir, "dup_results.png"),
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.3,
        )
        rel_results_dir = os.path.relpath(results_dir, _REPO_ROOT)
        
        _cprint(
            "Results files saved.\n"
            f"Directory: {rel_results_dir}\n\n"
            "Files include:\n"
            "  - results_all.csv: all HOG odds/log-odds and p-values\n"
            "  - fltrd_hits.csv: significant hits after occupancy/p-value filters\n"
            "  - results_summary.txt: text summary of run settings and statistics\n"
            "  - loss_stats_dists.png: permulated mean/stddev histograms for loss\n"
            "  - dup_stats_dists.png: permulated mean/stddev histograms for duplication\n"
            "  - loss_results.png: loss LOR true distribution vs permulated average\n"
            "  - dup_results.png: duplication LOR true distribution vs permulated average\n"
            "  - results.pkl: full results object (only when save_pickle=True)"
        )


def odds_ratio_test(
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
):
    """Run the full odds ratio test."""

    if results_dir is not None:
        if fg_name is None:
            raise ValueError(
                "Please provide a descriptive name for your test group, \n"
                "e.g. fg_name = 'orbweavers'. "
                "You may also name your background, e.g.\n"
                "bg_name = 'non-orbweavers'"
            )

    time = datetime.now()
    time_fmtd = time.strftime("%Y-%m-%d at %H:%M:%S")

    foreground_list_filename = _resolve_repo_path(foreground_list_filename)
    background_list_filename = _resolve_repo_path(background_list_filename)
    hog_node_genes_tsv = _resolve_repo_path(hog_node_genes_tsv)
    genecount_csv = _resolve_repo_path(genecount_csv)
    buscos_filename = (
        _resolve_repo_path(buscos_filename) if correct_for_buscos else None
    )
    permulations_tip_values_csv = _resolve_repo_path(permulations_tip_values_csv)
    results_dir = _resolve_repo_path(results_dir)

    # Generate the genecount file if not provided
    if genecount_csv is None:
        genecount_csv = orthogroup_gene_count.main(hog_node_genes_tsv)

    # Get the true odds ratio results
    true_odds = OddsRatioResults(
        genecount_csv=genecount_csv,
        hog_node_genes_tsv=hog_node_genes_tsv,
        time=time_fmtd,
        foreground_list_filename=foreground_list_filename,
        background_list_filename=background_list_filename,
        buscos_filename=buscos_filename,
    )

    # Determine permulation assignments source.
    # Priority: explicit in-memory values > provided CSV > random permulation.
    permulation_tip_values = load_permulation_tip_values_from_csv(
        permulations_tip_values_csv
    )

    # Cap to the requested number of reps when the CSV contains more rows than needed.
    if permulation_reps <= 0:
        raise ValueError("permulation_reps must be a positive integer.")

    permulation_rows_available = len(permulation_tip_values)
    effective_permulation_reps = min(permulation_reps, permulation_rows_available)
    sampled_from_available = False
    if effective_permulation_reps < permulation_rows_available:
        permulation_tip_values = random.sample(
            permulation_tip_values, k=effective_permulation_reps
        )
        sampled_from_available = True

    if results_dir is not None:
        unique_results_dir = _unique_results_dir(
            results_dir,
            time,
            min_occ,
            max_occ,
            permulation_reps,
        )
        os.makedirs(unique_results_dir, exist_ok=True)
        true_odds.results_df.to_csv(f"{unique_results_dir}/results_all.csv", index=True)

    permulation_test_results = PermulationTestResults(
        true_odds=true_odds,
        permulation_tip_values=permulation_tip_values,
        min_occ=min_occ,
        max_occ=max_occ,
        alpha=alpha,
        permulation_reps=effective_permulation_reps,
        permulation_rows_available=permulation_rows_available,
        sampled_from_available=sampled_from_available,
        species_of_interest=species_of_interest,
    )

    if results_dir is not None:
        permulation_test_results.save_results_files(
            results_dir=unique_results_dir,
            save_pickle=save_pickle,
            fg_name=fg_name,
            bg_name=bg_name,
        )

    else:
        _cprint(
            "No results directory provided, so plots and tables will not be saved to files.\n"
            "To save results files, provide a path to a results directory using the results_dir argument."
        )
        if save_pickle:
            _cprint(
                "save_pickle is set to True, but no results directory provided, so the full results object will not be saved as a pickle file.\n"
                "To save the results object for future use, set save_pickle to True and provide a path to a results directory using the results_dir argument."
            )

        # Show plots interactively even when not saving results to disk.
        display_fg_name = fg_name if fg_name is not None else "foreground"
        display_bg_name = bg_name if bg_name is not None else "background"

        permulation_test_results.plot_permulation_stats(
            display_fg_name, display_bg_name
        )

        permulation_test_results.plot_permulation_results(
            display_fg_name, display_bg_name
        )
        plt.show()

    return permulation_test_results


# Backward-compatible alias for old pickle files
PermutationTestResults = PermulationTestResults
