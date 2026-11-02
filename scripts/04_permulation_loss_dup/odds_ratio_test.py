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
from typing import Dict, Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from tqdm.auto import tqdm
from odds_ratio_test_helpers import (
    _cprint,
    _emit,
    _resolve_repo_path,
    _unique_results_dir,
    benjamini_hochberg,
    calculate_odds,
    calculate_odds_batch,
    calculate_odds_batch_weighted,
    directional_q_values,
    load_flagged_cells,
    drop_empty_cols,
    filter_for_sp_of_interest,
    load_permulation_tip_values_from_csv,
    occupancy_filter,
    permulated_foreground_matrices,
)
import odds_ratio_test_plotting as plotting
import orthogroup_gene_count


# Set the random seed for reproducibility
random.seed(42)

_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_SRC_DIR, os.pardir, os.pardir))
CONSOLE_PRINT_WIDTH = 96

# Permulated log odds ratios within this distance of the true value count as ties. Batch and
# per-species arithmetic differ by ~1e-14, so exact float comparison would split genuine ties arbitrarily.
TIE_TOLERANCE = 1e-9


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
        flagged_cells=None,
    ):
        """Initialize the inputs for the odds ratio calculations

        flagged_cells: optional DataFrame (HOG, species, n_flagged) of gene copies to treat as
        contaminated. They are removed from the gene counts, and a species with no copy left in a
        HOG is treated as missing (not absent) in that HOG only.
        """
        self.flagged_cells = flagged_cells
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

        # HOGs with a missing (flagged) cell weight that species 0, so they need the weighted calculation
        if len(self.flagged_hog_idx):
            rows = self.flagged_hog_idx
            for test in ("loss", "dup"):
                lor_arr = getattr(self, f"{test}_lor_arr")
                lor_arr[rows] = calculate_odds_batch_weighted(
                    self.foreground_bool_arr,
                    self.background_bool_arr,
                    getattr(self, f"{test}_bool_mat")[rows],
                    self.weight_matrix(test, rows),
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

        # Flagged gene copies (suspected contamination) are removed from the counts; a species with no
        # copy left in a HOG is missing there (not absent) and gets weight 0 in that HOG's tests
        missing_mat = np.zeros(genecount_df.shape, dtype=bool)
        self.occ_unflagged_arr = genecount_df.to_numpy().astype(bool).sum(axis=1)
        genecount_df, missing_mat = self._apply_flagged_cells(genecount_df, missing_mat)

        # Convert counts df to numpy
        genecount_mat = genecount_df.to_numpy()

        # Matrix where 1=species is present in the HOG, 0=species has no gene in the HOG
        occ_bool_mat = genecount_mat.astype(bool).astype(float)

        # Calculate the total occupancy of each HOG (number of species present)
        total_occ_arr = occ_bool_mat.sum(axis=1).reshape(genecount_mat.shape[0], 1)

        # Matrix where 1=species is absent from the HOG, 0=species is present (or missing)
        loss_bool_mat = np.asarray((genecount_mat == 0) & ~missing_mat).astype(float)

        # Matrix where 1=gene with more than 1 copy, 0=gene with 1 or no copies
        dup_bool_mat = np.asarray(genecount_mat > 1).astype(float)

        self.genecount_df = genecount_df
        self.missing_mat = missing_mat
        self.flagged_hog_idx = np.where(missing_mat.any(axis=1))[0]
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

    def _apply_flagged_cells(self, genecount_df, missing_mat):
        """Subtract flagged copies from the counts and mark cells left with no copies as missing."""
        flagged = self.flagged_cells
        if flagged is None or len(flagged) == 0:
            return genecount_df, missing_mat

        hog_pos = {hog: i for i, hog in enumerate(genecount_df.index)}
        sp_pos = {sp: j for j, sp in enumerate(genecount_df.columns)}
        counts = genecount_df.to_numpy().copy()
        n_applied = n_skipped = n_clipped = 0
        for hog, species, n_flagged in flagged[["HOG", "species", "n_flagged"]].itertuples(index=False):
            if hog not in hog_pos or species not in sp_pos:
                n_skipped += 1
                continue
            i, j = hog_pos[hog], sp_pos[species]
            if counts[i, j] == 0:
                n_skipped += 1  # nothing to remove: the species has no copy in this HOG
                continue
            if n_flagged > counts[i, j]:
                n_clipped += 1
            counts[i, j] = max(counts[i, j] - n_flagged, 0)
            missing_mat[i, j] = counts[i, j] == 0
            n_applied += 1
        _cprint(
            f"Flagged cells: {n_applied} applied ({int(missing_mat.sum())} left missing), "
            f"{n_skipped} skipped (not in the analysis or no copy), {n_clipped} had more flags than copies"
        )
        return pd.DataFrame(counts, index=genecount_df.index, columns=genecount_df.columns), missing_mat

    def weight_matrix(self, test, rows):
        """(len(rows), n_species) species weights for a test: the BUSCO weight, 0 for missing cells."""
        busco = getattr(self, f"{test}_busco_arr", None)
        base = np.ones(len(self.all_species_arr)) if busco is None else np.asarray(busco, dtype=float)
        return base[None, :] * (~self.missing_mat[rows])

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
        if self.flagged_cells is not None and len(self.flagged_cells):
            results_df["Missing flagged cells"] = self.missing_mat.sum(axis=1)
            results_df["Occupancy (before flagging)"] = self.occ_unflagged_arr
        results_df = results_df.set_index("HOG")

        return results_df


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

    def _count_species_of_interest_hits(self, species_name):
        return {
            key: filter_for_sp_of_interest(
                self.results_fltrd_dfs[key],
                self.true_odds.genecount_df,
                species_name,
            )
            for key in ("loss_fg", "loss_bg", "dup_fg", "dup_bg")
        }

    def _permulation_significant_mask(self, df: pd.DataFrame) -> pd.Series:
        """Return a robust boolean mask for permulation-significant rows."""
        col = "Significant by avgd thresholds"
        if col not in df.columns:
            return pd.Series(False, index=df.index)

        vals = df[col]
        if pd.api.types.is_bool_dtype(vals):
            return vals.fillna(False)

        # Column stores labels like "loss_fg"/"dup_bg" or empty strings.
        return vals.fillna("").astype(str).ne("")

    def _report_flag_effect(self, loss_idx, dup_idx):
        """Print how flagging changed which HOGs are tested."""
        occ_before = self.true_odds.occ_unflagged_arr
        occ_after = self.true_odds.total_occ_arr.ravel()
        tested_loss_before = (occ_before >= self.min_occ) & (occ_before <= self.max_occ)
        tested_dup_before = occ_before >= self.min_occ
        n_hogs = len(self.true_odds.flagged_hog_idx)
        _cprint(
            f"Flagging touched {n_hogs} HOGs; occupancy changed in {int((occ_before != occ_after).sum())}. "
            f"Tested for loss: {len(loss_idx)} (was {int(tested_loss_before.sum())}), "
            f"for duplication: {len(dup_idx)} (was {int(tested_dup_before.sum())})."
        )

    def _count_permulation_significant_hits(self):
        return {
            key: self.results_fltrd_dfs[key][
                self._permulation_significant_mask(self.results_fltrd_dfs[key])
            ].shape[0]
            for key in ("loss_fg", "loss_bg", "dup_fg", "dup_bg")
        }

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

        # Foreground/background labelings for all permulations as species-by-permulation matrices
        foreground_mat, background_mat = permulated_foreground_matrices(
            self.true_odds.all_species_arr, self.permulation_tip_values, species_incl_idx
        )
        if not foreground_mat.any():
            raise ValueError(
                "No species names from permulation tips matched analysis species names."
            )

        # Only the HOGs that are actually tested are evaluated: loss within min_occ-max_occ,
        # duplication at min_occ and above. HOGs outside those ranges get no p-value (NaN).
        occ = self.true_odds.total_occ_arr.ravel()
        loss_idx = np.where((occ >= self.min_occ) & (occ <= self.max_occ))[0]
        dup_idx = np.where(occ >= self.min_occ)[0]

        loss_mat = self.true_odds.loss_bool_mat[loss_idx]
        dup_mat = self.true_odds.dup_bool_mat[dup_idx]
        loss_busco = getattr(self.true_odds, "loss_busco_arr", None)
        dup_busco = getattr(self.true_odds, "dup_busco_arr", None)
        true_loss = self.true_odds.loss_lor_arr[loss_idx]
        true_dup = self.true_odds.dup_lor_arr[dup_idx]

        # Tested HOGs with a missing (flagged) cell: positions within loss_idx / dup_idx, and their weights
        flagged_hogs = self.true_odds.flagged_hog_idx
        loss_flag_pos = np.where(np.isin(loss_idx, flagged_hogs))[0]
        dup_flag_pos = np.where(np.isin(dup_idx, flagged_hogs))[0]
        loss_flag_w = self.true_odds.weight_matrix("loss", loss_idx[loss_flag_pos])
        dup_flag_w = self.true_odds.weight_matrix("dup", dup_idx[dup_flag_pos])
        if len(flagged_hogs):
            self._report_flag_effect(loss_idx, dup_idx)

        # Permulations are processed in chunks to bound memory (n_hogs x chunk matrices)
        chunk_size = 100
        for start in tqdm(range(0, self.permulation_reps, chunk_size)):
            stop = min(start + chunk_size, self.permulation_reps)
            fg_chunk = foreground_mat[:, start:stop]
            bg_chunk = background_mat[:, start:stop]

            new_loss = calculate_odds_batch(fg_chunk, bg_chunk, loss_mat, busco_arr=loss_busco)
            new_dup = calculate_odds_batch(fg_chunk, bg_chunk, dup_mat, busco_arr=dup_busco)
            if len(loss_flag_pos):
                new_loss[loss_flag_pos] = calculate_odds_batch_weighted(
                    fg_chunk, bg_chunk, loss_mat[loss_flag_pos], loss_flag_w
                )
            if len(dup_flag_pos):
                new_dup[dup_flag_pos] = calculate_odds_batch_weighted(
                    fg_chunk, bg_chunk, dup_mat[dup_flag_pos], dup_flag_w
                )

            # Count, per HOG, the permulations whose log odds ratio is at least as extreme as the
            # true one in each direction. Ties count as "as extreme" (inclusive), the standard
            # permutation-test convention.
            counters_hogs[loss_idx, 0] += (new_loss >= true_loss - TIE_TOLERANCE).sum(axis=1)
            counters_hogs[loss_idx, 1] += (new_loss <= true_loss + TIE_TOLERANCE).sum(axis=1)
            counters_hogs[dup_idx, 2] += (new_dup >= true_dup - TIE_TOLERANCE).sum(axis=1)
            counters_hogs[dup_idx, 3] += (new_dup <= true_dup + TIE_TOLERANCE).sum(axis=1)

            # Mean and standard deviation of each permulated distribution over the tested HOGs
            self.means_loss[start:stop] = new_loss.mean(axis=0)
            self.stddevs_loss[start:stop] = new_loss.std(axis=0)
            self.means_dup[start:stop] = new_dup.mean(axis=0)
            self.stddevs_dup[start:stop] = new_dup.std(axis=0)

        # Compare each permulated mean with the true mean in the direction of the alternative
        counter_mean_loss = int((self.means_loss <= self.true_mean_loss).sum())
        counter_mean_dup = int((self.means_dup >= self.true_mean_dup).sum())

        # Add-one estimator: the observed labeling counts as one more draw from the null, so a
        # p-value is never 0 and cannot go below 1 / (permulation_reps + 1)
        pval_mean_loss = (counter_mean_loss + 1) / (self.permulation_reps + 1)
        pval_mean_dup = (counter_mean_dup + 1) / (self.permulation_reps + 1)
        pvals_hogs = (counters_hogs + 1) / (self.permulation_reps + 1)
        # HOGs outside a test's occupancy range were not tested for it
        pvals_hogs[np.setdiff1d(np.arange(pvals_hogs.shape[0]), loss_idx), 0:2] = np.nan
        pvals_hogs[np.setdiff1d(np.arange(pvals_hogs.shape[0]), dup_idx), 2:4] = np.nan
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

        # Benjamini-Hochberg q-values, computed only over the HOGs actually tested (loss:
        # min_occ-max_occ; duplication: min_occ and above; NaN elsewhere). The two directions of a
        # test are two tails of one comparison, so each test's fg and bg p-values are stacked and
        # corrected together.
        occ = self.results_df["Occupancy"]
        tested = {
            "loss": (occ >= self.min_occ) & (occ <= self.max_occ),
            "duplication": occ >= self.min_occ,
        }
        for test_name, in_universe in tested.items():
            fg_col = f"P-value {test_name} more likely in fg"
            bg_col = f"P-value {test_name} more likely in bg"
            stacked = np.concatenate(
                [
                    self.results_df.loc[in_universe, fg_col].to_numpy(),
                    self.results_df.loc[in_universe, bg_col].to_numpy(),
                ]
            )
            qvals = benjamini_hochberg(stacked)
            n_tested = int(in_universe.sum())
            for group, block in (("fg", qvals[:n_tested]), ("bg", qvals[n_tested:])):
                qcol = f"Q-value {test_name} more likely in {group}"
                self.results_df[qcol] = np.nan
                self.results_df.loc[in_universe, qcol] = block

        # Pre-specified directions (loss in the background, duplication in the foreground), each
        # corrected within its own tail
        self.results_df = self.results_df.join(
            directional_q_values(self.results_df, self.min_occ, self.max_occ)
        )

        # Getting average stats from across the 10000 permulated distributions
        self.loss_mean_av = np.mean(self.means_loss)
        self.loss_stddev_av = np.mean(self.stddevs_loss)
        self.dup_mean_av = np.mean(self.means_dup)
        self.dup_stddev_av = np.mean(self.stddevs_dup)

    # Maybe switch which one is filtered first?
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

        # alpha is the two-sided level; each direction is tested at alpha / 2
        tail_alpha = self.alpha / 2

        loss_fg_df = df[
            (df["Occupancy"] >= min_occ)
            & (df["Occupancy"] <= max_occ)
            & (df["P-value loss more likely in fg"] <= tail_alpha)
        ].copy()

        loss_fg_df["Significant by avgd thresholds"] = np.where(
            loss_fg_df["Log odds ratio of loss"] > self.loss_ci_av[1], "loss_fg", None
        )

        loss_bg_df = df[
            (df["Occupancy"] >= min_occ)
            & (df["Occupancy"] <= max_occ)
            & (df["P-value loss more likely in bg"] <= tail_alpha)
        ].copy()

        loss_bg_df["Significant by avgd thresholds"] = np.where(
            loss_bg_df["Log odds ratio of loss"] < self.loss_ci_av[0], "loss_bg", None
        )

        dup_fg_df = df[
            (df["Occupancy"] >= min_occ)
             & (df["P-value duplication more likely in fg"] <= tail_alpha)
            
        ].copy()

        dup_fg_df["Significant by avgd thresholds"] = np.where(
            dup_fg_df["Log odds ratio of duplication"] > self.dup_ci_av[1],
            "dup_fg",
            None,
        )

        dup_bg_df = df[
            (df["Occupancy"] >= min_occ)
            & (df["P-value duplication more likely in bg"] <= tail_alpha)
        ].copy()

        dup_bg_df["Significant by avgd thresholds"] = np.where(
            dup_bg_df["Log odds ratio of duplication"] < self.dup_ci_av[0],
            "dup_bg",
            None,
        )

        dfs = {
            "loss_fg": loss_fg_df,
            "loss_bg": loss_bg_df,
            "dup_fg": dup_fg_df,
            "dup_bg": dup_bg_df
        }

        loss_fg_df["Significant by permulation"] = "loss_fg"
        loss_bg_df["Significant by permulation"] = "loss_bg"
        dup_fg_df["Significant by permulation"] = "dup_fg"
        dup_bg_df["Significant by permulation"] = "dup_bg"

        df = pd.concat(list(dfs.values())).drop_duplicates()

        df_all = (
            df.reset_index()
            .groupby("HOG", as_index=False)
            .agg(lambda s: ", ".join(pd.unique(s.dropna().astype(str))))
            .set_index("HOG")
            )

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
        subplot_titles=False,
        hist_color="blue",
        hist_alpha=0.3,
        edgecolor=None,
        legend_fontsize=10,
        axis_label_fontsize=12,
        xlim=None,
        ylim=None,
        binwidth=None,
    ):
        return plotting.plot_permulation_stats(
            self,
            test,
            fg_name=fg_name,
            bg_name=bg_name,
            include_stddev=include_stddev,
            title=title,
            subplot_titles=subplot_titles,
            hist_color=hist_color,
            hist_alpha=hist_alpha,
            edgecolor=edgecolor,
            legend_fontsize=legend_fontsize,
            axis_label_fontsize=axis_label_fontsize,
            xlim=xlim,
            ylim=ylim,
            binwidth=binwidth,
        )

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
        return plotting.plot_permulation_results(
            self,
            test,
            fg_name=fg_name,
            bg_name=bg_name,
            gaussfit_color=gaussfit_color,
            avpermulation_color=avpermulation_color,
            hist_color=hist_color,
            thresholds_color=thresholds_color,
            hist_alpha=hist_alpha,
            edgecolor=edgecolor,
            bins=bins,
            title=title,
            legend_fontsize=legend_fontsize,
            textbox_fontsize=textbox_fontsize,
            axis_label_fontsize=axis_label_fontsize,
        )

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
        legend_fontsize=10,
        axis_label_fontsize=12,
    ):
        return plotting.plot_permulation_results_layered(
            self,
            test,
            fg_name,
            bg_name=bg_name,
            gaussfit_color=gaussfit_color,
            avpermulation_color=avpermulation_color,
            hist_color=hist_color,
            thresholds_color=thresholds_color,
            bins=bins,
            title=title,
            legend_fontsize=legend_fontsize,
            axis_label_fontsize=axis_label_fontsize,
        )

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

    def _emit_count_summary(self, file_obj, counts_dict):
        for label in ("loss_fg", "loss_bg", "dup_fg", "dup_bg"):
            _emit(f"  - {label}: {counts_dict[label]}\n", file_obj)

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
            flagged = self.true_odds.flagged_cells
            if flagged is not None and len(flagged):
                _emit(
                    f"Flagged cells: {len(flagged)} species-by-HOG cells in {len(self.true_odds.flagged_hog_idx)} "
                    f"HOGs treated as contaminated ({int(self.true_odds.missing_mat.sum())} left missing)",
                    file_obj,
                )
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
                f"User-defined significance threshold: {self.alpha} (two-sided; {self.alpha / 2} per direction)\n"
                f"Permulation-derived alpha threshold, loss: {self._fmt_ci(self.loss_ci_av, ndigits=3)}\n"
                f"Permulation-derived alpha threshold, duplication: {self._fmt_ci(self.dup_ci_av, ndigits=3)}\n\n"
                f"Significant HOGs (all): {self.all_hits_count}\n",
                file_obj,
            )
            self._emit_count_summary(file_obj, self.counts_hits)

            _emit("\nPre-specified directions (BH within each tail):", file_obj)
            for label, qcol in (
                (
                    "loss more likely in non-orb-weavers",
                    "Q-value loss more likely in bg (directional)",
                ),
                (
                    "duplication more likely in orb-weavers",
                    "Q-value duplication more likely in fg (directional)",
                ),
            ):
                q = self.results_df[qcol].dropna()
                _emit(
                    f"  - {label}: {len(q)} HOGs tested, {int((q < 0.05).sum())} with q < 0.05, "
                    f"smallest q {q.min():.4f} ({q.idxmin()})",
                    file_obj,
                )

            if self.species_of_interest is not None:
                sp_of_int_counts = self._count_species_of_interest_hits(
                    self.species_of_interest
                )
                _emit(
                    "Significant HOGs, {} present:\n".format(self.species_of_interest),
                    file_obj,
                )
                self._emit_count_summary(file_obj, sp_of_int_counts)

                _emit(
                    "Significant HOGs, using permulation p-values:\n",
                    file_obj,
                )
                perm_pval_counts = self._count_permulation_significant_hits()
                self._emit_count_summary(file_obj, perm_pval_counts)

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
        resolved_path = _resolve_repo_path(filepath)
        with open(resolved_path, "rb") as f:
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

        figure_outputs = [
            (self.plot_permulation_stats("loss", fg_name, bg_name)[0], "loss_stats_dists.png"),
            (self.plot_permulation_stats("dup", fg_name, bg_name)[0], "dup_stats_dists.png"),
            (self.plot_permulation_results("loss", fg_name, bg_name)[0], "loss_results.png"),
            (self.plot_permulation_results("dup", fg_name, bg_name)[0], "dup_results.png"),
        ]
        for figure, filename in figure_outputs:
            plotting.save_figure(figure, os.path.join(results_dir, filename))
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
    dir_suffix=None,
    flagged_cells_file=None,
    flag_strengths=("strong",),
):
    """Run the full odds ratio test.

    flagged_cells_file: optional table of suspected-contamination gene copies (see
    load_flagged_cells); flagged copies are removed and emptied cells are treated as missing in
    that HOG only. flag_strengths: which strength labels of a gene-level table to use.
    """

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

    flagged_cells = (
        load_flagged_cells(flagged_cells_file, strengths=flag_strengths)
        if flagged_cells_file is not None
        else None
    )

    # Get the true odds ratio results
    true_odds = OddsRatioResults(
        genecount_csv=genecount_csv,
        hog_node_genes_tsv=hog_node_genes_tsv,
        time=time_fmtd,
        foreground_list_filename=foreground_list_filename,
        background_list_filename=background_list_filename,
        buscos_filename=buscos_filename,
        flagged_cells=flagged_cells,
    )

    # Permulated tip values are read from the CSV written by permulations.R.
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
            dir_suffix=dir_suffix,
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
            bg_name=bg_name or "background",
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
            "loss", display_fg_name, display_bg_name
        )
        permulation_test_results.plot_permulation_stats(
            "dup", display_fg_name, display_bg_name
        )

        permulation_test_results.plot_permulation_results(
            "loss", display_fg_name, display_bg_name
        )
        permulation_test_results.plot_permulation_results(
            "dup", display_fg_name, display_bg_name
        )
        plt.show()

    return permulation_test_results


# Backward-compatible alias for old pickle files
PermutationTestResults = PermulationTestResults
