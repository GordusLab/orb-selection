"""
Odds ratio analysis of gene loss and gain: testing whether
the odds that genes in an orthogroup are lost or duplicated
differs significantly between a test group and a reference
group in a phylogeny
"""

import sys
import os
import random
import pickle
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, skew
from tqdm.auto import tqdm
import seaborn as sns
import orthogroup_gene_count

# Set the random seed for reproducibility
random.seed(42)

plt.rcParams['font.family'] = 'Verdana'

def _unique_results_dir(base_dir: str) -> str:
    """Return a unique directory path by appending _1, _2, ... if needed."""
    if not os.path.exists(base_dir):
        return base_dir

    counter = 1
    while True:
        candidate = f"{base_dir}_{counter}"
        if not os.path.exists(candidate):
            return candidate
        counter += 1

def drop_empty_cols(df, print_txt=True):
    """
    Drops columns where all entries (ignoring headers) are 0,
    to get rid of species not included in the current node/hierarchy.
    """    

    # Print the number of columns before cleaning
    num_columns_before = df.shape[1]
    if print_txt:
        print(f"Number of columns before dropping empty columns: {num_columns_before}")

    # Drop columns where all entries (ignoring headers) are 0
    df_cleaned = df.loc[:, (df.ne(0)).any(axis=0)]

    # Print the number of columns after cleaning
    num_columns_after = df_cleaned.shape[1]

    if print_txt:
        print(f"Number of columns after dropping empty columns: {num_columns_after}")
        print("Species with no sequences in any orthogroup have been dropped.")

    return df_cleaned


def occupancy_filter(arr, minimum, maximum, total_occ_arr):
    """
    Function to filter an array of values according to whether the
    orthogroup meets a certain occupancy threshold.

    Args:
        arr, minimum, maximum, total_occ_arr

    Returns:
        fltrd_arr (numpy array): occupancy array filtered according 
        to thresholds provided
    """    

    if maximum is None:
        maximum = total_occ_arr.max()

    idx = np.asarray((total_occ_arr >= minimum) & (total_occ_arr <= maximum)).nonzero()[
        0
    ]
    fltrd_arr = arr[idx]

    return fltrd_arr


def filter_for_sp_of_interest(df, genecount_df, species_name):
    """Filter the DataFrame for HOGs that include a species of interest"""

    print(f"Filtering for presence of {species_name}\n")

    sp_of_interest_present = genecount_df[genecount_df[species_name] != 0]

    sp_of_interest_present_hogs = sp_of_interest_present.index.values
    sp_of_interest_present_hogs = set(sp_of_interest_present_hogs)
    df_fltrd_sp_of_int = df[df.index.isin(sp_of_interest_present_hogs)]

    sp_of_interest_hogs_count = len(df_fltrd_sp_of_int)

    return df_fltrd_sp_of_int, sp_of_interest_hogs_count


def calculate_odds(
    foreground_bool_arr,
    background_bool_arr,
    test_bool_mat,
    busco_arr
):
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

    return odds_foreground_arr, odds_background_arr, log_odds_ratio_arr


class OddsRatioResults:
    """Class to hold the results of the odds ratio calculations"""

    def __init__(
        self,
        genecount_csv,
        hog_node_genes_tsv,
        test,
        time,
        foreground_list_filename,
        background_list_filename=None,
        buscos_filename=None
    ):
        """Initialize the inputs for the odds ratio calculations"""
        self.genecount_csv = genecount_csv
        self.hog_node_genes_tsv = hog_node_genes_tsv
        self.test = test
        self.time = time
        self.foreground_list_filename = foreground_list_filename
        self.background_list_filename = background_list_filename

        tests = ["duplication", "loss"]
        if test not in tests:
            raise ValueError("Invalid test. Expected one of: %s" % tests)

        self.foreground_list_arr = np.loadtxt(foreground_list_filename, dtype=str)

        if background_list_filename is not None:
            self.background_list_arr = np.loadtxt(background_list_filename, dtype=str)
        else:
            self.background_list_arr = None

        if buscos_filename is not None:
            self.buscos = pd.read_csv(buscos_filename, header=None)
            self.buscos.columns = ['Species', 'Single_copy_buscos', 'Duplicated_buscos','Total_buscos', 'Fraction_sc', 'Fraction_total']

        # Get the genecount arrays and matrices
        self._get_genecount_arrays()

        # Define the foreground and background arrays
        self._define_foreground()

        # Calculate the odds ratios and log odds ratios
        self.odds_foreground_arr, self.odds_background_arr, self.log_odds_ratio_arr = (
            calculate_odds(
                self.foreground_bool_arr,
                self.background_bool_arr,
                self.test_bool_mat,
                busco_arr=getattr(self, 'busco_arr', None)
            )
        )

        # Create the results DataFrame
        self.results_df = self.results_to_df()

    def __repr__(self):
        """String representation of the OddsRatioResults object"""
        return (
            f"OddsRatioResults(genecount_csv={self.genecount_csv}, "
            f"hog_node_genes_tsv={self.hog_node_genes_tsv}, "
            f"test={self.test}, foreground_list_filename={self.foreground_list_filename}, "
            f"background_list_filename={self.background_list_filename})"
        )

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
        if 'Total' in genecount_df.columns:
            genecount_df = genecount_df.drop(columns=['Total'])

        # Remove any empty columns (species with no genes in any HOGs)
        genecount_df = drop_empty_cols(genecount_df, print_txt=False)

        # Remove any species not in the foreground or background lists
        species_to_keep = self.foreground_list_arr
        if self.background_list_arr is not None:
            species_to_keep = np.concatenate(
                (species_to_keep, self.background_list_arr), axis=0
            )
        genecount_df = genecount_df[species_to_keep]

        # List of HOG IDs
        hog_list = genecount_df.index.values

        # List of all species
        all_species_arr = genecount_df.columns.to_numpy()

        # Get busco scores for each species
        if hasattr(self, 'buscos'):
            sc_busco_arr = self.buscos.set_index('Species').loc[all_species_arr]['Fraction_sc'].values
            total_busco_arr = self.buscos.set_index('Species').loc[all_species_arr]['Fraction_total'].values

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

        if self.test == "loss":
            self.test_bool_mat = loss_bool_mat

            # If testing for loss, the total fraction of complete buscos, whether single-copy or 
            # duplicated, should be a good proxy for how reliable loss calls will be
            if hasattr(self, 'buscos'):
                self.busco_arr = total_busco_arr

        elif self.test == "duplication":
            self.test_bool_mat = dup_bool_mat

            # If testing for duplication, the fraction of recovered buscos that are single-copy
            # should be a good proxy for how reliable duplication calls will be
            if hasattr(self, 'buscos'):
                self.busco_arr = sc_busco_arr

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

        print(
            f"{total_species_count} species total, {foreground_count} foreground, "
            f"{background_count} background"
        )

        self.total_species_count = total_species_count
        self.foreground_count = foreground_count
        self.background_count = background_count
        self.foreground_bool_arr = foreground_bool_arr
        self.background_bool_arr = background_bool_arr

    def results_to_df(self, occupancy_threshold=0, maximum=None):
        """Function to convert the results of the odds ratio calculations
        into a DataFrame for easier analysis and plotting"""

        # Create a DataFrame with the results
        results_df = pd.DataFrame(
            {
                "HOG": self.hog_list,
                "Occupancy": self.total_occ_arr.flatten(),
                f"Odds {self.test}, foreground": self.odds_foreground_arr.flatten(),
                f"Odds {self.test}, background": self.odds_background_arr.flatten(),
                "Log odds ratio": self.log_odds_ratio_arr.flatten(),
            }
        )
        results_df = results_df.set_index("HOG")

        # Filter the results based on the occupancy threshold
        if occupancy_threshold > 0:
            results_df = results_df[
                results_df["Occupancy"] >= occupancy_threshold
            ].copy()

        # If a maximum occupancy is specified, filter the results further
        if maximum is not None:
            results_df = results_df[
                (results_df["Occupancy"] >= occupancy_threshold)
                & (results_df["Occupancy"] <= maximum)
            ].copy()
        else:
            maximum = "none"

        return results_df

    def filter_for_bootstrap_hits(self, minimum, maximum, ci, alternative="two-tailed"):
        """Filter the DataFrame based on the occupancy and log odds ratio"""
        df = self.results_df

        print(
            f"Filtering log odds ratio results df for occupancy >= {minimum}, <= {maximum}\n"
            f"and log odds ratio threshold {ci}\n"
        )

        if alternative == "greater":
            df_fltrd = df[
                (df["Occupancy"] >= minimum)
                & (df["Occupancy"] <= maximum)
                & (df["Log odds ratio"] > ci[1])
            ]
        elif alternative == "less":
            df_fltrd = df[
                (df["Occupancy"] >= minimum)
                & (df["Occupancy"] <= maximum)
                & (df["Log odds ratio"] < ci[0])
            ]
        else:
            df_fltrd = df[
                (df["Occupancy"] >= minimum)
                & (df["Occupancy"] <= maximum)
                & ((df["Log odds ratio"] < ci[0]) | (df["Log odds ratio"] > ci[1]))
            ]

        total_hits = len(df_fltrd)

        return df_fltrd, total_hits


def define_foreground_random(species_incl_idx, total_species_count, foreground_count):
    """Function to randomly define the species designated as
    foreground and background (with the same number of each
    as the true number) for bootstrapping."""

    # Create new (empty for now) foreground and background arrays
    new_foregrounds = np.zeros(total_species_count)
    new_backgrounds = np.zeros(total_species_count)

    # Randomly select indices for the new foreground
    new_foregrounds_idx = np.random.choice(
        species_incl_idx, foreground_count, replace=False
    )
    new_foregrounds[new_foregrounds_idx] = 1

    # The rest of the species are background
    new_backgrounds_idx = np.setdiff1d(
        species_incl_idx, new_foregrounds_idx, assume_unique=True
    )
    new_backgrounds[new_backgrounds_idx] = 1

    return new_foregrounds, new_backgrounds


class BootstrapTestResults:
    """Class to perform the bootstrapping test for the odds ratio"""

    def __init__(
        self,
        true_odds,
        occupancy_threshold=0,
        maximum=None,
        alternative="less",
        a=0.05,
        bootstrap_reps=10000,
        species_of_interest=None,
    ):
        """Initialize the bootstrapping test for a given odds ratio results object"""

        self.true_odds = true_odds

        self.occupancy_threshold = occupancy_threshold

        alts = ["less", "greater"]
        if alternative not in alts:
            raise ValueError(
                "Invalid alternative hypothesis. Expected one of: %s" % alts
            )

        self.alternative = alternative
        self.bootstrap_reps = bootstrap_reps
        self.species_of_interest = species_of_interest

        if maximum is not None:
            self.maximum = maximum
        else:
            self.maximum = self.true_odds.total_species_count

        self.a = a
        self.z_crit = norm.ppf(self.a)

        # Filter the true log odds ratios based on the occupancy threshold
        self.true_fltrd_log_odds_ratios = occupancy_filter(
            self.true_odds.log_odds_ratio_arr,
            self.occupancy_threshold,
            self.maximum,
            self.true_odds.total_occ_arr,
        )

        # Calculate the skew, mean, and standard deviation of the true log odds ratios
        self.true_skew = skew(self.true_fltrd_log_odds_ratios)[0]
        self.true_mean = np.mean(self.true_fltrd_log_odds_ratios)
        self.true_stddev = np.std(self.true_fltrd_log_odds_ratios)

        # Run the bootstrapping test
        self._run_bootstrap()

        # Get the list of HOGs which are significantly different between
        # the test groups according to the bootstrapping thresholds
        self._get_hits_df()
        self._get_hits_hogs()

        # Print the results of the bootstrapping test
        self.print_bootstrap_results()

    def __repr__(self):
        """String representation of the BootstrapTestResults object"""
        return (
            f"BootstrapTestResults(true_odds={self.true_odds}, "
            #            f"stat={self.stat}, "
            f"occupancy_threshold={self.occupancy_threshold}, "
            f"maximum={self.maximum}, alternative={self.alternative}, "
            f"bootstrap_reps={self.bootstrap_reps})"
        )

    def _bootstrap_iter(self, i, counters, species_incl_idx):
        """Function to perform a single iteration of the bootstrapping test"""

        # Create a new set of foreground and background arrays for bootstrapping
        new_foregrounds, new_backgrounds = define_foreground_random(
            species_incl_idx,
            self.true_odds.total_species_count,
            self.true_odds.foreground_count,
        )

        # Recalculate the odds ratios and log odds ratios for the new foreground/background arrays
        new_log_odds_ratio = calculate_odds(
            new_foregrounds,
            new_backgrounds,
            self.true_odds.test_bool_mat,
            busco_arr=getattr(self, 'busco_arr', None)
        )[
            2
        ]  # only return the log odds ratio array

        new_log_odds_ratio_fltrd = occupancy_filter(
            new_log_odds_ratio,
            self.occupancy_threshold,
            self.maximum,
            self.true_odds.total_occ_arr,
        )

        # Calculate the statistics for the new bootstrapped distribution
        new_skew = skew(new_log_odds_ratio_fltrd)[0]
        new_mean = np.mean(new_log_odds_ratio_fltrd)
        new_stddev = np.std(new_log_odds_ratio_fltrd)

        if self.alternative == "greater":
            if new_skew > self.true_skew:
                counters["sk"] += 1
            if new_mean > self.true_mean:
                counters["mn"] += 1
            if new_stddev > self.true_stddev:
                counters["sd"] += 1

        else:
            if new_skew < self.true_skew:
                counters["sk"] += 1
            if new_mean < self.true_mean:
                counters["mn"] += 1
            if new_stddev < self.true_stddev:
                counters["sd"] += 1

        self.means[i] = new_mean
        self.stddevs[i] = new_stddev
        self.skews[i] = new_skew
        self.cis[i] = [
            new_mean + self.z_crit * new_stddev,  # the z-crit value is negative
            new_mean - self.z_crit * new_stddev,
        ]

        return counters

    def _run_bootstrap(self):
        """Creating 10,000 test log odds ratio distributions with the set
        of species defined as "foreground" chosen randomly (but still the
        same number of foreground  as the true number)"""



        print("\nLAUNCHING BOOTSTRAPPING TEST\n")

        if self.maximum != self.true_odds.total_species_count:
            print(f"** Maximum occupancy set to {self.maximum} **\n")
        else:
            pass

        if self.occupancy_threshold > 0:
            print(f"** Minimum occupancy set to {self.occupancy_threshold} **\n")

        if self.alternative == "greater":
            print(
                "Counting bootstrapped distributions in which statistics\n"
                "derived from randomly assigned test groups EXCEED\n"
                "the true distribution's statistics (right-tailed)...\n"
            )
        else:
            print(
                "Counting bootstrapped distributions in which statistics\n"
                "derived from randomly assigned test groups are SMALLER\n"
                "than the true distribution's statistics (left-tailed)...\n"
            )

        # Initialize lists to store the results of the bootstrapping
        self.cis = np.zeros((self.bootstrap_reps, 2))
        self.means = np.zeros(self.bootstrap_reps)
        self.stddevs = np.zeros(self.bootstrap_reps)
        self.skews = np.zeros(self.bootstrap_reps)

        # Initialize counters for how often the bootstrapped
        # distributions exceed the true values
        counters = {"mn": 0, "sd": 0, "sk": 0}

        # Get a list of the indices of the species being compared,
        # in case not all species are included in the analysis
        species_incl_arr = (
            self.true_odds.foreground_bool_arr + self.true_odds.background_bool_arr
        ).astype(bool)
        species_incl_idx = species_incl_arr.nonzero()[0]

        for i in tqdm(range(self.bootstrap_reps)):
            # Perform a single iteration of the bootstrapping test
            counters = self._bootstrap_iter(i, counters, species_incl_idx)

        p_vals = {
            "mn": counters["mn"] / 10000,
            "sd": counters["sd"] / 10000,
            "sk": counters["sk"] / 10000,
        }

        print(f"\nBootstrapping counter for MEAN: {str(counters['mn'])}")
        print(f"Bootstrapping counter for STD DEV: {str(counters['sd'])}")
        print(f"Bootstrapping counter for SKEW: {str(counters['sk'])}\n")

        self.p_values = p_vals

        # Getting average stats from across the 10000 bootstrapped distributions
        self.mean_av = np.mean(self.means)
        self.stddev_av = np.mean(self.stddevs)
        self.skew_av = np.mean(self.skews)
        self.ci_av = np.mean(self.cis, axis=0)

    def plot_bootstrap_stats(self, fg_name, bg_name="background"):
        """Plotting the bootstrapping means, standard deviations
        and alpha thresholds to ensure the results are relatively
        tightly distributed"""

        fig, axs = plt.subplots(1, 3, figsize=(15, 5))

        if self.alternative == "greater":
            alt = "right-tailed"
        else:
            alt = "left-tailed"

        # Making the text bold deletes spaces
        fg_name = fg_name.replace(" ", r"\ ")
        bg_name = bg_name.replace(" ", r"\ ")

        fig.suptitle(
            rf"$\bf{{Bootstrapped\ distribution\ stats\ for\ gene\ {self.true_odds.test}, {fg_name}\ vs. {bg_name}}}$" + "\n"
            f"Maximum occupancy = {self.maximum}, "
            f"minimum occupancy = {self.occupancy_threshold}, "
            f"{alt}",
            fontsize=16,
        )

        sns.histplot(
            data=self.skews,
            kde=True,
            bins=50,
            stat="density",
            ax=axs[0],
            legend=False,
        )

        axs[0].set_title("Skews")
        axs[0].set(xlabel="skew", ylabel="Density")
        axs[0].axvline(x=self.skew_av, linestyle="dotted")

        sns.histplot(
            data=self.means,
            kde=True,
            bins=50,
            stat="density",
            ax=axs[1],
            legend=False,
        )
        axs[1].set_title("Means")
        axs[1].set(xlabel="mean", ylabel="Density")
        axs[1].axvline(x=self.mean_av, linestyle="dotted")

        sns.histplot(
            data=self.stddevs,
            kde=True,
            bins=50,
            stat="density",
            ax=axs[2],
            legend=False,
        )

        axs[2].set_title("Standard deviations")
        axs[2].set(xlabel="stddev", ylabel="Density")
        axs[2].axvline(x=self.stddev_av, linestyle="dotted")

        plt.tight_layout()

        return fig, axs

    def plot_bootstrap_results(
        self,
        fg_name,
        bg_name="background",
        gaussfit_color="blue",
        avbootstrap_color="red",
        hist_color="red",
        thresholds_color="darkred",
        bins=100,
        title=True
    ):
        """Function to plot the results of the bootstrapping test"""
        
        # Making the text bold deletes spaces
        fg_name = fg_name.replace(" ", r"\ ")
        bg_name = bg_name.replace(" ", r"\ ")

        fig, ax = plt.subplots(figsize=(6, 5))

        if title:
            fig.suptitle(
                rf"$\bf{{Log\ odds\ ratio\ of\ gene\ {self.true_odds.test}, {fg_name}\ vs. {bg_name}}}$" + "\n"
                f"Maximum occupancy = {self.maximum}, "
                f"minimum occupancy = {self.occupancy_threshold}",
                fontsize=14,
            )

        # Histogram of the true log odds ratios, filtered for occupancy
        ax.hist(
            self.true_fltrd_log_odds_ratios,
            bins=bins,
            density=True,
            color=hist_color,
            alpha=0.3,
            edgecolor=hist_color
        )

        print(hist_color)

        x = np.linspace(
            self.true_fltrd_log_odds_ratios.min(),
            self.true_fltrd_log_odds_ratios.max(),
            100,
        )

        # Gaussian fit to the true log odds ratios
        ax.plot(
            x,
            norm.pdf(
                x,
                np.mean(self.true_fltrd_log_odds_ratios),
                np.std(self.true_fltrd_log_odds_ratios),
            ),
            color=gaussfit_color,
            linestyle="--",
            label="Gaussian fit to\ntrue distribution",
        )

        # Normal distribution using the average bootstrapped stats
        ax.plot(
            x,
            norm.pdf(x, self.mean_av, self.stddev_av),
            color=avbootstrap_color,
            linestyle="--",
            label="Average BS'd\ndistribution",
        )

        # Bootstrap-derived confidence intervals
        ax.axvline(
            x=self.ci_av[0],
            label=f"Mean BS'd\nthresholds for\nalpha={self.a}",
            linestyle="dotted",
            color=thresholds_color
        )

        ax.axvline(
            x=self.ci_av[1],
            linestyle="dotted",
            color=thresholds_color
        )

        ax.text(
            0.35,
            0.95,
            f"BS'd mean = {self.mean_av:.2f}\n"
            f"True mean = {self.true_mean:.2f}\n\n"
            f"BS'd std. dev. = {self.stddev_av:.2f}\n"
            f"True std. dev. = {self.true_stddev:.2f}\n\n"
            f"BS'd skew = {self.skew_av:.2f}\n"
            f"True skew = {self.true_skew:.2f}\n",
            transform=ax.transAxes,
            fontsize=10,
            ha="left",
            va="top",
        )

        plt.xlabel("Log odds ratio", fontsize=12, fontweight="bold")
        plt.ylabel("Density", fontsize=12, fontweight="bold")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        plt.legend(fontsize=10, loc="upper right")
        plt.tight_layout()

        return fig, ax

    def plot_bootstrap_results_layered(
        self,
        fg_name,
        bg_name="background",
        gaussfit_color="blue",
        avbootstrap_color="red",
        hist_color="red",
        thresholds_color="darkred",
        bins=100,
        title=True
    ):
        """Function to create 4 sequential plots with elements layered on top of each other.
        
        Returns 4 figures showing progressive buildup:
        1. Average bootstrapped distribution
        2. + Bootstrap-derived thresholds (with BS'd stats)
        3. + Histogram of true log odds ratios (with BS'd stats)
        4. + Gaussian fit to the histogram (with True and BS'd stats)
        """
        
        # Making the text bold deletes spaces
        fg_name = fg_name.replace(" ", r"\ ")
        bg_name = bg_name.replace(" ", r"\ ")

        # Define x-axis range
        x = np.linspace(
            self.true_fltrd_log_odds_ratios.min(),
            self.true_fltrd_log_odds_ratios.max(),
            100,
        )

        # ----------- Compute a common y-limit so all panels share the same scale -----------
        avg_pdf = norm.pdf(x, self.mean_av, self.stddev_av)
        true_pdf = norm.pdf(
            x,
            np.mean(self.true_fltrd_log_odds_ratios),
            np.std(self.true_fltrd_log_odds_ratios),
        )
        hist_vals, _ = np.histogram(
            self.true_fltrd_log_odds_ratios, bins=bins, density=True
        )
        hist_max = hist_vals.max() if hist_vals.size else 0.0
        y_max = max(avg_pdf.max(), true_pdf.max(), hist_max) * 1.05
        if y_max == 0:
            y_max = 1  # fallback to avoid zero-height axis
        # -------------------------------------------------------------------------------

        title_str = (
            rf"$\bf{{Log\ odds\ ratio\ of\ gene\ {self.true_odds.test}, {fg_name}\ vs. {bg_name}}}$" + "\n"
            f"Maximum occupancy = {self.maximum}, "
            f"minimum occupancy = {self.occupancy_threshold}"
        )

        figs = []
        axes = []

        # ========== PLOT 1: Average bootstrapped distribution ==========
        fig1, ax1 = plt.subplots(figsize=(8, 6))
        
        if title:
            fig1.suptitle(title_str, fontsize=14)

        ax1.plot(
            x,
            norm.pdf(x, self.mean_av, self.stddev_av),
            color=avbootstrap_color,
            linewidth=2.5,
            label="Average BS'd\ndistribution",
        )
        ax1.fill_between(
            x,
            norm.pdf(x, self.mean_av, self.stddev_av),
            alpha=0.2,
            color=avbootstrap_color,
        )
        
        ax1.set_xlabel("Log odds ratio", fontsize=14, fontweight="bold")
        ax1.set_ylabel("Density", fontsize=14, fontweight="bold")
        ax1.set_ylim(bottom=0, top=y_max)
        ax1.set_xlim(x.min(), x.max())
        plt.setp(ax1.get_xticklabels(), fontsize=13)
        plt.setp(ax1.get_yticklabels(), fontsize=13)
        ax1.legend(fontsize=13, loc="upper right", ncol=1, labelspacing=0.8, handlelength=1.5)
        plt.tight_layout()
        figs.append(fig1)
        axes.append(ax1)

        # ========== PLOT 2: + Bootstrap-derived thresholds ==========
        fig2, ax2 = plt.subplots(figsize=(8, 6))
        
        if title:
            fig2.suptitle(title_str, fontsize=14)

        ax2.plot(
            x,
            norm.pdf(x, self.mean_av, self.stddev_av),
            color=avbootstrap_color,
            linewidth=2.5,
            label="Average BS'd\ndistribution",
        )
        ax2.fill_between(
            x,
            norm.pdf(x, self.mean_av, self.stddev_av),
            alpha=0.2,
            color=avbootstrap_color,
            zorder=0,
        )

        ax2.axvline(
            x=self.ci_av[0],
            label=f"Mean BS'd\nthresholds for\nalpha={self.a}",
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )
        ax2.axvline(
            x=self.ci_av[1],
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )

        ax2.text(
            0.03,
            0.95,
            f"BS'd mean = {self.mean_av:.2f}\n"
            f"BS'd std. dev. = {self.stddev_av:.2f}\n"
            f"BS'd skew = {self.skew_av:.2f}\n",
            transform=ax2.transAxes,
            fontsize=12,
            ha="left",
            va="top",
        )

        ax2.set_xlabel("Log odds ratio", fontsize=14, fontweight="bold")
        ax2.set_ylabel("Density", fontsize=14, fontweight="bold")
        ax2.set_ylim(bottom=0, top=y_max)
        ax2.set_xlim(x.min(), x.max())
        plt.setp(ax2.get_xticklabels(), fontsize=13)
        plt.setp(ax2.get_yticklabels(), fontsize=13)
        ax2.legend(fontsize=13, loc="upper right", ncol=1, labelspacing=0.8, handlelength=1.5)
        plt.tight_layout()
        figs.append(fig2)
        axes.append(ax2)

        # ========== PLOT 3: + Histogram ==========
        fig3, ax3 = plt.subplots(figsize=(8, 6))
        
        if title:
            fig3.suptitle(title_str, fontsize=14)

        ax3.plot(
            x,
            norm.pdf(x, self.mean_av, self.stddev_av),
            color=avbootstrap_color,
            linewidth=2.5,
            label="Average BS'd\ndistribution",
        )
        ax3.fill_between(
            x,
            norm.pdf(x, self.mean_av, self.stddev_av),
            alpha=0.2,
            color=avbootstrap_color,
            zorder=0,
        )

        ax3.axvline(
            x=self.ci_av[0],
            label=f"Mean BS'd\nthresholds for\nalpha={self.a}",
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )
        ax3.axvline(
            x=self.ci_av[1],
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )

        ax3.hist(
            self.true_fltrd_log_odds_ratios,
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
            f"BS'd mean = {self.mean_av:.2f}\n"
            f"BS'd std. dev. = {self.stddev_av:.2f}\n"
            f"BS'd skew = {self.skew_av:.2f}\n",
            transform=ax3.transAxes,
            fontsize=12,
            ha="left",
            va="top",
        )

        ax3.set_xlabel("Log odds ratio", fontsize=14, fontweight="bold")
        ax3.set_ylabel("Density", fontsize=14, fontweight="bold")
        ax3.set_xlim(x.min(), x.max())
        ax3.set_ylim(bottom=0, top=y_max)
        plt.setp(ax3.get_xticklabels(), fontsize=13)
        plt.setp(ax3.get_yticklabels(), fontsize=13)
        ax3.legend(fontsize=13, loc="upper right", ncol=1, labelspacing=0.8, handlelength=1.5)
        plt.tight_layout()
        figs.append(fig3)
        axes.append(ax3)

        # ========== PLOT 4: + Gaussian fit ==========
        fig4, ax4 = plt.subplots(figsize=(8, 6))
        
        if title:
            fig4.suptitle(title_str, fontsize=14)

        ax4.plot(
            x,
            norm.pdf(x, self.mean_av, self.stddev_av),
            color=avbootstrap_color,
            linewidth=2.5,
            label="Average BS'd\ndistribution",
        )
        ax4.fill_between(
            x,
            norm.pdf(x, self.mean_av, self.stddev_av),
            alpha=0.2,
            color=avbootstrap_color,
            zorder=0,
        )

        ax4.axvline(
            x=self.ci_av[0],
            label=f"Mean BS'd\nthresholds for\nalpha={self.a}",
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )
        ax4.axvline(
            x=self.ci_av[1],
            linestyle="dotted",
            color=thresholds_color,
            linewidth=2,
        )

        ax4.hist(
            self.true_fltrd_log_odds_ratios,
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
                np.mean(self.true_fltrd_log_odds_ratios),
                np.std(self.true_fltrd_log_odds_ratios),
            ),
            color=gaussfit_color,
            linestyle="--",
            linewidth=2,
            label="Gaussian fit to\ntrue distribution",
        )

        ax4.text(
            0.03,
            0.95,
            f"BS'd mean = {self.mean_av:.2f}\n"
            f"BS'd std. dev. = {self.stddev_av:.2f}\n"
            f"BS'd skew = {self.skew_av:.2f}\n\n"
            f"True mean = {self.true_mean:.2f}\n"
            f"True std. dev. = {self.true_stddev:.2f}\n"
            f"True skew = {self.true_skew:.2f}\n",
            transform=ax4.transAxes,
            fontsize=12,
            ha="left",
            va="top",
        )

        ax4.set_xlabel("Log odds ratio", fontsize=14, fontweight="bold")
        ax4.set_ylabel("Density", fontsize=14, fontweight="bold")
        ax4.set_xlim(x.min(), x.max())
        ax4.set_ylim(bottom=0, top=y_max)
        plt.setp(ax4.get_xticklabels(), fontsize=13)
        plt.setp(ax4.get_yticklabels(), fontsize=13)
        ax4.legend(fontsize=13, loc="upper right", ncol=1, labelspacing=0.8, handlelength=1.5)
        plt.tight_layout()
        figs.append(fig4)
        axes.append(ax4)

        return figs, axes

    def _get_hits_df(self):
        """Filter the dataframe for HOGs that have LORs exceeding the
        alpha values determined by bootstrapping, and filter for those
        which include a species of interest if specified."""

        # Filter the DataFrame for significant hits
        df_fltrd, total_hits = self.true_odds.filter_for_bootstrap_hits(
            self.occupancy_threshold, self.maximum, self.ci_av
        )

        if self.species_of_interest is not None:
            # Filter out the hits that do not include the
            # species of interest if specified
            df_fltrd_sp_of_int, sp_of_int_hits = filter_for_sp_of_interest(
                df_fltrd, self.true_odds.genecount_df, self.species_of_interest
            )

            self.results_fltrd_df = df_fltrd_sp_of_int
            self.all_hits_count = total_hits
            self.sp_of_int_hits_count = sp_of_int_hits

        else:
            self.results_fltrd_df = df_fltrd
            self.all_hits_count = total_hits
            self.sp_of_int_hits_count = None

    def _get_hits_hogs(self):
        """Get list of HOGs which are significant according to
        the bootstrapping test"""

        self.hits_hogs_list = self.results_fltrd_df.index.tolist()

    def print_bootstrap_results(self, fname=sys.stdout):
        """Function to print the results of the bootstrapping test"""

        if fname is not sys.stdout:
            fname = open(fname, 'w')

        if self.alternative == "greater":
            alt = "right-tailed"
        else:
            alt = "left-tailed"

        if self.maximum == self.true_odds.total_species_count:
            maximum = "no"
        else:
            maximum = self.maximum


        print(
            "*********************** RESULTS ***********************\n\n"
            f"Bootstrapping test with {self.bootstrap_reps} repetitions "
            f"for {self.true_odds.test} ({alt})\n"
            f"with minimum occupancy *{self.occupancy_threshold}* and "
            f"maximum occupancy *{maximum}* \n"
            f"Analysis run on {self.true_odds.time}"
        , file=fname)

        print(f"Foreground list: {self.true_odds.foreground_list_filename}", file=fname)
        if self.true_odds.background_list_filename is not None:
            print(f"Background list: {self.true_odds.background_list_filename}", file=fname)
        print(f"Gene count file: {self.true_odds.genecount_csv}", file=fname)
        print(f"Hierarchical orthogroup file: {self.true_odds.hog_node_genes_tsv}\n", file=fname)

        print(
            f"Total number of HOGs in node: {len(self.true_odds.hog_list)}\n"
            f"Total number of HOGs within occupancy threshold: "
            f"{len(self.true_fltrd_log_odds_ratios)}\n"
            f"Total species: {self.true_odds.total_species_count}\n"
            f"Foreground count: {self.true_odds.foreground_count}\n"
            f"Background count: {self.true_odds.background_count}\n"
            f"True mean: {self.true_mean}\n"
            f"True standard deviation: {self.true_stddev}\n"
            f"True skew: {self.true_skew}\n"
        , file=fname)

        print(
            "** BOOTSTRAPPING P-VALUES ** \n\n"
            f"Probability that the null is true for MEAN: {self.p_values['mn']}\n"
            f"Probability that the null is true for STANDARD DEVIATION: {self.p_values['sd']}\n"
            f"Probability that the null is true for SKEW: {self.p_values['sk']}\n"
        , file=fname)

        print(
            f"Bootstrapped average mean: {self.mean_av}\n"
            f"Bootstrapped average standard deviation: {self.stddev_av}\n"
            f"Bootstrapped average skew: {self.skew_av}\n"
            f"User-defined significance threshold: {self.a}\n"
            f"Bootstrap-derived alpha threshold: {self.ci_av}\n\n"
            "Total HOGs with significantly different LORs between\n"
            f"foreground and background (two-tailed): {self.all_hits_count}\n"
        , file=fname)

        if self.species_of_interest is not None:
            print(
                "Total HOGs with significantly different LORs between\n"
                f"foreground and background (two-tailed, {self.species_of_interest} "
                f"present): {self.sp_of_int_hits_count}\n"
            , file=fname)

    def save_pickle_file(self, fname):
        """
        Saves the bootstrap results object to a pickle file.
        """

        with open(fname, 'wb') as file:
            pickle.dump(self, file)

    @classmethod
    def load_from_pickle(cls, filepath: str) -> 'BootstrapTestResults':
        """Load results object from a pickle file."""
        with open(filepath, 'rb') as f:
            return pickle.load(f)

    def save_results_files(self, results_dir, fg_name, bg_name="background"):
        """
        Takes in bootstrap test results instance and saves all
        relevant plots and tables.
        """

        # Save the full table of odds and odds ratios to a csv
        self.true_odds.results_df.to_csv(
            f"{results_dir}/{self.true_odds.test}_odds_results_all.csv",
            index=True
        )

        filename = (
            f"{results_dir}/{self.true_odds.test}"
            f"_occ{self.occupancy_threshold}-"
            f"{self.maximum}"
        )

        if self.alternative == "greater":
            alt = "_RT"
        else:
            alt = "_LT"

        # Save the table filtered to HOGs considered "hits"
        if self.species_of_interest is not None:
            self.results_fltrd_df.to_csv(
                (
                    filename +
                    f"{alt}" +
                    f"_{self.species_of_interest}" +
                    "_fltrd_bootstrap_hits.csv"
                    ),
                index=True
            )
        else:
            self.results_fltrd_df.to_csv(
                (
                    filename + 
                    f"{alt}" +
                    "_fltrd_bootstrap_hits.csv"
                    ),
                index=True
            )

        # Save a text file summarizing results from the analysis
        self.print_bootstrap_results(
            fname = filename +
                f"{alt}" +
                "_results_summary.txt"
            )
        
        # Save the bootstrapping results object as a pickle file
        self.save_pickle_file(
            fname = filename +
            f"{alt}" +
            ".pkl"
        )

        #### Save figures ####

        # Bootstrapping statistics
        fig, _ = self.plot_bootstrap_stats(fg_name, bg_name)

        fig.savefig(
            (
                filename + 
                f"{alt}" +
                "_bootstrap_stats_dists.png"
                ),
            dpi=300
        )

        # True distribution overlaid with average bootstrapped distribution
        fig, _ = self.plot_bootstrap_results(fg_name, bg_name)

        fig.savefig(
            (
                filename + 
                f"{alt}" +
                "_bootstrap_results.png"
                ),
            dpi=300
        )

        print(
            f"Results files saved to {results_dir}\n\n"
            "Files include: \n"
            "\t 1. [test]_bootstrap_results.png: True LORs distribution\n" 
            "\t\tvs. average bootstrapped distribution\n"
            "\t 2. [test]_bootstrap_stats_dists.png: Histograms of the\n"
            "\t\tmeans, standard deviations, and skews of all 10,000\n"
            "\t\tbootstrapped LOR distributions\n"
            "\t 3. [test]_results_summary.txt: Text file summarizing results\n"
            "\t 4. [test].pkl: Pickle file storing the BoostrapTestResults\n"
            "\t\tinstance generated by this analysis\n"
            "\t 5. [test]_[species]_fltrd_bootstrap_hits.csv: Results table \n"
            "\t\tfiltered for occupancy, species of interest (if specified),\n"
            "\t\tand surpassing bootstrapping-derived significance thresholds\n"
            "\t 6. All odds and log odds ratios (not filtered for occupancy,\n"
            "\t\tspecies, or significance)"
            )

def odds_ratio_test(
    foreground_list_filename,
    hog_node_genes_tsv,
    test,
    genecount_csv=None,
    occupancy_threshold=0,
    max_occ=None,
    alternative="less",  # or "greater"
    alpha=0.05,
    bootstrap_reps=10000,
    background_list_filename=None,
    species_of_interest=None,
    results_dir=None,
    fg_name=None,
    bg_name=None,
    buscos_filename=None
):
    """Run the full odds ratio test"""

    if results_dir is not None: 
        if fg_name is None: 
            raise ValueError(
                "Please provide a descriptive name for your test group, \n"
                "e.g. fg_name = 'orbweavers'. "
                "You may also name your background, e.g.\n"
                "bg_name = 'non-orbweavers'")

    time = datetime.now()
    time_fmtd = time.strftime("%Y-%m-%d at %H:%M:%S")
    date_short = time.strftime("%b%d")  # e.g., "Jan29"

    # Generate the genecount file if not provided
    if genecount_csv is None:
        genecount_csv = orthogroup_gene_count.main(hog_node_genes_tsv)

    # Get the true odds ratio results
    true_odds = OddsRatioResults(
        genecount_csv=genecount_csv,
        hog_node_genes_tsv=hog_node_genes_tsv,
        test=test,
        time=time_fmtd,
        foreground_list_filename=foreground_list_filename,
        background_list_filename=background_list_filename, 
        buscos_filename=buscos_filename
    )

    # Run bootstrapping test
    bootstrap_test_results = BootstrapTestResults(
        true_odds=true_odds,
        occupancy_threshold=occupancy_threshold,
        maximum=max_occ,
        alternative=alternative,
        a=alpha,
        bootstrap_reps=bootstrap_reps,
        species_of_interest=species_of_interest,
    )

    bootstrap_test_results.print_bootstrap_results()

    if results_dir is not None:
        base_results_dir = f"{results_dir}/Results_{date_short}"
        unique_results_dir = _unique_results_dir(base_results_dir)
        os.makedirs(unique_results_dir, exist_ok=True)
        bootstrap_test_results.save_results_files(
            results_dir=unique_results_dir,
            fg_name=fg_name,
            bg_name=bg_name
        )

    return bootstrap_test_results
