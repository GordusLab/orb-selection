"""
Odds ratio analysis of gene loss and gain: testing whether
the odds that genes in an orthogroup are missing, in single copy,
or in multiple copies differs significantly between orb-weavers
and non-orb-weavers"""

import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, skew  # , skewnorm
from tqdm.auto import tqdm
import seaborn as sns
import src.id_converter as id_converter
from src.orthogroup_filter import drop_empty_cols


def occupancy_filter(arr, minimum, maximum, total_occ_arr):
    """Function to filter an array of values according to whether the
    orthogroup meets a certain occupancy threshold."""

    if maximum is None:
        maximum = total_occ_arr.max()

    idx = np.asarray((total_occ_arr >= minimum) & (total_occ_arr <= maximum)).nonzero()[
        0
    ]
    fltrd_arr = arr[idx]

    return fltrd_arr


def filter_for_udiv(df, genecount_df):
    """Filter the DataFrame for HOGs that include U. diversus"""

    udiv_present = genecount_df[genecount_df["Uloborus_diversus"] != 0]

    udiv_present_hogs = udiv_present.index.values
    udiv_present_hogs = set(udiv_present_hogs)
    df_ud_fltrd = df[df.index.isin(udiv_present_hogs)]

    udiv_hogs_count = len(df_ud_fltrd)

    return df_ud_fltrd, udiv_hogs_count


def calculate_odds(
    orb_bool_arr, nonorb_bool_arr, test_bool_mat, orb_count, nonorb_count
):
    """Function to calculate the odds ratio and log odds ratio"""

    # I don't know why this is necessary but my kernel crashes without it
    orb_bool_arr = orb_bool_arr.reshape(orb_bool_arr.size, 1)
    nonorb_bool_arr = nonorb_bool_arr.reshape(nonorb_bool_arr.size, 1)

    # Calculate the number of orb-weavers and non-orb-weavers that are
    # [present, missing, duplicated, single copy]

    ## orb-weavers yes
    orb_yes_arr = np.matmul(test_bool_mat, orb_bool_arr)

    ## non-orbweavers yes
    nonorb_yes_arr = np.matmul(test_bool_mat, nonorb_bool_arr)

    ## orb-weavers no
    orb_no_arr = orb_count - orb_yes_arr

    ## non-orb-weavers no
    nonorb_no_arr = nonorb_count - nonorb_yes_arr

    # Use the Haldane-Anscombe correction to account for any 0 entries
    # when calculating the odds ratios
    orb_yes_arr += 0.5
    nonorb_yes_arr += 0.5
    orb_no_arr += 0.5
    nonorb_no_arr += 0.5

    # Calculate odds missing for orb & non-orb, odds ratio, and log odds ratio

    ## odds
    odds_orb_arr = orb_yes_arr / orb_no_arr
    odds_nonorb_arr = nonorb_yes_arr / nonorb_no_arr

    ## odds ratio
    odds_ratio_arr = odds_orb_arr / odds_nonorb_arr

    ## log odds ratio
    log_odds_ratio_arr = np.log(odds_ratio_arr)

    return odds_orb_arr, odds_nonorb_arr, log_odds_ratio_arr


class OddsRatioResults:
    """Class to hold the results of the odds ratio calculations"""

    def __init__(
        self,
        genecount_csv,
        hog_node_genes_tsv,
        test,
        orb_list_filename,
        nonorb_list_filename=None,
    ):
        """Initialize the inputs for the odds ratio calculations"""
        self.genecount_csv = genecount_csv
        self.hog_node_genes_tsv = hog_node_genes_tsv
        self.test = test
        self.orb_list_filename = orb_list_filename
        self.nonorb_list_filename = nonorb_list_filename

        tests = ["duplication", "loss"]
        if test not in tests:
            raise ValueError("Invalid test. Expected one of: %s" % tests)

        self.orb_list_arr = np.loadtxt(orb_list_filename, dtype=str)

        if nonorb_list_filename is not None:
            self.nonorb_list_arr = np.loadtxt(nonorb_list_filename, dtype=str)
        else:
            self.nonorb_list_arr = None

        # Get the genecount arrays and matrices
        self._get_genecount_arrays()

        # Define the orb and non-orb arrays
        self._define_orb()

        # Calculate the odds ratios and log odds ratios
        self.odds_orb_arr, self.odds_nonorb_arr, self.log_odds_ratio_arr = (
            calculate_odds(
                self.orb_bool_arr,
                self.nonorb_bool_arr,
                self.test_bool_mat,
                self.orb_count,
                self.nonorb_count,
            )
        )

        # Create the results DataFrame
        self.results_df = self.results_to_df()

    def __repr__(self):
        """String representation of the OddsRatioResults object"""
        return (
            f"OddsRatioResults(genecount_csv={self.genecount_csv}, "
            f"hog_node_genes_tsv={self.hog_node_genes_tsv}, "
            f"test={self.test}, orb_list_filename={self.orb_list_filename}, "
            f"nonorb_list_filename={self.nonorb_list_filename})"
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

        # Remove any species not in this node
        genecount_df = drop_empty_cols(genecount_df)

        # List of HOG IDs
        hog_list = genecount_df.index.values

        # List of all spiders
        all_species_arr = genecount_df.columns.to_numpy()

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
        elif self.test == "duplication":
            self.test_bool_mat = dup_bool_mat

    def _define_orb(self):
        """Function to make a numpy vector corresponding to whether
        each species is an orb-weaver or not: 1 = orb, 0 = non-orb.
        If the non-orbweavers are not specified, any animal that is
        not specified to be OW is considered non-OW."""

        # assign True to index of orbweavers in the spiders array
        orb_bool_arr = np.isin(self.all_species_arr, self.orb_list_arr)

        total_spider_count = len(orb_bool_arr)
        orb_count = orb_bool_arr.sum()

        # If no background is specified, all spiders not identified as
        # orb-weavers are considered non-orb-weavers
        if self.nonorb_list_arr is None:
            nonorb_bool_arr = np.isin(
                self.all_species_arr, self.orb_list_arr, invert=True
            )

        # If the non-orb-weavers are specified
        else:
            nonorb_bool_arr = np.isin(self.all_species_arr, self.nonorb_list_arr)

        nonorb_count = nonorb_bool_arr.sum()

        nonorb_bool_arr = nonorb_bool_arr.reshape(nonorb_bool_arr.size, 1).astype(float)
        orb_bool_arr = orb_bool_arr.reshape(orb_bool_arr.size, 1).astype(float)

        print(
            f"{total_spider_count} spiders total, {orb_count} orb-weavers, "
            f"{nonorb_count} non-orb-weavers"
        )

        self.total_spider_count = total_spider_count
        self.orb_count = orb_count
        self.nonorb_count = nonorb_count
        self.orb_bool_arr = orb_bool_arr
        self.nonorb_bool_arr = nonorb_bool_arr

    def results_to_df(self, occupancy_threshold=0, maximum=None):
        """Function to convert the results of the odds ratio calculations
        into a DataFrame for easier analysis and plotting"""

        # Create a DataFrame with the results
        results_df = pd.DataFrame(
            {
                "HOG": self.hog_list,
                "Occupancy": self.total_occ_arr.flatten(),
                f"Odds {self.test}, orb": self.odds_orb_arr.flatten(),
                f"Odds {self.test}, non-orb": self.odds_nonorb_arr.flatten(),
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

    def plot_log_odds_ratios(self, occupancy_threshold=0, maximum=None):
        """Function to plot the log odds ratio results for a certain test,
        given occupancy thresholds (no bootstrapping)"""

        fig, ax = plt.subplots()

        fig.suptitle(
            f"Log odds ratio of {self.test}, orb:non-orb\n"
            f"Maximum occupancy = {maximum}, "
            f"minimum occupancy = {occupancy_threshold}",
            fontsize=16,
        )

        fltrd_log_odds_ratios = occupancy_filter(
            self.log_odds_ratio_arr, occupancy_threshold, maximum, self.total_occ_arr
        )

        # Plotting the LORs
        sns.histplot(
            data=fltrd_log_odds_ratios,
            kde=False,
            bins=100,
            stat="density",
            ax=ax,
            legend=False,
        )

        ax.set(xlabel="Log odds ratio", ylabel="Density")

        plt.tight_layout()

        plt.show()

        return fig, ax

    def filter_for_bootstrap_hits(
        self, minimum, maximum, alpha, alternative="two-tailed"
    ):
        """Filter the DataFrame based on the occupancy and log odds ratio"""
        alpha = float(np.abs(alpha))
        df = self.results_df

        print(
            f"Filtering for occupancy >= {minimum}, <= {maximum}, "
            f"and log odds ratio threshold {alpha}"
        )

        if alternative == "greater":
            df_fltrd = df[
                (df["Occupancy"] >= minimum)
                & (df["Occupancy"] <= maximum)
                & (df["Log odds ratio"] > alpha)
            ]
        elif alternative == "less":
            df_fltrd = df[
                (df["Occupancy"] >= minimum)
                & (df["Occupancy"] <= maximum)
                & (df["Log odds ratio"] < -alpha)
            ]
        else:
            df_fltrd = df[
                (df["Occupancy"] >= minimum)
                & (df["Occupancy"] <= maximum)
                & ((df["Log odds ratio"] < -alpha) | (df["Log odds ratio"] > alpha))
            ]

        total_hits = len(df_fltrd)

        return df_fltrd, total_hits


def define_orb_random(species_incl_idx, total_spider_count, orb_count):
    """Function to randomly define the spiders designated as
    orb-weaver and non-orb-weaver (with the same number of each
    as the true number) for bootstrapping."""

    # Create new (empty for now) orb and non-orb arrays
    new_orbs = np.zeros(total_spider_count)
    new_nonorbs = np.zeros(total_spider_count)

    # Randomly select indices for the new orb-weavers
    new_orbs_idx = np.random.choice(species_incl_idx, orb_count, replace=False)
    new_orbs[new_orbs_idx] = 1

    # The rest of the species are non-orb-weavers
    new_nonorbs_idx = np.setdiff1d(species_incl_idx, new_orbs_idx, assume_unique=True)
    new_nonorbs[new_nonorbs_idx] = 1

    return new_orbs, new_nonorbs


class BootstrapTestResults:
    """Class to perform the bootstrapping test for the odds ratio"""

    def __init__(
        self,
        true_odds,
        stat,
        occupancy_threshold=0,
        maximum=None,
        alternative="less",
        a=0.05,
        bootstrap_reps=10000,
    ):
        """Initialize the bootstrapping test for a given odds ratio results object"""

        self.true_odds = true_odds

        stats = ["mean", "stddev", "skew"]
        if stat not in stats:
            raise ValueError("Invalid bootstrapping stat. Expected one of: %s" % stats)

        self.stat = stat
        self.occupancy_threshold = occupancy_threshold
        self.alternative = alternative
        self.bootstrap_reps = bootstrap_reps

        if maximum is not None:
            self.maximum = maximum
        else:
            self.maximum = self.true_odds.total_spider_count

        # Determine the p-value threshold based on the alternative hypothesis
        if self.alternative == "greater":
            self.a = a
            self.z_crit = norm.ppf(self.a)
        else:
            self.a = 1-a
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

    def __repr__(self):
        """String representation of the BootstrapTestResults object"""
        return (
            f"BootstrapTestResults(true_odds={self.true_odds}, "
            f"stat={self.stat}, occupancy_threshold={self.occupancy_threshold}, "
            f"maximum={self.maximum}, alternative={self.alternative}, "
            f"bootstrap_reps={self.bootstrap_reps})"
        )

    def _bootstrap_iter(self, i, counter, species_incl_idx):
        """Function to perform a single iteration of the bootstrapping test"""

        # Create a new set of orb and non-orb arrays for bootstrapping
        new_orbs, new_nonorbs = define_orb_random(
            species_incl_idx,
            self.true_odds.total_spider_count,
            self.true_odds.orb_count,
        )

        # Recalculate the odds ratios and log odds ratios for the new orb/non-orb arrays
        new_log_odds_ratio = calculate_odds(
            new_orbs,
            new_nonorbs,
            self.true_odds.test_bool_mat,
            self.true_odds.orb_count,
            self.true_odds.nonorb_count,
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
        new_std_dev = np.std(new_log_odds_ratio_fltrd)

        if self.stat == "mean":
            true_stat = self.true_mean
            new_stat = new_mean
        elif self.stat == "stddev":
            true_stat = self.true_stddev
            new_stat = new_std_dev
        else:
            true_stat = self.true_skew
            new_stat = new_skew

        if self.alternative == "greater":
            if new_stat > true_stat:
                counter += 1
            else:
                pass

        else:
            if new_stat < true_stat:
                counter += 1
            else:
                pass

        self.alphas[i] = new_mean + self.z_crit * new_std_dev
        self.means[i] = new_mean
        self.stddevs[i] = new_std_dev
        self.skews[i] = new_skew

        return counter

    def _run_bootstrap(self):
        """Creating 10,000 test log odds ratio distributions with the set
        of spiders defined as "orb-weavers" chosen randomly (but still the
        same number of orb weavers as the true number)"""

        print("\n\nLAUNCHING BOOTSTRAPPING TEST\n")

        # Set the random seed for reproducibility
        random.seed(144)

        if self.maximum != self.true_odds.total_spider_count:
            print(f"** Maximum occupancy set to {self.maximum} **\n")
        else:
            pass

        if self.occupancy_threshold > 0:
            print(f"** Minimum occupancy set to {self.occupancy_threshold} **\n")

        print(f"Skew of LOR distribution: {str(self.true_skew)}")
        print(f"Standard deviation of LOR distribution: {str(self.true_stddev)}")
        print(f"Mean of LOR distribution: {str(self.true_mean)}\n")

        # Initialize lists to store the results of the bootstrapping
        self.alphas = np.zeros(10000)
        self.means = np.zeros(10000)
        self.stddevs = np.zeros(10000)
        self.skews = np.zeros(10000)

        # Initialize counters for how often the bootstrapped
        # distributions exceed the true values
        counter = 0

        # Get a list of the indices of the spiders being compared,
        # in case not all spiders are included in the analysis
        species_incl_arr = (
            self.true_odds.orb_bool_arr + self.true_odds.nonorb_bool_arr
        ).astype(bool)
        species_incl_idx = species_incl_arr.nonzero()[0]

        for i in tqdm(range(self.bootstrap_reps)):
            # Perform a single iteration of the bootstrapping test
            counter = self._bootstrap_iter(i, counter, species_incl_idx)

        p = counter / 10000
        print("Bootstrapping counter: " + str(counter))
        print(f"Probability that null is true: {p:.2f}\n")

        self.p_value = p

        # Getting average stats from across the 10000 bootstrapped distributions
        self.alpha_av = np.mean(self.alphas)
        self.mean_av = np.mean(self.means)
        self.stddev_av = np.mean(self.stddevs)
        self.skew_av = np.mean(self.skews)

        print("Average bootstrapping mean: " + str(self.mean_av))
        print("Average bootstrapping stddev: " + str(self.stddev_av))
        print("Average bootstrapping alpha=0.05: " + str(self.alpha_av) + "\n")

    def plot_bootstrap_stats(self, hist_txt):
        """Plotting the bootstrapping means, standard deviations
        and alpha thresholds to ensure the results are relatively
        tightly distributed"""

        fig, axs = plt.subplots(1, 3, figsize=(15, 5))

        fig.suptitle(
            rf"$\bf{{Bootstrapped\ distribution\ stats\ for\ {hist_txt}}}$" + "\n"
            f"Maximum occupancy = {self.maximum}, "
            f"minimum occupancy = {self.occupancy_threshold}",
            fontsize=16,
        )

        sns.histplot(
            data=self.alphas,
            kde=True,
            bins=50,
            stat="density",
            ax=axs[0],
            legend=False,
        )
        axs[0].set_title(f"Alpha={self.a} threshold")
        axs[0].set(xlabel="alpha", ylabel="Density")
        axs[0].axvline(x=self.alpha_av, linestyle="dotted")

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
        plt.show()

        return fig, axs

    def plot_bootstrap_results(
        self,
        plot_name,
        gaussfit_color="blue",
        avbootstrap_color="red",
        histcolor="lightblue",
    ):
        """Function to plot the results of the bootstrapping test"""

        fig, ax = plt.subplots(figsize=(8, 6))

        fig.suptitle(
            rf"$\bf{{Log\ odds\ ratio\ of\ {plot_name}}}$" + "\n"
            f"Maximum occupancy = {self.maximum}, "
            f"minimum occupancy = {self.occupancy_threshold}",
            fontsize=16,
        )

        sns.histplot(
            data=self.true_fltrd_log_odds_ratios,
            kde=False,
            bins=50,
            stat="density",
            ax=ax,
            legend=False,
            color=histcolor,
        )

        x = np.linspace(
            self.true_fltrd_log_odds_ratios.min(),
            self.true_fltrd_log_odds_ratios.max(),
            100,
        )

        ax.plot(
            x,
            norm.pdf(
                x,
                np.mean(self.true_fltrd_log_odds_ratios),
                np.std(self.true_fltrd_log_odds_ratios),
            ),
            color=gaussfit_color,
            linestyle="--",
            label="Gaussian fit to true\ndistribution",
        )

        ax.plot(
            x,
            norm.pdf(x, self.mean_av, self.stddev_av),
            color=avbootstrap_color,
            linestyle="--",
            label="Average bootstrapped\ndistribution",
        )

        ax.axvline(
            x=self.alpha_av,
            label="Mean bootstrapped\nthreshold for alpha=0.05,\nleft-tailed",
            linestyle="dotted",
        )

        ax.text(
            0.03,
            0.95,
            f"Bootstrapped skew = {self.skew_av:.2f}\n"
            f"True skew = {self.true_skew:.2f}\n"
            f"Bootstrapped p-value = {self.p_value:.2f}",
            transform=ax.transAxes,
            fontsize=10,
            ha="left",
            va="top",
        )

        ax.set(xlabel="Log odds ratio", ylabel="Density")

        plt.legend(bbox_to_anchor=(1.1, 1.05), fontsize=10)
        plt.tight_layout()
        plt.show()

        return fig, ax

    def get_hits_df(self, includes_udiv=True):
        """Filter the dataframe for HOGs that have LORs exceeding the
        alpha values determined by bootstrapping, and filter for those
        which include Uloborus diversus if specified."""

        # Filter the DataFrame for significant hits
        df_fltrd, total_hits = self.true_odds.filter_for_bootstrap_hits(
            self.occupancy_threshold, self.maximum, self.alpha_av
        )

        # Filter out the hits that do not include u.diversus
        df_fltrd_ud, udiv_hits = filter_for_udiv(df_fltrd, self.true_odds.genecount_df)

        if includes_udiv:
            return df_fltrd_ud, total_hits, udiv_hits

        else:
            return df_fltrd, total_hits, udiv_hits

    def get_hits_hogs(self, includes_udiv=True):
        """Get list of HOGs which are significant according to
        the bootstrapping test"""

        df_fltrd, _, _ = self.get_hits_df(includes_udiv=includes_udiv)

        hogs_list = df_fltrd.index.tolist

        return hogs_list

    def get_hits_annots(
        self,
        results_dir=None,
        includes_udiv=True,
    ):
        """Find and/or export significant results with annotations"""

        minimum = self.occupancy_threshold
        maximum = self.maximum
        test = self.true_odds.test

        # Filter the DataFrame for significant hits
        df_fltrd, total_hits, udiv_hits = self.get_hits_df(includes_udiv=includes_udiv)

        # Get annotations and save, or just return filtered dataframe and counts of hits
        if results_dir is not None:

            filename = f"{results_dir}/{test}_sig_hog_list_\
            {minimum}-{maximum}_occ_udiv_annots.csv"

            print(
                f"************** Exporting significant hits and annotations for {test}, "
                f"occ {minimum} to {maximum} **************\n"
            )

            merged_df = id_converter.main(df_fltrd, self.true_odds.hog_node_genes_tsv)
            merged_df.to_csv(filename, index=False)
            print(f"Merged results file with annotations saved to {filename}.\n\n")
            print(f"Total significant hits: {total_hits}\n")

            if includes_udiv:
                print(f", including U. diversus: {udiv_hits}\n")

            return merged_df, total_hits, udiv_hits

        else:
            print(
                "No annotations saved, just returning filtered "
                "DataFrame and counts of hits.\n"
            )
            print(
                f"Total significant hits: {total_hits}, "
                f"including U. diversus: {udiv_hits}\n"
            )
            return df_fltrd, total_hits, udiv_hits

    def print_bootstrap_results(self, total_hits, udiv_hits):
        """Function to print the results of the bootstrapping test"""

        if self.alternative == "greater":
            self.alternative = "right-tailed"
        else:
            self.alternative = "left-tailed"

        if self.maximum == self.true_odds.total_spider_count:
            maximum = "no"
        else:
            maximum = self.maximum

        print(
            "*********************** RESULTS ***********************\n\n"
            f"Bootstrapping test with {self.bootstrap_reps} repetitions "
            f"for {self.true_odds.test} ({self.alternative})\n"
            f"with minimum occupancy *{self.occupancy_threshold}* and "
            f"maximum occupancy *{maximum}* \n"
        )

        print(f"Orb-weaver list: {self.true_odds.orb_list_filename}")
        if self.true_odds.nonorb_list_filename is not None:
            print(f"Non-orb-weaver list: {self.true_odds.nonorb_list_filename}")
        print(f"Gene count file: {self.true_odds.genecount_csv}")
        print(f"Hierarchical orthogroup file: {self.true_odds.hog_node_genes_tsv}\n")

        print(
            f"Total number of HOGs in node: {len(self.true_odds.hog_list)}\n"
            f"Total number of HOGs within occupancy threshold: "
            f"{len(self.true_fltrd_log_odds_ratios)}\n"
            f"Total spiders: {self.true_odds.total_spider_count}\n"
            f"Orb-weavers: {self.true_odds.orb_count}\n"
            f"Non-orb-weavers: {self.true_odds.nonorb_count}\n"
            f"True skew: {self.true_skew}\n"
            f"True mean: {self.true_mean}\n"
            f"True standard deviation: {self.true_stddev}\n"
            f"Bootstrapping p-value for {self.stat}: {self.p_value}\n\n"
        )

        print(
            f"Bootstrapped average skew: {self.skew_av}\n"
            f"Bootstrapped average mean: {self.mean_av}\n"
            f"Bootstrapped average standard deviation: {self.stddev_av}\n"
            f"Significance threshold: {self.a}\n"
            f"Bootstrap-derived alpha threshold: {self.alpha_av}\n\n"
            "Total HOGs with significantly different LORs between "
            f"orb and non-orb (two-tailed): {total_hits}\n\n"
            "Total HOGs with significantly different LORs between "
            f"orb and non-orb (two-tailed, U. diversus present): {udiv_hits}\n"
        )


def get_locs_list(df):
    """Get list of LOCs which are significant according to the
    bootstrapping test"""

    locs_list = df["LOC"].dropna().unique()

    return locs_list


def odds_ratio_test(
    genecount_csv,
    orb_list_filename,
    hog_node_genes_tsv,
    test,
    results_dir=None,
    occupancy_threshold=0,
    max_occ=None,
    bootstrap_stat="mean",  # or "stddev" or "skew"
    alternative="less",
    alpha=0.05,
    bootstrap_reps=10000,
    nonorb_list_filename=None,
    includes_udiv=True,
):
    """Run the full odds ratio test"""

    # Get the true odds ratio results
    true_odds = OddsRatioResults(
        genecount_csv=genecount_csv,
        hog_node_genes_tsv=hog_node_genes_tsv,
        test=test,
        orb_list_filename=orb_list_filename,
        nonorb_list_filename=nonorb_list_filename,
    )

    # Run bootstrapping test
    bootstrap_test_results = BootstrapTestResults(
        true_odds=true_odds,
        stat=bootstrap_stat,
        occupancy_threshold=occupancy_threshold,
        maximum=max_occ,
        alternative=alternative,
        a=alpha,
        bootstrap_reps=bootstrap_reps,
    )

    # Get the hits and annotations
    _, hits, uhits = bootstrap_test_results.get_hits_annots()

    # Print the results of the bootstrapping test
    bootstrap_test_results.print_bootstrap_results(hits, uhits)

    # Save the results to a CSV file if specified
    if results_dir is not None:
        true_odds.results_df.to_csv(
            f"{results_dir}/{test}_odds_results_all.csv",
            index=True,
        )

        fig, _ = true_odds.plot_log_odds_ratios(
            occupancy_threshold=occupancy_threshold, maximum=max_occ
        )

        fig.savefig(
            f"{results_dir}/LOR_dists_{test}_occ_min_"
            f"{occupancy_threshold}_max_{max_occ}.png",
            dpi=300,
        )

        bootstrap_test_results.get_hits_annots(
            results_dir=results_dir, includes_udiv=includes_udiv
        )

    return bootstrap_test_results
