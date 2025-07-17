"""
Odds ratio analysis of gene loss and gain: testing whether
the odds that genes in an orthogroup are missing, in single copy,
or in multiple copies differs significantly between orb-weavers
and non-orb-weavers"""

import copy
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

    idx = np.asarray((total_occ_arr >= minimum) & (total_occ_arr <= maximum)).nonzero()[
        0
    ]
    fltrd_arr = arr[idx]

    return fltrd_arr


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

    def get_genecount_arrays(self):
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

        if self.test is not None:
            if self.test == "loss":
                self.test_bool_mat = loss_bool_mat
            elif self.test == "dup":
                self.test_bool_mat = dup_bool_mat

    def define_orb(self):
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
            f"{total_spider_count} spiders total, {orb_count} orb-weavers,\
            {nonorb_count} non-orb-weavers"
        )

        self.total_spider_count = total_spider_count
        self.orb_count = orb_count
        self.nonorb_count = nonorb_count
        self.orb_bool_arr = orb_bool_arr
        self.nonorb_bool_arr = nonorb_bool_arr

    def calculate_odds(self):
        """Function to calculate the odds ratio and log odds ratio"""

        # I don't know why this is necessary but my kernel crashes without it
        self.orb_bool_arr = self.orb_bool_arr.reshape(self.orb_bool_arr.size, 1)
        self.nonorb_bool_arr = self.nonorb_bool_arr.reshape(
            self.nonorb_bool_arr.size, 1
        )

        # Calculate the number of orb-weavers and non-orb-weavers that are
        # [present, missing, duplicated, single copy]

        ## orb-weavers yes
        orb_yes_arr = np.matmul(self.test_bool_mat, self.orb_bool_arr)

        ## non-orbweavers yes
        nonorb_yes_arr = np.matmul(self.test_bool_mat, self.nonorb_bool_arr)

        ## orb-weavers no
        orb_no_arr = self.orb_count - orb_yes_arr

        ## non-orb-weavers no
        nonorb_no_arr = self.nonorb_count - nonorb_yes_arr

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

        self.odds_orb_arr = odds_orb_arr
        self.odds_nonorb_arr = odds_nonorb_arr
        self.log_odds_ratio_arr = log_odds_ratio_arr

    def results_to_df(
        self, occupancy_threshold=0, maximum=None, save_csv=False, results_dir=None
    ):
        """Function to convert the results of the odds ratio calculations
        into a DataFrame for easier analysis and plotting"""

        # Create a DataFrame with the results
        results_df = pd.DataFrame(
            {
                "HOG": self.hog_list,
                "occupany": self.total_occ_arr.flatten(),
                f"odds {self.test}, orb": self.odds_orb_arr.flatten(),
                f"odds {self.test}, non-orb": self.odds_nonorb_arr.flatten(),
                "log_odds_ratio": self.log_odds_ratio_arr.flatten(),
            }
        )
        results_df = results_df.set_index("HOG")

        # Filter the results based on the occupancy threshold
        if occupancy_threshold > 0:
            results_df = results_df[
                results_df["occupany"] >= occupancy_threshold
            ].copy()

        # If a maximum occupancy is specified, filter the results further
        if maximum is not None:
            results_df = results_df[
                (results_df["occupany"] >= occupancy_threshold)
                & (results_df["occupany"] <= maximum)
            ].copy()
        else:
            maximum = "none"

        # Save the results to a CSV file if specified
        if save_csv:
            results_df.to_csv(
                f"{results_dir}/{self.test}_odds_results_occ_min_{occupancy_threshold}_max_{maximum}.csv",
                index=True,
            )

        return results_df

    def plot_log_odds_ratios(
        self,
        results_dir,
        occupancy_threshold=0,
        maximum=None,
        save_fig=False,
    ):
        """Function to plot the log odds ratio results for a certain test,
        given occupancy thresholds (no bootstrapping)"""

        fig, ax = plt.subplots()
        fig.suptitle(f"Log odds ratio of {self.test}, orb:non-orb")

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

        ax.set_title(f"Occupancy >= {occupancy_threshold}, max occupancy {maximum}")
        ax.set(xlabel="Log odds ratio", ylabel="Density")

        plt.tight_layout()
        if save_fig:
            plt.savefig(
                f"{results_dir}/LOR_dists_{self.test}_occ_min_{occupancy_threshold}_max_{maximum}.png",
                dpi=300,
            )

        plt.show()

        return fig, ax


class BootstrapTest:
    """Class to perform the bootstrapping test for the odds ratio"""

    def __init__(
        self,
        odds_ratio_results,
        stat,
        occupancy_threshold=0,
        maximum=None,
        alternative="less",
        bootstrap_reps=10000,
    ):
        """Initialize the bootstrapping test for a given odds ratio results object"""

        self.odds_ratio_results = odds_ratio_results
        self.stat = stat
        self.occupancy_threshold = occupancy_threshold
        self.maximum = maximum
        self.alternative = alternative
        self.bootstrap_reps = bootstrap_reps

        # Unpack the odds ratio results
        self.hog_node_genes_tsv = odds_ratio_results.hog_node_genes_tsv
        self.genecount_csv = odds_ratio_results.genecount_csv
        self.true_orb_bool_arr = odds_ratio_results.orb_bool_arr
        self.true_nonorb_bool_arr = odds_ratio_results.nonorb_bool_arr
        self.test_bool_mat = odds_ratio_results.test_bool_mat
        self.orb_count = odds_ratio_results.orb_count
        self.nonorb_count = odds_ratio_results.nonorb_count
        self.total_spider_count = odds_ratio_results.total_spider_count
        self.total_occ_arr = odds_ratio_results.total_occ_arr
        self.true_log_odds_ratio_arr = odds_ratio_results.log_odds_ratio_arr

        # Determine the p-value threshold based on the alternative hypothesis
        if self.alternative == "greater":
            self.a = 0.05
            self.z_crit = norm.ppf(0.95)
        else:
            self.a = 0.95
            self.z_crit = norm.ppf(0.05)

        # Filter the true log odds ratios based on the occupancy threshold
        self.true_fltrd_log_odds_ratios = occupancy_filter(
            self.true_log_odds_ratio_arr,
            self.occupancy_threshold,
            self.maximum,
            self.total_occ_arr,
        )

        # Calculate the skewness, mean, and standard deviation of the true log odds ratios
        self.true_skewness = skew(self.true_fltrd_log_odds_ratios)[0]
        self.true_mean = np.mean(self.true_fltrd_log_odds_ratios)
        self.true_stddev = np.std(self.true_fltrd_log_odds_ratios)

    def define_orb_random(self):
        """Function to define the orb-weaver and non-orb-weaver arrays for bootstrapping"""

        # Get a list of the indices of the spiders being compared,
        # in case not all spiders are included in the analysis
        species_incl_arr = (self.true_orb_bool_arr + self.true_nonorb_bool_arr).astype(
            bool
        )
        species_incl_idx = species_incl_arr.nonzero()[0]

        # Set the random seed for reproducibility
        random.seed(42)

        # Create new orb and non-orb arrays based on the specified counts
        new_orbs = np.zeros(self.total_spider_count)
        new_nonorbs = np.zeros(self.total_spider_count)

        # Randomly select indices for the new orb-weavers
        new_orbs_idx = np.random.choice(species_incl_idx, self.orb_count, replace=False)
        new_orbs[new_orbs_idx] = 1

        # The rest of the species are non-orb-weavers
        new_nonorbs_idx = np.setdiff1d(
            species_incl_idx, new_orbs_idx, assume_unique=True
        )
        new_nonorbs[new_nonorbs_idx] = 1

        return new_orbs, new_nonorbs

    def bootstrap_iter(self, i, counter):
        """Function to perform a single iteration of the bootstrapping test"""

        # Create a new set of orb and non-orb arrays for bootstrapping
        new_orbs, new_nonorbs = self.define_orb_random()

        dummy_odds_ratio_results = copy.deepcopy(self.odds_ratio_results)

        # Create a new OddsRatioResults object with the new orb and non-orb arrays
        dummy_odds_ratio_results.orb_bool_arr = new_orbs
        dummy_odds_ratio_results.nonorb_bool_arr = new_nonorbs

        dummy_odds_ratio_results.calculate_odds()

        new_log_odds_ratio = dummy_odds_ratio_results.log_odds_ratio_arr
        new_log_odds_ratio_fltrd = occupancy_filter(
            new_log_odds_ratio,
            self.occupancy_threshold,
            self.maximum,
            self.total_occ_arr,
        )

        # Calculate the statistics for the new bootstrapped distribution
        new_skewness = skew(new_log_odds_ratio_fltrd)[0]
        new_mean = np.mean(new_log_odds_ratio_fltrd)
        new_std_dev = np.std(new_log_odds_ratio_fltrd)

        if self.stat == "mean":
            true_stat = self.true_mean
            new_stat = new_mean
        elif self.stat == "stddev":
            true_stat = self.true_stddev
            new_stat = new_std_dev
        else:
            true_stat = self.true_skewness
            new_stat = new_skewness

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
        self.skews[i] = new_skewness

        return counter

    def run_bootstrap(self):
        """Creating 10,000 test log odds ratio distributions with the set
        of spiders defined as "orb-weavers" chosen randomly (but still the
        same number of orb weavers as the true number)"""

        print("Launching bootstrapping test\n")

        if self.maximum is not None:
            print(f"** Maximum occupancy set to {self.maximum} **\n")
        else:
            self.maximum = self.total_spider_count

        if self.occupancy_threshold > 0:
            print(f"** Minimum occupancy set to {self.occupancy_threshold} **\n")

        print(f"Skewness of LOR distribution: {str(self.true_skewness)}")
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

        for i in tqdm(range(self.bootstrap_reps)):
            # Perform a single iteration of the bootstrapping test
            counter = self.bootstrap_iter(i, counter)

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

    def plot_bootstrap_stats(self, hist_txt, results_dir):
        """Plotting the bootstrapping means, standard deviations
        and alpha thresholds to ensure the results are relatively
        tightly distributed"""

        fig, axs = plt.subplots(1, 3, figsize=(9.5, 2.5))

        fig.suptitle(f"Bootstrapped distribution stats for {hist_txt}")
        fig.text(
            0.5,
            0.96,
            f"Maximum occupancy = {self.maximum}, \
                 minimum occupancy = {self.occupancy_threshold}",
            ha="center",
            va="top",
            fontsize=12,
        )

        fig.subplots_adjust(top=0.90)

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
        plt.savefig(f"{results_dir}/{hist_txt}_bootstrap_stats.png", dpi=300)

    def display_bootstrap_res(self, plot_name, results_dir):
        """Function to display the results of the bootstrapping test"""

        colors = ["blue", "red"]
        colors2 = ["lightblue", "lightcoral"]

        fig, ax = plt.subplots()
        fig.suptitle(f"Log odds ratio of {plot_name}")
        fig.text(
            0.5,
            0.96,
            f"Maximum occupancy = {self.maximum}, \
                 minimum occupancy = {self.occupancy_threshold}",
            ha="center",
            va="top",
            fontsize=12,
        )

        sns.histplot(
            data=self.true_fltrd_log_odds_ratios,
            kde=False,
            bins=50,
            stat="density",
            ax=ax,
            legend=False,
            color=colors2[0],
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
            color=colors[0],
            linestyle="--",
            label="Gaussian fit to true distribution",
        )

        ax.plot(
            x,
            norm.pdf(x, self.mean_av, self.stddev_av),
            color=colors[1],
            linestyle="--",
            label="Average bootstrapped \ndistribution",
        )

        ax.axvline(
            x=self.alpha_av,
            label="Mean bootstrapped threshold for alpha=0.05, \nleft-tailed",
            linestyle="dotted",
        )

        ax.text(
            f"Bootstrapped skew = {self.skew_av:.2f}\n \
                True skew = {self.true_skewness:.2f}",
            xy=(0.05, 0.95),
            xycoords="axes fraction",
            fontsize=10,
            ha="left",
            va="top",
        )
        ax.text(
            f"Bootstrapped p-value = {self.p_value:.2f}",
            xy=(0.05, 0.90),
            xycoords="axes fraction",
            fontsize=10,
            ha="left",
            va="top",
        )

        ax.set(xlabel="Log odds ratio", ylabel="Density")

        plt.legend(bbox_to_anchor=(1.1, 1.05))
        plt.tight_layout()
        plt.savefig(f"{results_dir}/{plot_name}_bootstrap_results.png", dpi=300)
        plt.show()

    def filter_for_hits(self, df, minimum, maximum, alpha):
        """Filter the DataFrame based on the occupancy and log odds ratio"""
        if self.alternative == "greater":
            df_fltrd = df[
                (df["Occupancy"] >= minimum)
                & (df["Occupancy"] <= maximum)
                & (df["Log odds ratio"] > alpha)
            ]
        else:
            df_fltrd = df[
                (df["Occupancy"] >= minimum)
                & (df["Occupancy"] <= maximum)
                & (df["Log odds ratio"] < alpha)
            ]

        self.total_hits = len(df_fltrd)
        print("Total no. hits: " + str(self.total_hits) + "\n")

        return df_fltrd, self.total_hits

    def filter_for_udiv(self, df_fltrd):
        """Filter the DataFrame for hits that include U. diversus"""

        genecount_df = self.odds_ratio_results.genecount_df
        udiv_present = genecount_df[genecount_df["Uloborus_diversus"] != 0]

        udiv_present_hogs = udiv_present.index.values
        udiv_present_hogs = set(udiv_present_hogs)
        df_ud_fltrd = df_fltrd[df_fltrd.index.isin(udiv_present_hogs)]

        self.udiv_hits = len(df_ud_fltrd)
        print("No. hits after filtering for U. div: " + str(self.udiv_hits) + "\n")

        return df_ud_fltrd, self.udiv_hits

    def get_hits_annots(
        self,
        results_dir,
        udiv=True,
        merge_and_save=False,
    ):
        """Export significant results with annotations"""

        df = self.odds_ratio_results.results_to_df(
            occupancy_threshold=self.occupancy_threshold,
            save_csv=merge_and_save,
            results_dir=results_dir,
        )

        alpha = self.alpha_av
        minimum = self.occupancy_threshold
        maximum = self.maximum
        test = self.odds_ratio_results.test

        filename = f"{results_dir}/{test}_sig_hog_list_\
            {minimum}-{maximum}_occ_udiv_annots.csv"

        print(
            f"************** Exporting significant hits and annotations for {test}, \
                occ {minimum} to {maximum} **************\n"
        )

        # Filter the DataFrame for significant hits
        df_fltrd, self.total_hits = self.filter_for_hits(df, minimum, maximum, alpha)

        if udiv:
            # Filter out the hits that include u.diversus
            df_fltrd, self.udiv_hits = self.filter_for_udiv(df_fltrd)

        # Get annotations and save, or just return filtered dataframe and counts of hits
        if merge_and_save:
            merged_df = id_converter.main(df_fltrd, self.hog_node_genes_tsv)
            merged_df.to_csv(filename, index=False)
            print(f"Merged results file with annotations saved to {filename}.\n\n")

            return merged_df, self.total_hits, self.udiv_hits

        else:
            print(
                "No annotations saved, just returning filtered DataFrame and counts of hits.\n\n"
            )
            return df_fltrd, self.total_hits, self.udiv_hits

    #### FIX
    # def save_bootstrap_stats(
    #     self,
    #     results_dir,
    #     udiv=True,
    #     merge_and_save=False
    # ):
    #     """Save the bootstrap stats to a text file"""

    #     df, self.total_hits, ud_hits = self.get_hits_annots(
    #         results_dir,
    #         udiv,
    #         merge_and_save
    #     )

    #     stats_df = pd.DataFrame(
    #         {
    #             "thresholds": self.occupancy_threshold,
    #             "mean of true LOR dist": self.true_mean,
    #             "std dev of true LOR dist": self.true_stddev,
    #             "skew of true LOR dist": self.true_skewness,
    #             "avg mean of BSd LOR dists": self.mean_av,
    #             "avg std dev of BSd LOR dists": self.stddev_av,
    #             "BSd alpha threshold": self.alpha_av,
    #             "p-value for mean": self.p_value,
    #             "p-value for std dev": self.p_value,
    #             "p-value for skew": self.p_value,
    #             "total no. significant hog_list": self.total_hits,
    #             "No. significant hog_list inc. U. div": ud_hits,
    #         }

    #     )

    #     filename = f"{results_dir}/{test}_bs_results_stats.tsv"

    #     stats_df.to_csv(filename, sep="\t", index=False)

    #     if alternative:
    #         tail = "right-tailed"
    #     else:
    #         tail = "left-tailed"

    #     comments = (
    #         f"# Orbweaver list used: {orb_list_filename}\n"
    #         + f"# Genecount file used: {genecount_csv}\n"
    #         + f"# This analysis includes {str(args['total_spider_count'])} total spiders, \
    #         with {str(args['orb_count'])} designated as orbweavers, and \
    #         {str(args['nonorb_count'])} designated as non-orbweavers.\n"
    #         + f"# Maximum occupancy set to {str(results['max'])}\n"
    #         + f"# Bootsrapping of log odds ratio distribution of {test}, \
    #         orbweavers to non-orbweavers, {tail}"
    #     )

    #     with open(filename, "r+", encoding="utf-8") as f:
    #         content = f.read()
    #         f.seek(0, 0)
    #         f.write(comments.rstrip("\r\n") + "\n" + content)
    #         f.close()

    def print_bootstrap_results(self):
        """Function to print the results of the bootstrapping test"""

        if self.alternative == "greater":
            self.alternative = "right-tailed"
        else:
            self.alternative = "left-tailed"

        if self.maximum is None:
            self.maximum = "no"

        print(
            f"Bootstrapping test with {self.bootstrap_reps} repetitions for \
              {self.odds_ratio_results.test} ({self.alternative})\n \
              with {self.occupancy_threshold} minimum occupancy and \
              {self.maximum} maximum occupancy \n\n"
        )

        print(f"Orb-weaver list: {self.odds_ratio_results.orb_list_filename}\n")
        if self.odds_ratio_results.nonorb_list_filename is not None:
            print(
                f"Non-orb-weaver list: {self.odds_ratio_results.nonorb_list_filename}\n"
            )
        print(f"Hierarchical orthogroup file: {self.hog_node_genes_tsv}\n\n")

        print(
            f"Total spiders: {self.total_spider_count}\n"
            f"Orb-weavers: {self.orb_count}, Non-orb-weavers: {self.nonorb_count}\n"
            f"True skewness: {self.true_skewness}\n"
            f"True mean: {self.true_mean}\n"
            f"True standard deviation: {self.true_stddev}\n\n"
        )

        print(
            f"Bootstrapped average skewness: {self.skew_av}\n"
            f"Bootstrapped average mean: {self.mean_av}\n"
            f"Bootstrapped average standard deviation: {self.stddev_av}\n"
            f"Bootstrapped average alpha threshold: {self.alpha_av}\n"
            f"Bootstrapping p-value: {self.p_value}\n"
            f"Total significant hits: {self.total_hits}\n"
            f"Total significant hits including U. diversus: \
                {self.udiv_hits}\n"
        )


### NEED TO GO THROUGH THIS
def odds_ratio_test(
    genecount_csv,
    orb_list_filename,
    hog_node_genes_tsv,
    results_dir,
    test,
    occupancy_threshold=0,
    max_occ=None,
    alternative="less",
    bootstrap_reps=10000,
    nonorb_list_filename=None,
):
    """Run the full odds ratio test"""

    # Initialize the OddsRatioResults object
    odds_ratio_results = OddsRatioResults(
        genecount_csv=genecount_csv,
        hog_node_genes_tsv=hog_node_genes_tsv,
        test=test,
        orb_list_filename=orb_list_filename,
        nonorb_list_filename=nonorb_list_filename,
    )

    # Get the gene count arrays and define orb-weavers
    odds_ratio_results.get_genecount_arrays()
    odds_ratio_results.define_orb()
    odds_ratio_results.calculate_odds()

    # Run bootstrapping test
    bootstrap_test = BootstrapTest(
        odds_ratio_results=odds_ratio_results,
        stat="mean",  # or "stddev" or "skew"
        occupancy_threshold=occupancy_threshold,
        maximum=max_occ,
        alternative=alternative,
        bootstrap_reps=bootstrap_reps,
    )
    bootstrap_test.run_bootstrap()
    bootstrap_test.print_bootstrap_results()
    bootstrap_test.plot_bootstrap_stats(
        hist_txt=f"{odds_ratio_results.test}_odds_ratio",
        results_dir=results_dir,
    )
    bootstrap_test.display_bootstrap_res(
        plot_name=f"{odds_ratio_results.test}_odds_ratio",
        results_dir=results_dir,
    )
    bootstrap_test.get_hits_annots(
        results_dir=results_dir,
        udiv=True,
        merge_and_save=True,
    )
    # bootstrap_test.save_bootstrap_stats(
    #     results_dir=results_dir,
    #     udiv=True,
    #     merge_and_save=True,
    # )

    return bootstrap_test
