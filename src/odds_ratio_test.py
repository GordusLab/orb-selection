"""
Odds ratio analysis of gene loss and gain: testing whether
the odds that genes in an orthogroup are missing, in single copy,
or in multiple copies differs significantly between orb-weavers
and non-orb-weavers"""

import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, skew
from tqdm.auto import tqdm
import seaborn as sns
import src.id_converter as id_converter
from src.orthogroup_filter import drop_empty_cols


def get_genecount_arrays(genecount_csv, test):
    """Function to table of counts of genes per species in each orthogroup
    at the specified node of my family tree, and matrices of 1s and 0s
    corresponding to whether a gene is missing, duplicated, or single copy
    for each species."""

    # Counts of genes per species in each HOG
    genecount_df = pd.read_csv(genecount_csv, index_col="HOG", sep="\t")
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
    if test == "loss":
        test_bool_mat = np.asarray(genecount_mat == 0).astype(float)

    # Matrix where 1=gene with more than 1 copy, 0=gene with 1 or no copies
    elif test == "dup":
        test_bool_mat = np.asarray(genecount_mat > 1).astype(float)

    # Matrix where 1=gene in single copy, 0=gene missing or duplicated
    else:
        test_bool_mat = np.asarray(genecount_mat == 1).astype(float)

    args = {
        "genecount_df": genecount_df,
        "hog_list": hog_list,
        "all_species_arr": all_species_arr,
        "total_occ_arr": total_occ_arr,
        "test_bool_mat": test_bool_mat,
    }
    return args


def define_orb(args, orb_arr, nonorb_arr=None):
    """Function to make a numpy vector corresponding to whether
    each species is an orb-weaver or not: 1 = orb, 0 = non-orb.
    If the non-orbweavers are not specified, any animal that is
    not specified to be OW is considered non-OW."""

    # assign True to index of orbweavers in the spiders array
    orb_bool_arr = np.isin(args["all_species_arr"], orb_arr)

    total_spider_count = len(orb_bool_arr)
    orb_count = orb_bool_arr.sum()

    # If no background is specified, all spiders not identified as
    # orb-weavers are considered non-orb-weavers
    if nonorb_arr is None:
        nonorb_bool_arr = np.isin(args["all_species_arr"], orb_arr, invert=True)

    # If the non-orb-weavers are specified
    else:
        nonorb_bool_arr = np.isin(args["all_species_arr"], nonorb_arr)

    nonorb_count = nonorb_bool_arr.sum()

    nonorb_bool_arr = nonorb_bool_arr.reshape(nonorb_bool_arr.size, 1).astype(float)
    orb_bool_arr = orb_bool_arr.reshape(orb_bool_arr.size, 1).astype(float)

    print(
        f"{total_spider_count} spiders total, {orb_count} orb-weavers,\
           {nonorb_count} non-orb-weavers"
    )

    args["total_spider_count"] = total_spider_count
    args["orb_count"] = orb_count
    args["nonorb_count"] = nonorb_count
    args["orb_bool_arr"] = orb_bool_arr
    args["nonorb_bool_arr"] = nonorb_bool_arr

    return args


def get_args(genecount_csv, test, orb_arr, nonorb_arr=None):
    """Function to generate all arguments and arrays needed for LOR test"""

    genecount_args = get_genecount_arrays(genecount_csv, test)
    all_args = define_orb(genecount_args, orb_arr, nonorb_arr)
    return all_args


def calculate_odds(args):
    """Function to calculate the odds ratio and log odds ratio"""

    # I don't know why this is necessary but my kernel crashes without it
    args["orb_bool_arr"] = args["orb_bool_arr"].reshape(args["orb_bool_arr"].size, 1)
    args["nonorb_bool_arr"] = args["nonorb_bool_arr"].reshape(
        args["nonorb_bool_arr"].size, 1
    )

    # Calculate the number of orb-weavers and non-orb-weavers that are
    # [present, missing, duplicated, single copy]

    ## orb-weavers yes
    orb_yes_arr = np.matmul(args["test_bool_mat"], args["orb_bool_arr"])

    ## non-orbweavers yes
    nonorb_yes_arr = np.matmul(args["test_bool_mat"], args["nonorb_bool_arr"])

    ## orb-weavers no
    orb_no_arr = args["orb_count"] - orb_yes_arr

    ## non-orb-weavers no
    nonorb_no_arr = args["nonorb_count"] - nonorb_yes_arr

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

    odds_arrs = {}

    odds_arrs["odds_orb_arr"] = odds_orb_arr
    odds_arrs["odds_nonorb_arr"] = odds_orb_arr
    odds_arrs["log_odds_ratio_arr"] = log_odds_ratio_arr

    return odds_arrs


def occupancy_filter(arr, minimum, maximum, total_occ_arr):
    """Function to filter an array of values according to whether the
    orthogroup meets a certain occupancy threshold."""

    idx = np.asarray((total_occ_arr >= minimum) & (total_occ_arr <= maximum)).nonzero()[
        0
    ]
    fltrd_arr = arr[idx]

    return fltrd_arr


def plot_bootstrap_stats(
    num_thresholds,
    alpha_thresh,
    alphas,
    means,
    stddevs,
    alpha_av,
    mean_av,
    stddev_av,
    k,
    threshold,
    hist_txt,
    maximum,
    results_dir,
):
    """Plotting the bootstrapping means, standard deviations
    and alpha thresholds to ensure the results are relatively
    tightly distributed"""

    fig, axs = plt.subplots(
        num_thresholds, 3, figsize=(9.5, 2.5 * num_thresholds), sharex="col"
    )
    fig.suptitle(
        f"Bootstrapped distribution stats for {hist_txt} (maximum occupancy = {maximum})"
    )
    fig.subplots_adjust(top=0.90)

    sns.histplot(
        data=alphas,
        kde=True,
        bins=50,
        stat="density",
        ax=axs[k][0],
        legend=False,
    )
    axs[k][0].set_title(f"Alpha={alpha_thresh} threshold, occ ≥ {threshold}")
    axs[k][0].set(xlabel="alpha", ylabel="Density")
    axs[k][0].axvline(x=alpha_av, linestyle="dotted")

    sns.histplot(
        data=means,
        kde=True,
        bins=50,
        stat="density",
        ax=axs[k][1],
        legend=False,
    )
    axs[k][1].set_title(f"Means, occ ≥ {threshold}")
    axs[k][1].set(xlabel="mean", ylabel="Density")
    axs[k][1].axvline(x=mean_av, linestyle="dotted")

    sns.histplot(
        data=stddevs,
        kde=True,
        bins=50,
        stat="density",
        ax=axs[k][2],
        legend=False,
    )

    axs[k][2].set_title(f"Std devs, occ ≥ {threshold}")
    axs[k][2].set(xlabel="stddev", ylabel="Density")
    axs[k][2].axvline(x=stddev_av, linestyle="dotted")

    plt.tight_layout()
    plt.savefig(f"{results_dir}/bs_stats_{hist_txt}.png", dpi=300)


def bootstrap(
    thresholds,
    results_dir,
    args,
    odds,
    hist_txt=None,
    maximum=None,
    right_tailed=False,
    bootstrap_reps=10000,
):
    """Creating 10,000 test log odds ratio distributions with the set
    of spiders defined as "orb-weavers" chosen randomly (but still the
    same number of orb weavers as the true number)"""

    # Set the random seed for reproducibility
    random.seed(42)

    # Used to determine the p-value threshold for each bootstrap distribution
    if right_tailed:
        a = 0.95
    else:
        a = 0.05

    z_crit = norm.ppf(a)

    print("Launching bootstrapping test\n")

    if max is not None:
        print(f"** Maximum occupancy set to {maximum} **\n")
    else:
        maximum = args["total_spider_count"]

    bootstrap_means = []
    bootstrap_stddevs = []
    bootstrap_alphas = []

    pvals_skew = []
    pvals_mean = []
    pvals_stddev = []

    fltrd_log_odds_ratios = []

    # Get a list of the indices of the spiders being compared,
    # in case not all spiders are included in the analysis
    species_incl_arr = (args["orb_bool_arr"] + args["nonorb_bool_arr"]).astype(bool)
    species_incl_idx = species_incl_arr.nonzero()[0]

    for k, threshold in enumerate(thresholds):

        print(f"OCCUPANCY >= {threshold} RESULTS... \n")

        fltrd_log_odds_ratios = occupancy_filter(
            odds["log_odds_ratio_arr"], threshold, max, args["total_occ_arr"]
        )
        skewness = skew(fltrd_log_odds_ratios)
        mean = np.mean(fltrd_log_odds_ratios)
        stddev = np.std(fltrd_log_odds_ratios)

        print(f"Skewness of LOR distribution: {str(skewness[0])}")
        print(f"Standard deviation of LOR distribution: {str(stddev)}")
        print(f"Mean of LOR distribution: {str(mean)}\n")

        # Initialize counters for how often the bootstrapped distributions exceed the true values
        counter_skew = 0
        counter_mean = 0
        counter_stddev = 0

        alphas = np.zeros(10000)
        means = np.zeros(10000)
        stddevs = np.zeros(10000)

        for i in tqdm(range(bootstrap_reps)):

            # set all spiders to "non-orb" (0)
            orbs = np.zeros(args["total_spider_count"])

            # Generate a list of orbweaver indices chosen from the length
            # of the total list of spiders **included in this comparison**
            # these will be the new "orbweavers"
            orbs_idx = np.random.choice(
                species_incl_idx, args["orb_count"], replace=False
            )

            # set the randomized orbweaver indices to 1
            orbs[orbs_idx] = 1

            # the rest of the "included species" are non-orbs
            nonorbs = np.zeros(args["total_spider_count"])
            nonorbs_idx = np.setdiff1d(species_incl_idx, orbs_idx, assume_unique=True)

            # set the randomized orbweaver indices to 1
            nonorbs[nonorbs_idx] = 1

            bootstrap_args = {
                "orb_count": args["orb_count"],
                "nonorb_count": args["nonorb_count"],
                "orb_bool_arr": orbs,
                "nonorb_bool_arr": nonorbs,
                "test_bool_mat": args["test_bool_mat"],
            }

            new_log_odds_ratio = calculate_odds(bootstrap_args)["log_odds_ratio_arr"]
            new_log_odds_ratio_fltrd = occupancy_filter(
                new_log_odds_ratio, threshold, max, args["total_occ_arr"]
            )

            new_skewness = skew(new_log_odds_ratio_fltrd)
            new_mean = np.mean(new_log_odds_ratio_fltrd)
            new_std_dev = np.std(new_log_odds_ratio_fltrd)

            # NEED TO CHECK ON THESE CONDITIONS
            if right_tailed:  # action store true
                if new_skewness > skewness:
                    counter_skew += 1
                else:
                    pass

                if new_mean > mean:
                    counter_mean += 1
                else:
                    pass

            else:
                if new_skewness < skewness:
                    counter_skew += 1
                else:
                    pass

                if new_mean < mean:
                    counter_mean += 1
                else:
                    pass

            # same hypothesis for stddev either way
            if new_std_dev > stddev:
                counter_stddev += 1
            else:
                pass

            alphas[i] = new_mean + z_crit * new_std_dev
            means[i] = new_mean
            stddevs[i] = new_std_dev

        prob_skew = "{:.10f}".format(counter_skew / 10000)
        print("Bootstrapping counter for SKEW: " + str(counter_skew))
        print("Probability that null is true for SKEW: " + prob_skew + "\n")

        prob_mean = "{:.10f}".format(counter_mean / 10000)
        print("Bootstrapping counter for MEAN: " + str(counter_mean))
        print("Probability that null is true for MEAN: " + prob_mean + "\n")

        prob_stddev = "{:.10f}".format(counter_stddev / 10000)
        print("Bootstrapping counter for STANDARD DEVIATION: " + str(counter_stddev))
        print(
            "Probability that null is true for STANDARD DEVIATION: "
            + prob_stddev
            + "\n"
        )

        # Save the p-values
        pvals_mean.append(prob_mean)
        pvals_skew.append(prob_skew)
        pvals_stddev.append(prob_stddev)

        # Getting average stats from across the 10000 bootstrapped distributions
        alpha_av = np.mean(alphas)
        mean_av = np.mean(means)
        stddev_av = np.mean(stddevs)

        print("Average bootstrapping mean: " + str(mean_av))
        print("Average bootstrapping stddev: " + str(stddev_av))
        print("Average bootstrapping alpha=0.05: " + str(alpha_av) + "\n")

        if hist_txt is not None:
            # Plot the bootstrapped distributions and their stats
            print("Plotting bootstrapped distributions and their stats... \n")

            plot_bootstrap_stats(
                len(thresholds),
                a,
                alphas,
                means,
                stddevs,
                alpha_av,
                mean_av,
                stddev_av,
                k,
                threshold,
                hist_txt,
                maximum,
                results_dir,
            )

        else:
            pass

        # Save the average stats across all the bootstrapped distribution to use as significance
        # thresholds for the true distribution
        bootstrap_means.append(mean_av)
        bootstrap_stddevs.append(stddev_av)
        bootstrap_alphas.append(alpha_av)

        fltrd_log_odds_ratios.append(fltrd_log_odds_ratios)

    results = {
        "bootstrap_means": bootstrap_means,
        "bootstrap_stddevs": bootstrap_stddevs,
        "bootstrap_alphas": bootstrap_alphas,
        "pvals_mean": pvals_mean,
        "pvals_skew": pvals_skew,
        "pvals_stddev": pvals_stddev,
        "thresholds": thresholds,
        "max": maximum,
        "fltrd_log_odds_ratios": fltrd_log_odds_ratios,
    }

    return results


def display_bootstrap_res(results, plot_name, results_dir):
    """Function to display the results of the bootstrapping test"""

    colors = ["blue", "red"]
    colors2 = ["lightblue", "lightcoral"]

    a = len(results["thresholds"])
    fig, axs = plt.subplots(1, a)
    fig.set_figwidth(a * 4)
    fig.set_figheight(4)
    fig.suptitle(f"Log odds ratio of {plot_name}; max occupancy = {results['max']}")

    for i, threshold in enumerate(results["thresholds"]):
        bs_av = results["bootstrap_means"][i]
        bs_sd = results["bootstrap_stddevs"][i]
        bs_al = results["bootstrap_alphas"][i]

        fltrd_log_odds_ratios = results["fltrd_log_odds_ratioss"][i]

        sns.histplot(
            data=fltrd_log_odds_ratios,
            kde=False,
            bins=50,
            stat="density",
            ax=axs[i],
            legend=False,
            color=colors2[0],
        )
        # skewness = skew(fltrd_log_odds_ratios)[0]
        x = np.linspace(fltrd_log_odds_ratios.min(), fltrd_log_odds_ratios.max(), 100)

        axs[i].plot(
            x,
            norm.pdf(x, np.mean(fltrd_log_odds_ratios), np.std(fltrd_log_odds_ratios)),
            color=colors[0],
            linestyle="--",
            label="Gaussian fit",
        )

        axs[i].plot(
            x,
            norm.pdf(x, bs_av, bs_sd),
            color=colors[1],
            linestyle="--",
            label="Average bootstrapped \ndistribution",
        )

        axs[i].axvline(
            x=bs_al, label="Threshold for alpha=0.05, \nleft-tailed", linestyle="dotted"
        )
        # axs[i].text(xmin,ymax,f"Skewness = {skewness}")

        axs[i].set_title(f"Occupancy >= {threshold}")
        axs[i].set(xlabel="Log odds ratio", ylabel="Density")

    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.tight_layout()
    plt.savefig(f"{results_dir}/bs_results_{plot_name}.png", dpi=300)
    plt.show()


def results_to_df(
    args,
    test,
    results_dir,
    save_csv=False,
    threshold=None,
):
    """Function to convert the results of the odds ratio test into a pandas DataFrame
    and save it to a CSV file if specified."""

    results_df = pd.DataFrame(
        {
            "HOG": args["hog_list"],
            "Occupancy": args["total_occ_arr"].flatten(),
            f"Odds {test}, orb": args["odds_orb_arr"].flatten(),
            f"Odds {test}, non-orb": args["odds_nonorb_arr"].flatten(),
            "Log odds ratio": args["log_odds_ratio_arr"].flatten(),
        }
    )
    results_df = results_df.set_index("HOG")

    if threshold:
        results_df = results_df[
            (results_df["Occupancy"] >= threshold[0])
            & (results_df["Occupancy"] <= threshold[1])
        ]

        if save_csv:
            results_df.to_csv(
                f"{results_dir}/odds_ratio_results_{test}_{threshold[0]}-{threshold[1]}.csv",
                index=True,
            )
        else:
            pass

    else:
        if save_csv:
            results_df.to_csv(
                f"{results_dir}/odds_ratio_results_{test}_all.csv", index=True
            )
        else:
            pass

    return results_df


def get_hits_annots(
    args,
    results,
    test,
    hog_node_genes_tsv,
    results_dir,
    udiv=False,
    merge_and_save=False,
):
    """Export significant results with annotations"""

    df = results_to_df(
        args,
        test,
        results_dir,
        save_csv=merge_and_save,
    )

    dfs = []
    total_hits_list = []
    udiv_hits_list = []

    for i, minimum in enumerate(tqdm(results["thresholds"], leave=False)):

        alpha = results["bootstrap_alphas"][i]

        filename = f"{results_dir}/{test}_sig_hog_list_{minimum}-{results["max"]}_udiv_annots.csv"

        print(
            f"************** Exporting significant hits and annotations for {test}, \
                occ {minimum} to {results["max"]} **************\n"
        )

        if test == "sc":
            df_fltrd = df[
                (df["Occupancy"] >= minimum)
                & (df["Occupancy"] <= results["max"])
                & (df["Log odds ratio"] > alpha)
            ]
        else:
            df_fltrd = df[
                (df["Occupancy"] >= minimum)
                & (df["Occupancy"] <= results["max"])
                & (df["Log odds ratio"] < alpha)
            ]

        total_hits = len(df_fltrd)
        print("Total no. hits: " + str(total_hits) + "\n")
        total_hits_list.append(total_hits)

        # Filter out the hits that include u.diversus
        if udiv:
            # the udiv present function needs unfiltered LOR df
            udiv_present = args["genecount_df"].index[
                args["genecount_df"]["Uloborus_diversus"] != 0
            ]
            df_ud = df.loc[udiv_present]

            if test == "sc":
                df_ud_fltrd = df_ud[
                    (df_ud["Occupancy"] >= minimum)
                    & (df_ud["Occupancy"] <= results["max"])
                    & (df_ud["Log odds ratio"] > alpha)
                ]
            else:
                df_ud_fltrd = df_ud[
                    (df_ud["Occupancy"] >= minimum)
                    & (df_ud["Occupancy"] <= results["max"])
                    & (df_ud["Log odds ratio"] < alpha)
                ]

            udiv_hits = len(df_ud_fltrd)
            print("No. hits after filtering for U. div: " + str(udiv_hits) + "\n")
            udiv_hits_list.append(udiv_hits)

            # Get annotations and save, or just return filtered dataframe and counts of hits
            if merge_and_save:
                df_ud_fltrd.to_csv(filename)
                merged_df = id_converter.main(filename, hog_node_genes_tsv)
                merged_df.to_csv(filename, index=False)
                print(f"Merged results file with annotations saved to {filename}.\n\n")

                dfs.append(merged_df)
            else:
                dfs.append(df_ud_fltrd)

        # or don't
        else:
            if merge_and_save:
                df_fltrd.to_csv(filename)
                merged_df = id_converter.main(filename, hog_node_genes_tsv)
                merged_df.to_csv(filename, index=False)
                print(f"Merged results file with annotations saved to {filename}.\n\n")

                dfs.append(merged_df)

            else:
                dfs.append(df_fltrd)

    if udiv:
        return dfs, total_hits_list, udiv_hits_list
    else:
        return dfs, total_hits_list


def save_test_stats(
    args,
    results,
    test,
    orb_list_file,
    genecount_csv,
    hog_node_genes_tsv,
    results_dir,
    right_tailed,
):
    """Save the test stats to a tsv file"""

    num_thr = len(results["thresholds"])
    means = [np.mean(results["fltrd_log_odds_ratioss"][i]) for i in range(num_thr)]
    skews = [skew(results["fltrd_log_odds_ratioss"][i])[0] for i in range(num_thr)]
    stddevs = [np.std(results["fltrd_log_odds_ratioss"][i]) for i in range(num_thr)]

    hits = get_hits_annots(
        args,
        results,
        test,
        hog_node_genes_tsv,
        results_dir,
        udiv=True,
        merge_and_save=False,
    )

    d = {
        "thresholds": results["thresholds"],
        "mean of true LOR dist": means,
        "std dev of true LOR dist": stddevs,
        "skew of true LOR dist": skews,
        "avg mean of BSd LOR dists": results["bootstrap_means"],
        "avg std dev of BSd LOR dists": results["bootstrap_stddevs"],
        "BSd alpha threshold": results["bootstrap_alphas"],
        "p-value for mean": results["pvals_mean"],
        "p-value for std dev": results["pvals_stddev"],
        "p-value for skew": results["pvals_skew"],
        "total no. significant hog_list": hits[1],
        "No. significant hog_list inc. U. div": hits[2],
    }

    stats_df = pd.DataFrame(data=d)

    filename = f"{results_dir}/{test}_bs_results_stats.tsv"

    stats_df.to_csv(filename, sep="\t", index=False)

    if right_tailed:
        tail = "right-tailed"
    else:
        tail = "left-tailed"

    comments = (
        f"# Orbweaver list used: {orb_list_file}\n"
        + f"# Genecount file used: {genecount_csv}\n"
        + f"# This analysis includes {str(args['total_spider_count'])} total spiders, \
        with {str(args['orb_count'])} designated as orbweavers, and \
        {str(args['nonorb_count'])} designated as non-orbweavers.\n"
        + f"# Maximum occupancy set to {str(results['max'])}\n"
        + f"# Bootsrapping of log odds ratio distribution of {test}, \
        orbweavers to non-orbweavers, {tail}"
    )

    with open(filename, "r+", encoding="utf-8") as f:
        content = f.read()
        f.seek(0, 0)
        f.write(comments.rstrip("\r\n") + "\n" + content)
        f.close()


def plot_log_odds_ratios(
    thresholds, log_odds_ratio_arr, test, maximum, total_occ_arr, results_dir
):
    """Function to plot the log odds ratio results for a certain test,
    given occupancy thresholds (no bootstrapping)"""

    a = len(thresholds)
    fig, axs = plt.subplots(1, len(thresholds))
    fig.set_figwidth(3.5 * a)
    fig.set_figheight(4)
    fig.suptitle(f"Log odds ratio of {test}, orb:non-orb")

    for i, threshold in enumerate(thresholds):
        fltrd_log_odds_ratios = occupancy_filter(
            log_odds_ratio_arr, threshold, maximum, total_occ_arr
        )

        # Plotting the LORs
        sns.histplot(
            data=fltrd_log_odds_ratios,
            kde=False,
            bins=100,
            stat="density",
            ax=axs[i],
            legend=False,
        )

        axs[i].set_title(f"Occupancy >= {threshold}")
        axs[i].set(xlabel="Log odds ratio", ylabel="Density")

    plt.tight_layout()
    plt.savefig(f"{results_dir}/LOR_dists_{test}.png", dpi=300)
    plt.show()


def display_and_save(
    results,
    args,
    test,
    results_dir,
    orb_list_file,
    genecount_csv,
    hog_node_genes_tsv,
    right_tailed=False,
):
    """Display results and save test stats"""

    hits = get_hits_annots(
        args,
        results,
        test,
        hog_node_genes_tsv,
        results_dir,
        udiv=True,
        merge_and_save=True,
    )

    display_bootstrap_res(results, test, results_dir)

    save_test_stats(
        hits,
        args,
        results,
        test,
        orb_list_file,
        genecount_csv,
        results_dir,
        right_tailed,
    )


def odds_ratio_test(
    genecount_csv,
    orb_list_file,
    results_dir,
    test,
    thresholds,
    max_occ=None,
    right_tailed=False,
    bootstrap_reps=10000,
    nonorb_list_file=None,
):
    """Run the full odds ratio test"""

    tests = ["sc", "dup", "loss"]
    if test not in tests:
        raise ValueError("Invalid test. Expected one of: %s" % tests)

    orb_arr = np.loadtxt(orb_list_file, dtype=str)

    if nonorb_list_file is not None:
        nonorb_arr = np.loadtxt(nonorb_list_file, dtype=str)
    else:
        nonorb_arr = None

    args = get_args(genecount_csv, test, orb_arr, nonorb_arr)

    odds = calculate_odds(args)

    res = bootstrap(
        thresholds,
        results_dir,
        args,
        odds,
        hist_txt=test,
        maximum=max_occ,
        right_tailed=right_tailed,
        bootstrap_reps=bootstrap_reps,
    )

    return args, odds, res
