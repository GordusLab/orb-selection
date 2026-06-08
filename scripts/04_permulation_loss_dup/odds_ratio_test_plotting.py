"""Plotting helpers for the odds ratio permulation workflow."""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import norm

plt.rcParams["font.family"] = "Verdana"


def format_stat(value, ndigits=2):
    """Format a float while normalizing negative zero to positive zero."""
    rounded = round(float(value), ndigits)
    if rounded == 0.0:
        rounded = 0.0
    return f"{rounded:.{ndigits}f}"


def format_ci(ci_values, ndigits=3):
    """Format CI vectors like [lower, upper] with fixed precision."""
    ci_arr = np.asarray(ci_values, dtype=float).flatten()
    if ci_arr.size == 0:
        return "[]"
    formatted = ", ".join(f"{x:.{ndigits}f}" for x in ci_arr)
    return f"[{formatted}]"


def save_figure(fig, output_path):
    """Save a matplotlib figure with the repo's standard settings."""
    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.3,
    )


def get_permulation_plot_data(results, test, bins=100):
    """Return shared arrays used by permulation plotting functions."""
    true_vals = getattr(results, f"true_fltrd_{test}_lors")
    x = np.linspace(true_vals.min(), true_vals.max(), 100)
    avg_pdf = norm.pdf(
        x, getattr(results, f"{test}_mean_av"), getattr(results, f"{test}_stddev_av")
    )
    true_pdf = norm.pdf(
        x,
        getattr(results, f"true_mean_{test}"),
        getattr(results, f"true_stddev_{test}"),
    )
    hist_vals, _ = np.histogram(true_vals, bins=bins, density=True)
    hist_max = hist_vals.max() if hist_vals.size else 0.0
    y_max = max(avg_pdf.max(), true_pdf.max(), hist_max) * 1.05
    if y_max == 0:
        y_max = 1.0
    return true_vals, x, avg_pdf, true_pdf, y_max


def get_permulation_thresholds(results, test):
    """Return the confidence interval vector for the chosen test."""
    if test == "loss":
        return results.loss_ci_av
    if test == "dup":
        return results.dup_ci_av
    raise ValueError(f"Invalid test type: {test}. Must be 'loss' or 'dup'.")


def plot_average_permulation_curve(
    ax,
    x,
    avg_pdf,
    color,
    linestyle="-",
    fill=True,
):
    ax.plot(
        x,
        avg_pdf,
        color=color,
        linestyle=linestyle,
        zorder=4,
        label="Average permulated\ndistribution",
    )
    if fill:
        ax.fill_between(x, avg_pdf, alpha=0.2, color=color, zorder=0)


def plot_true_histogram(ax, true_vals, bins, hist_color, hist_alpha, edgecolor):
    ax.hist(
        true_vals,
        bins=bins,
        density=True,
        color=hist_color,
        alpha=hist_alpha,
        edgecolor=edgecolor,
        label="True distribution",
        zorder=3,
    )


def plot_true_gaussian_fit(ax, x, true_pdf, color):
    ax.plot(
        x,
        true_pdf,
        color=color,
        linestyle="--",
        zorder=4,
        label="Gaussian fit to\ntrue distribution",
    )


def plot_stat_histogram(ax, values, binwidth, hist_color, hist_alpha, edgecolor):
    sns.histplot(
        data=values,
        binwidth=binwidth,
        stat="count",
        ax=ax,
        legend=False,
        color=hist_color,
        alpha=hist_alpha,
        edgecolor=edgecolor,
    )


def plot_threshold_lines(results, ax, test, thresholds_color):
    ci_av = get_permulation_thresholds(results, test)
    ax.axvline(
        x=ci_av[0],
        label=f"Mean permulated\nthresholds for\nalpha={results.alpha}",
        linestyle="dotted",
        color=thresholds_color,
        zorder=4,
    )
    ax.axvline(
        x=ci_av[1],
        linestyle="dotted",
        color=thresholds_color,
        zorder=4,
    )


def style_density_axes(
    ax,
    x,
    y_max,
    legend_fontsize=13,
    axis_label_fontsize=14,
    tick_fontsize=13,
    show_legend=True,
):
    ax.set_xlabel("Log odds ratio", fontsize=axis_label_fontsize, fontweight="bold")
    ax.set_ylabel("Density", fontsize=axis_label_fontsize, fontweight="bold")
    ax.set_ylim(bottom=0, top=y_max)
    ax.set_xlim(x.min(), x.max())
    plt.setp(ax.get_xticklabels(), fontsize=tick_fontsize)
    plt.setp(ax.get_yticklabels(), fontsize=tick_fontsize)
    if show_legend:
        ax.legend(
            fontsize=legend_fontsize,
            loc="upper right",
            ncol=1,
            labelspacing=0.8,
            handlelength=1.5,
        )


def style_stat_axes(
    ax,
    xlabel="Means",
    axis_label_fontsize=12,
    legend_fontsize=10,
    xlim=None,
    ylim=None,
    show_legend=True,
):
    ax.set(xlabel=xlabel, ylabel="Count")
    ax.xaxis.label.set_fontsize(axis_label_fontsize)
    ax.xaxis.label.set_fontweight("bold")
    ax.yaxis.label.set_fontsize(axis_label_fontsize)
    ax.yaxis.label.set_fontweight("bold")
    if show_legend:
        ax.legend(fontsize=legend_fontsize)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)


def plot_permulation_stats_vertical_panels(
    panel_specs,
    title="Permulated means",
    xlim=None,
    binwidth=0.03,
    fig_size=(9, 10),
    axis_label_fontsize=14,
    tick_fontsize=10,
    title_fontsize=15,
    panel_title_fontsize=13,
    legend_fontsize=9,
    show_legend=True,
    ylims=None,
):
    """Plot permulation means as stacked vertical panels.

    panel_specs should be an iterable of (panel_label, values, perm_mean, true_mean, color).
    """

    fig, axs = plt.subplots(len(panel_specs), 1, figsize=fig_size, sharex=True)
    axs = np.atleast_1d(axs)

    if xlim is None:
        all_values = np.concatenate([np.asarray(spec[1]) for spec in panel_specs])
        xlim = (float(np.min(all_values)), float(np.max(all_values)))

    bins = np.arange(xlim[0], xlim[1] + binwidth, binwidth)
    panel_ymax = 0.0

    for ax, (panel_label, values, perm_mean, true_mean, color) in zip(axs, panel_specs):
        counts, _, _ = ax.hist(
            values,
            bins=bins,
            histtype="stepfilled",
            linewidth=1.3,
            color=color,
            edgecolor=color,
            alpha=0.22,
        )
        panel_ymax = max(panel_ymax, float(np.max(counts)) if np.size(counts) else 0.0)
        ax.axvline(
            perm_mean,
            color="black",
            linestyle=":",
            linewidth=1.8,
            alpha=0.95,
        )
        ax.axvline(
            true_mean,
            color=color,
            linestyle="--",
            linewidth=1.8,
            alpha=0.95,
        )
        ax.set_title(panel_label, fontsize=panel_title_fontsize)
        ax.set_ylabel("Count", fontsize=axis_label_fontsize, fontweight="bold")
        ax.tick_params(axis="both", labelsize=tick_fontsize)
        ax.set_xlim(xlim)

        if show_legend:
            legend_handles = [
                Line2D([], [], color=color, linewidth=8, alpha=0.22, label="Histogram"),
                Line2D(
                    [],
                    [],
                    color="black",
                    linestyle=":",
                    linewidth=1.8,
                    label="Permulated mean",
                ),
                Line2D(
                    [],
                    [],
                    color=color,
                    linestyle="--",
                    linewidth=1.8,
                    label="True mean",
                ),
            ]
            ax.legend(handles=legend_handles, fontsize=legend_fontsize, frameon=True, loc="upper right")

    axs[-1].set_xlabel("Means", fontsize=axis_label_fontsize, fontweight="bold")

    shared_ymax = panel_ymax * 1.05 if panel_ymax > 0 else 1.0
    if ylims is not None:
        if isinstance(ylims, (tuple, list)) and len(ylims) == 2 and isinstance(ylims[0], (tuple, list)):
            shared_ymax = max(shared_ymax, max(float(y[1]) for y in ylims))
        elif isinstance(ylims, (tuple, list)) and len(ylims) == 2:
            shared_ymax = max(shared_ymax, float(ylims[1]))

    for ax in axs:
        ax.set_ylim(0, shared_ymax)

    if title:
        fig.suptitle(title, fontsize=title_fontsize)

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig, axs


def plot_permulation_hist_overlay(
    series_specs,
    title="Permulated means: overlaid histograms",
    xlim=None,
    binwidth=0.03,
    fig_size=(9, 6),
    axis_label_fontsize=14,
    tick_fontsize=10,
    title_fontsize=14,
    legend_fontsize=9,
    show_legend=True,
):
    """Plot multiple permulation histograms overlaid on one axis.

    series_specs should be an iterable of (series_label, values, true_mean, color).
    """

    fig, ax = plt.subplots(figsize=fig_size)

    if xlim is None:
        all_values = np.concatenate([np.asarray(spec[1]) for spec in series_specs])
        xlim = (float(np.min(all_values)), float(np.max(all_values)))

    bins = np.arange(xlim[0], xlim[1] + binwidth, binwidth)

    for label, values, true_mean, color in series_specs:
        ax.hist(
            values,
            bins=bins,
            histtype="stepfilled",
            linewidth=1.3,
            color=color,
            edgecolor=color,
            alpha=0.22,
            label=label,
        )
        ax.axvline(
            true_mean,
            color=color,
            linestyle="--",
            linewidth=1.8,
            alpha=0.95,
            label=f"{label} true mean",
        )

    ax.set_xlim(xlim)
    ax.set_xlabel("Means", fontsize=axis_label_fontsize, fontweight="bold")
    ax.set_ylabel("Count", fontsize=axis_label_fontsize, fontweight="bold")
    ax.set_title(title, fontsize=title_fontsize)
    ax.tick_params(axis="both", labelsize=tick_fontsize)

    if show_legend:
        handles, labels = ax.get_legend_handles_labels()
        hist_handles, hist_labels = [], []
        mean_handles, mean_labels = [], []
        for handle, label in zip(handles, labels):
            if label.endswith(" true mean"):
                mean_handles.append(handle)
                mean_labels.append(label)
            else:
                hist_handles.append(handle)
                hist_labels.append(label)

        spacer = Line2D([], [], linestyle="none")
        ordered_handles = hist_handles + [spacer] + mean_handles
        ordered_labels = hist_labels + [""] + mean_labels

        ax.legend(
            ordered_handles,
            ordered_labels,
            fontsize=legend_fontsize,
            ncol=1,
            frameon=True,
            handlelength=1.8,
            labelspacing=0.5,
            borderpad=0.6,
            loc="upper right",
        )

    fig.tight_layout()
    return fig, ax

def plot_permulation_stats(
    results,
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
    show_legend=True,
):
    # Normalize test parameter: allow "loss", "dup", or "duplication"
    if test == "duplication":
        test = "dup"
    elif test not in ("loss", "dup"):
        raise ValueError(f"Invalid test type: {test}. Must be 'loss', 'dup', or 'duplication'.")
    
    # Use "dup" internally for attribute lookups, "duplication" for display
    test_name_display = "duplication" if test == "dup" else test
    
    ncols = 2 if include_stddev else 1
    fig_width = 12 if include_stddev else 6.5


    fig, axs = plt.subplots(1, ncols, figsize=(fig_width, 5))
    axs = np.atleast_1d(axs)
    means = getattr(results, f"means_{test}")
    stddevs = getattr(results, f"stddevs_{test}")
    true_mean = getattr(results, f"true_mean_{test}")
    true_stddev = getattr(results, f"true_stddev_{test}")
    if test == "loss":
        maximum = results.max_occ
        if binwidth is None:
            binwidth = 0.05
    else:
        maximum = results.true_odds.total_species_count
        if binwidth is None:
            binwidth = 0.01

    if title:
        fig.suptitle(
            f"Permulated (null) distribution stats for gene {test_name_display},\n"
            f"{fg_name} vs. {bg_name}\n"
            f"Maximum occupancy = {maximum}, minimum occupancy = {results.min_occ}",
            fontsize=16,
        )

    plot_stat_histogram(axs[0], means, binwidth, hist_color, hist_alpha, edgecolor)

    if subplot_titles:
        axs[0].set_title("permulated means")
    axs[0].axvline(
        x=getattr(results, f"{test}_mean_av"),
        linestyle="dotted",
        color="black",
        label="Avg. permulated mean",
    )
    axs[0].axvline(
        x=true_mean,
        linestyle="--",
        color="salmon",
        label="True mean",
    )
    style_stat_axes(
        axs[0],
        xlabel="Means",
        axis_label_fontsize=axis_label_fontsize,
        legend_fontsize=legend_fontsize,
        xlim=xlim,
        ylim=ylim,
        show_legend=show_legend,
    )

    if include_stddev:
        plot_stat_histogram(axs[1], stddevs, binwidth, hist_color, hist_alpha, edgecolor)

        if subplot_titles:
            axs[1].set_title("Standard deviations")

        axs[1].axvline(
            x=getattr(results, f"{test}_stddev_av"),
            linestyle="dotted",
            color="black",
            label="Avg. permulated stddev",
        )
        axs[1].axvline(
            x=true_stddev,
            linestyle="--",
            color="salmon",
            label="True stddev",
        )
        style_stat_axes(
            axs[1],
            xlabel="Standard deviations",
            axis_label_fontsize=axis_label_fontsize,
            legend_fontsize=legend_fontsize,
            xlim=xlim,
            ylim=ylim,
            show_legend=show_legend,
        )

    plt.tight_layout()
    return fig, axs


def plot_permulation_results(
    results,
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
    show_legend=True,
    show_textbox=True,
):
    # Normalize test parameter: allow "loss", "dup", or "duplication"
    if test == "duplication":
        test = "dup"
    elif test not in ("loss", "dup"):
        raise ValueError(f"Invalid test type: {test}. Must be 'loss', 'dup', or 'duplication'.")
    
    # Use "dup" internally for attribute lookups, "duplication" for display
    test_name_display = "duplication" if test == "dup" else test
    
    # Compute maximum occupancy based on test type
    if test == "loss":
        maximum = results.max_occ
    else:
        maximum = results.true_odds.total_species_count
    
    fig, ax = plt.subplots(figsize=(6, 5))
    true_vals, x, avg_pdf, true_pdf, y_max = get_permulation_plot_data(results, test, bins=bins)

    if title:
        fig.suptitle(
            f"Log odds ratio of gene {test_name_display}, {fg_name} vs. {bg_name}\n"
            f"Maximum occupancy = {maximum}, minimum occupancy = {results.min_occ}",
            fontsize=14,
        )

    plot_true_histogram(ax, true_vals, bins, hist_color, hist_alpha, edgecolor)
    plot_true_gaussian_fit(ax, x, true_pdf, gaussfit_color)
    plot_average_permulation_curve(
        ax,
        x,
        avg_pdf,
        avpermulation_color,
        linestyle="--",
        fill=False,
    )
    plot_threshold_lines(results, ax, test, thresholds_color)

    style_density_axes(
        ax,
        x,
        y_max,
        legend_fontsize=legend_fontsize,
        axis_label_fontsize=axis_label_fontsize,
        show_legend=show_legend,
    )

    if show_textbox:
        legend = ax.get_legend()
        text_bbox = {"boxstyle": "round,pad=0.2"}
        if legend is not None:
            frame = legend.get_frame()
            frame_alpha = frame.get_alpha()
            if frame_alpha is None:
                frame_alpha = frame.get_facecolor()[-1]
            text_bbox.update(
                {
                    "facecolor": frame.get_facecolor(),
                    "edgecolor": frame.get_edgecolor(),
                    "linewidth": frame.get_linewidth(),
                    "alpha": frame_alpha,
                }
            )
        else:
            text_bbox.update(
                {
                    "facecolor": "white",
                    "edgecolor": "0.8",
                    "linewidth": 1.0,
                    "alpha": 0.8,
                }
            )

        ax.text(
            0.03,
            0.95,
            f"Permulated mean = {format_stat(getattr(results, f'{test}_mean_av'))}\n"
            f"True mean = {format_stat(getattr(results, f'true_mean_{test}'))}\n\n"
            f"Permulated std. dev. = {format_stat(getattr(results, f'{test}_stddev_av'))}\n"
            f"True std. dev. = {format_stat(getattr(results, f'true_stddev_{test}'))}",
            transform=ax.transAxes,
            fontsize=textbox_fontsize,
            ha="left",
            va="top",
            bbox=text_bbox,
            zorder=5,
        )

    plt.tight_layout()
    return fig, ax


def plot_permulation_results_layered(
    results,
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
    show_legend=True,
    show_textbox=True,
):
    # Normalize test parameter: allow "loss", "dup", or "duplication"
    if test == "duplication":
        test = "dup"
    elif test not in ("loss", "dup"):
        raise ValueError(f"Invalid test type: {test}. Must be 'loss', 'dup', or 'duplication'.")
    
    # Use "dup" internally for attribute lookups, "duplication" for display
    test_name_display = "duplication" if test == "dup" else test
    
    # Compute maximum occupancy based on test type
    if test == "loss":
        maximum = results.max_occ
    else:
        maximum = results.true_odds.total_species_count
    
    true_vals, x, avg_pdf, true_pdf, y_max = get_permulation_plot_data(results, test, bins=bins)

    title_str = (
        f"Log odds ratio of gene {test_name_display}, {fg_name} vs. {bg_name}\n"
        f"Maximum occupancy = {maximum}, minimum occupancy = {results.min_occ}"
    )

    figs = []
    axes = []

    fig1, ax1 = plt.subplots(figsize=(8, 6))
    if title:
        fig1.suptitle(title_str, fontsize=14)
    plot_average_permulation_curve(ax1, x, avg_pdf, avpermulation_color)
    style_density_axes(
        ax1,
        x,
        y_max,
        legend_fontsize=legend_fontsize,
        axis_label_fontsize=axis_label_fontsize,
        show_legend=show_legend,
    )
    plt.tight_layout()
    figs.append(fig1)
    axes.append(ax1)

    fig2, ax2 = plt.subplots(figsize=(8, 6))
    if title:
        fig2.suptitle(title_str, fontsize=14)
    plot_average_permulation_curve(ax2, x, avg_pdf, avpermulation_color)
    plot_threshold_lines(results, ax2, test, thresholds_color)
    if show_textbox:
        ax2.text(
            0.03,
            0.95,
            f"Permulated mean = {format_stat(getattr(results, f'{test}_mean_av'))}\n"
            f"Permulated std. dev. = {format_stat(getattr(results, f'{test}_stddev_av'))}",
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
    style_density_axes(
        ax2,
        x,
        y_max,
        legend_fontsize=legend_fontsize,
        axis_label_fontsize=axis_label_fontsize,
        show_legend=show_legend,
    )
    plt.tight_layout()
    figs.append(fig2)
    axes.append(ax2)

    fig3, ax3 = plt.subplots(figsize=(8, 6))
    if title:
        fig3.suptitle(title_str, fontsize=14)
    plot_average_permulation_curve(ax3, x, avg_pdf, avpermulation_color)
    plot_threshold_lines(results, ax3, test, thresholds_color)
    plot_true_histogram(ax3, true_vals, bins, hist_color, 0.3, hist_color)
    if show_textbox:
        ax3.text(
            0.03,
            0.95,
            f"Permulated mean = {format_stat(getattr(results, f'{test}_mean_av'))}\n"
            f"Permulated std. dev. = {format_stat(getattr(results, f'{test}_stddev_av'))}",
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
    style_density_axes(
        ax3,
        x,
        y_max,
        legend_fontsize=legend_fontsize,
        axis_label_fontsize=axis_label_fontsize,
        show_legend=show_legend,
    )
    plt.tight_layout()
    figs.append(fig3)
    axes.append(ax3)

    fig4, ax4 = plt.subplots(figsize=(8, 6))
    if title:
        fig4.suptitle(title_str, fontsize=14)
    plot_average_permulation_curve(ax4, x, avg_pdf, avpermulation_color)
    plot_threshold_lines(results, ax4, test, thresholds_color)
    plot_true_histogram(ax4, true_vals, bins, hist_color, 0.3, hist_color)
    plot_true_gaussian_fit(ax4, x, true_pdf, gaussfit_color)
    if show_textbox:
        ax4.text(
            0.03,
            0.95,
            f"Permulated mean = {format_stat(getattr(results, f'{test}_mean_av'))}\n"
            f"Permulated std. dev. = {format_stat(getattr(results, f'{test}_stddev_av'))}\n\n"
            f"True mean = {format_stat(getattr(results, f'true_mean_{test}'))}\n"
            f"True std. dev. = {format_stat(getattr(results, f'true_stddev_{test}'))}",
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
    style_density_axes(
        ax4,
        x,
        y_max,
        legend_fontsize=legend_fontsize,
        axis_label_fontsize=axis_label_fontsize,
        show_legend=show_legend,
    )
    plt.tight_layout()
    figs.append(fig4)
    axes.append(ax4)

    return figs, axes
