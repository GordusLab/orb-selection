"""Plotting helpers for the odds ratio permulation workflow."""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import plot_split_helpers as psh
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


def _draw_y_axis_break_marks(ax_top, ax_bottom, dx=0.012, angle_deg=55.0):
    """Compatibility wrapper around shared split-axis break mark helper."""
    psh.draw_y_axis_break_marks(ax_top, ax_bottom, dx=dx, angle_deg=angle_deg)


def _compute_shared_split_tick_step(lower_ylim, upper_ylim, tick_count=4):
    """Compatibility wrapper around shared split-axis tick-step helper."""
    return psh.compute_shared_split_tick_step(
        lower_ylim,
        upper_ylim,
        tick_count=tick_count,
    )


def _split_axis_ticks(lower_ylim, upper_ylim, tick_count=4):
    """Compatibility wrapper around shared split-axis tick generator."""
    return psh.split_axis_ticks(
        lower_ylim,
        upper_ylim,
        tick_count=tick_count,
    )


def _resolve_split_height_ratios(lower_ylim, upper_ylim, split_y_height_ratios):
    """Compatibility wrapper around shared split-axis panel-height helper."""
    return psh.resolve_split_height_ratios(
        lower_ylim,
        upper_ylim,
        split_y_height_ratios,
    )


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
    split_y_axis=False,
    split_y_lims=None,
    split_y_height_ratios="auto",
    split_y_tick_count=4,
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

    if split_y_axis:
        if split_y_lims is None:
            raise ValueError(
                "split_y_lims is required when split_y_axis=True. "
                "Pass ((lower_min, lower_max), (upper_min, upper_max))."
            )
        if len(split_y_lims) != 2:
            raise ValueError(
                "split_y_lims must contain two ranges: "
                "((lower_min, lower_max), (upper_min, upper_max))."
            )

        resolved_height_ratios = _resolve_split_height_ratios(
            split_y_lims[0], split_y_lims[1], split_y_height_ratios
        )

        fig = plt.figure(figsize=(fig_width, 6.5))
        gs = fig.add_gridspec(
            2,
            ncols,
            height_ratios=resolved_height_ratios,
            hspace=0.05,
            wspace=0.25,
        )
        top_axes = []
        bottom_axes = []
        for col in range(ncols):
            ax_top = fig.add_subplot(gs[0, col])
            ax_bottom = fig.add_subplot(gs[1, col], sharex=ax_top)
            top_axes.append(ax_top)
            bottom_axes.append(ax_bottom)
        axs = np.atleast_1d(bottom_axes)
    else:
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

    if split_y_axis:
        means_axes = (top_axes[0], bottom_axes[0])
        for ax in means_axes:
            plot_stat_histogram(ax, means, binwidth, hist_color, hist_alpha, edgecolor)
            ax.axvline(
                x=getattr(results, f"{test}_mean_av"),
                linestyle="dotted",
                color="black",
                label="Avg. permulated mean",
            )
            ax.axvline(
                x=true_mean,
                linestyle="--",
                color="salmon",
                label="True mean",
            )

        ax_top, ax_bottom = means_axes
        if subplot_titles:
            ax_top.set_title("permulated means")
        ax_bottom.set_xlabel("Means", fontsize=axis_label_fontsize, fontweight="bold")
        ax_bottom.set_ylabel("Count", fontsize=axis_label_fontsize, fontweight="bold")
        if xlim is not None:
            ax_top.set_xlim(xlim)
            ax_bottom.set_xlim(xlim)
        ax_bottom.set_ylim(split_y_lims[0])
        ax_top.set_ylim(split_y_lims[1])

        lower_ticks, upper_ticks = _split_axis_ticks(
            split_y_lims[0], split_y_lims[1], tick_count=split_y_tick_count
        )
        ax_bottom.set_yticks(lower_ticks)
        ax_top.set_yticks(upper_ticks)

        ax_top.set_ylabel("")
        if show_legend:
            ax_top.legend(fontsize=legend_fontsize)
        ax_top.spines["bottom"].set_visible(False)
        ax_bottom.spines["top"].set_visible(False)
        ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        ax_bottom.tick_params(axis="x", which="both", top=False)
        _draw_y_axis_break_marks(ax_top, ax_bottom)
    else:
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
        if split_y_axis:
            std_axes = (top_axes[1], bottom_axes[1])
            for ax in std_axes:
                plot_stat_histogram(ax, stddevs, binwidth, hist_color, hist_alpha, edgecolor)
                ax.axvline(
                    x=getattr(results, f"{test}_stddev_av"),
                    linestyle="dotted",
                    color="black",
                    label="Avg. permulated stddev",
                )
                ax.axvline(
                    x=true_stddev,
                    linestyle="--",
                    color="salmon",
                    label="True stddev",
                )

            ax_top, ax_bottom = std_axes
            if subplot_titles:
                ax_top.set_title("Standard deviations")
            ax_bottom.set_xlabel(
                "Standard deviations", fontsize=axis_label_fontsize, fontweight="bold"
            )
            ax_bottom.set_ylabel("Count", fontsize=axis_label_fontsize, fontweight="bold")
            if xlim is not None:
                ax_top.set_xlim(xlim)
                ax_bottom.set_xlim(xlim)
            ax_bottom.set_ylim(split_y_lims[0])
            ax_top.set_ylim(split_y_lims[1])

            lower_ticks, upper_ticks = _split_axis_ticks(
                split_y_lims[0], split_y_lims[1], tick_count=split_y_tick_count
            )
            ax_bottom.set_yticks(lower_ticks)
            ax_top.set_yticks(upper_ticks)

            ax_top.set_ylabel("")
            if show_legend:
                ax_top.legend(fontsize=legend_fontsize)
            ax_top.spines["bottom"].set_visible(False)
            ax_bottom.spines["top"].set_visible(False)
            ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
            ax_bottom.tick_params(axis="x", which="both", top=False)
            _draw_y_axis_break_marks(ax_top, ax_bottom)
        else:
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

    if split_y_axis:
        fig.subplots_adjust(top=0.92, bottom=0.12)
    else:
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
