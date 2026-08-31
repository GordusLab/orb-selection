import numpy as np
from matplotlib import patheffects as pe
from matplotlib import pyplot as plt
from matplotlib import ticker

OMEGA_COLORS = ('salmon', 'steelblue', 'goldenrod')
OMEGA_STROKE_COLORS = ('darkred', 'darkblue', 'brown')
MIN_VISIBLE_PROPORTION = 0.006
TICK_EDGE_PADDING = 10 ** 0.12


def _visible_proportion(proportion):
    """Keep nonzero bars visible without representing zero as nonzero."""
    return max(proportion, MIN_VISIBLE_PROPORTION) if proportion > 0 else 0


def _plot_tiny_proportion_marker(ax, omega, proportion, color, alpha):
    """Mark bars that were raised to the minimum visible height."""
    if 0 < proportion < MIN_VISIBLE_PROPORTION:
        ax.plot(
            omega,
            MIN_VISIBLE_PROPORTION,
            marker='_',
            markersize=20,
            markeredgewidth=2,
            color=color,
            alpha=alpha,
            zorder=10,
        )
        print(f"plotting marker: omega={omega}, proportion={proportion}")


def _plot_distributions(ax, x, group, logbins, labels=False):
    for index, (color, stroke_color) in enumerate(
        zip(OMEGA_COLORS, OMEGA_STROKE_COLORS), start=1
        ):
        histogram_kwargs = {
            'bins': logbins,
            'histtype': 'stepfilled',
            'color': color,
            'alpha': 0.17,
            'path_effects': [pe.Stroke(linewidth=1, foreground=stroke_color)],
        }
        if labels:
            histogram_kwargs['label'] = f'ω{index} distribution'
        ax.hist(x[f'ω{index}_{group}'], **histogram_kwargs)


def _plot_mean_lines(ax, x, group, inverted=False, stroke_width=7.5):
    means = [x[f'ω{index}_{group}'].mean() for index in range(1, 4)]
    weights = [x[f'ω{index}_{group}_P'].mean() for index in range(1, 4)]

    for index, (mean, weight, color) in enumerate(
        zip(means, weights, OMEGA_COLORS), start=1
        ):
        line_kwargs = {
            'linewidth': 6,
            'color': color,
            'path_effects': [
                pe.withStroke(linewidth=stroke_width, foreground='white'),
                pe.Normal(),
            ],
        }
        if inverted:
            line_kwargs['ymin'] = 1 - weight
            line_kwargs['label'] = f'mean inferred ω{index}'
        else:
            line_kwargs['ymax'] = weight
        ax.axvline(mean, **line_kwargs)


def _style_mean_axis(ax):
    ax.yaxis.set_label_position('left')
    ax.yaxis.tick_left()
    for tick in ax.get_yticklabels():
        tick.set_fontweight('bold')
        tick.set_fontsize(14)


def _add_secondary_ylabel(fig):
    # Use a transparent overlay axis to place an independent right-side label.
    dummy_ax = fig.add_subplot(1, 1, 1)
    dummy_ax.set_xticks([])
    dummy_ax.set_yticks([])
    for side in ('left', 'top', 'right', 'bottom'):
        dummy_ax.spines[side].set_visible(False)
    dummy_ax.patch.set_visible(False)
    dummy_ax.yaxis.set_label_position('right')
    dummy_ax.set_ylabel(
        'number of orthogroups',
        labelpad=48,
        rotation=270,
        color='whitesmoke',
        alpha=0.6,
        fontsize=16,
        path_effects=[pe.withStroke(linewidth=1, foreground='black')],
    )


def plot_omega_distributions(
    x,
    result,
    top_title,
    bottom_title,
    numeral="",
    suptitle=None,
    filename=None,
    xlim_max=100000,
    xlim_min=0.001,
    shift_top_title=False,
    transparent=True,
):
    """Plot ω distributions and mean site proportions for test and reference."""
    logbins = np.geomspace(xlim_min, xlim_max, 100)

    y_limits = []

    for group in ['test', 'ref']:
        ω3_counts = np.histogram(x[f'ω3_{group}'], bins=logbins)[0]
        ω2_counts = np.histogram(x[f'ω2_{group}'], bins=logbins)[0]

        ymax = np.max([np.max(ω3_counts), np.max(ω2_counts)])
        pad = ymax / 10

        y_limits.append(ymax + 0.75 * pad)

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(6, 5))

    plt.subplots_adjust(hspace=0)
    plt.xlim(xlim_min, xlim_max)
    plt.xscale('log')
    plt.rcParams['font.family'] = 'Verdana'

    # Match histogram scales so the mirrored distributions remain comparable.
    axs[0].set_ylim(0, np.max(y_limits))
    axs[1].set_ylim(0, np.max(y_limits))
    axs[1].invert_yaxis()

    _plot_distributions(axs[0], x, 'test', logbins)
    _plot_distributions(axs[1], x, 'ref', logbins, labels=True)

    ax_avgs = axs[0].twinx()
    ax_avgs.set_facecolor('none')
    _plot_mean_lines(
        ax_avgs, x, 'test', stroke_width=7.5 if result == 'busted' else 7
    )
    ax_avgs.spines.right.set_visible(False)

    ax_avgs2 = axs[1].twinx()
    ax_avgs2.set_facecolor('none')
    ax_avgs2.invert_yaxis()
    _plot_mean_lines(
        ax_avgs2,
        x,
        'ref',
        inverted=True,
        stroke_width=7 if result == 'busted' else 7.5,
    )
    ax_avgs2.spines.right.set_visible(False)

    ax_avgs.axvline(1, linewidth=0.5, linestyle='dashed', color='k', alpha=0.5)
    ax_avgs2.axvline(1, linewidth=0.5, linestyle='dashed', color='k', alpha=0.5)

    if result == 'all':
        ax_avgs.set_title(
            top_title, x=0.2, y=0.7, fontsize=13, color='white',
            weight='bold', backgroundcolor='lightgray',
        )
        ax_avgs2.set_title(
            bottom_title, x=0.2, y=0.1, fontsize=13, color='white',
            weight='bold', backgroundcolor='lightgray',
        )
    elif result != 'busted':
        bottom_title_x = 0.185 if result in ('relaxed', 'intensified') else 0.2
        ax_avgs.set_title(
            top_title, x=0.75, y=0.7, fontsize=13, color='white',
            weight='bold', backgroundcolor='lightgray',
        )
        ax_avgs2.set_title(
            bottom_title, x=bottom_title_x, y=0.075, fontsize=13,
            color='silver', weight='bold', backgroundcolor='white',
        )

        if result == 'intensified':
            plt.text(
                0.8, -0.8, '$\\it{p}$≤0.05\n$\\it{k}$>1', fontsize=15,
                ha='left', va='center', transform=ax_avgs.transAxes,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'),
            )
        else:
            plt.text(
                0.8, -0.8, '$\\it{p}$≤0.05\n$\\it{k}$<1', fontsize=15,
                ha='left', va='center', transform=ax_avgs.transAxes,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'),
            )
    else:
        if shift_top_title:
            ax_avgs.set_title(
                top_title, x=0.8, y=0.6, fontsize=13, color='white',
                weight='bold', backgroundcolor='lightgray',
            )
        else:
            ax_avgs.set_title(
                top_title, x=0.8, y=0.7, fontsize=13, color='white',
                weight='bold', backgroundcolor='lightgray',
            )
        ax_avgs2.set_title(
            bottom_title, x=0.185, y=0.1, fontsize=13, color='silver',
            weight='bold', backgroundcolor='white',
        )

        plt.text(
            0.775, -0.75,
            '$\\it{p_1}$≤0.05\n$\\it{p_2}$>0.05\n$\\it{p_3}$≤0.05',
            fontsize=15, ha='left', va='center', transform=ax_avgs.transAxes,
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'),
        )

    for ax in [ax_avgs, ax_avgs2]:
        _style_mean_axis(ax)

    for axis in axs:
        axis.yaxis.set_label_position('right')
        axis.yaxis.tick_right()
        for tick in axis.get_yticklabels():
            tick.set_fontsize(14)

    for tick in axs[1].get_xticklabels():
        tick.set_fontsize(14)

    fig.supylabel('mean proportion of sites', x=0.01, weight='bold', fontsize=16)
    _add_secondary_ylabel(fig)

    for ax in fig.axes:
        ax.tick_params(direction='in')
        ax.tick_params(axis='x', which='major', pad=10)

    if suptitle is not None:
        y_pos = 0.96 if result == 'all' else 0.98
        fig.suptitle(f"{numeral}{suptitle}", y=y_pos, fontsize=16)

    if filename is not None:
        plt.savefig(filename, dpi=600, transparent=transparent, bbox_inches='tight')

    return fig, axs


def plot_omega_single_gene(
    df,
    gene,
    suptitle=None,
    subtitle=None,
    i="",
    offset_zero_w1=False,
    offset_zero_w2=False,
    k=False,
    filename=None,
    transparent=True,
    build_in=False,
    xlim_min=None,
    xlim_max=None,
    zorder_1ref=1,
    zorder_2ref=1,
    zorder_3ref=1,
    zorder_1test=1,
    zorder_2test=1,
    zorder_3test=1,
    axis_fontsize=14,
):
    """Plot inferred ω values and site proportions for one gene."""

    x = df.loc[gene]

    fig, ax = plt.subplots(figsize=(5.5, 5))

    # Collapse numerical near-zero estimates so they share the zero position.
    if x['ω1_ref'] < 1e-9:
        x['ω1_ref'] = 0

    if x['ω1_test'] < 1e-9:
        x['ω1_test'] = 0

    omega_values = np.asarray([
        x['ω1_ref'], x['ω2_ref'], x['ω3_ref'],
        x['ω1_test'], x['ω2_test'], x['ω3_test'],
    ])
    has_zero = np.any(omega_values == 0)

    if np.min([x['ω2_ref'], x['ω2_test']]) > 0:
        if np.min([x['ω1_ref'], x['ω1_test']]) > 0:
            thresh = np.min([x['ω1_ref'], x['ω1_test']]) * 0.1
        elif np.max([x['ω1_ref'], x['ω1_test']]) > 0:
            thresh = np.max([x['ω1_ref'], x['ω1_test']]) * 0.1
        else:
            thresh = np.min([x['ω2_ref'], x['ω2_test']]) * 0.1
    elif np.max([x['ω2_ref'], x['ω2_test']]) > 0:
        thresh = np.max([x['ω2_ref'], x['ω2_test']]) * 0.1
    else:
        thresh = 0.01


    plt.subplots_adjust(hspace=0)
    if has_zero:
        plt.xscale('symlog', linthresh=thresh)
    else:
        plt.xscale('log')
    ax.set_ylim(0, 1)
    plt.rcParams['font.family'] = 'Verdana'

    ref_path_effects_w1_offset = [
        pe.SimpleLineShadow(offset=(-2, 0), alpha=0.3, foreground='salmon'),
        pe.Normal(),
    ]
    ref_path_effects_w2_offset = [
        pe.SimpleLineShadow(offset=(-2, 0), alpha=0.3, foreground='steelblue'),
        pe.Normal(),
    ]
    ref_path_effects = [
        pe.withStroke(linewidth=22, foreground='white', offset=(-0.5, 1)),
        pe.Normal(),
    ]

    ax.vlines(
        x['ω1_ref'],
        0,
        _visible_proportion(x['ω1_ref_P']),
        linewidth=20,
        color='salmon',
        alpha=0.17,
        zorder=zorder_1ref,
        path_effects=(
            ref_path_effects if not offset_zero_w1 else ref_path_effects_w1_offset
        ),
    )
    _plot_tiny_proportion_marker(
        ax, x['ω1_ref'], x['ω1_ref_P'], 'salmon', alpha=0.17
    )
    ax.vlines(
        x['ω2_ref'],
        0,
        _visible_proportion(x['ω2_ref_P']),
        linewidth=20,
        color='steelblue',
        alpha=0.17,
        zorder=zorder_2ref,
        path_effects=(
            ref_path_effects if not offset_zero_w2 else ref_path_effects_w2_offset
        ),
    )
    _plot_tiny_proportion_marker(
        ax, x['ω2_ref'], x['ω2_ref_P'], 'steelblue', alpha=0.17
    )
    ax.vlines(
        x['ω3_ref'],
        0,
        _visible_proportion(x['ω3_ref_P']),
        linewidth=20,
        color='goldenrod',
        alpha=0.17,
        zorder=zorder_3ref,
        path_effects=ref_path_effects
    )
    _plot_tiny_proportion_marker(
        ax, x['ω3_ref'], x['ω3_ref_P'], 'goldenrod', alpha=0.3
    )

    test_alpha = 1
    test_effects = [] if build_in else [
        pe.withStroke(linewidth=22, foreground='white', offset=(-0.5, 1))
    ]

    ax.vlines(
        x['ω1_test'],
        0,
        _visible_proportion(x['ω1_test_P']),
        linewidth=20,
        color='salmon',
        alpha=test_alpha,
        zorder=zorder_1test,
        path_effects=test_effects,
    )
    _plot_tiny_proportion_marker(
        ax, x['ω1_test'], x['ω1_test_P'], 'salmon', alpha=1
    )
    ax.vlines(
        x['ω2_test'],
        0,
        _visible_proportion(x['ω2_test_P']),
        linewidth=20,
        color='steelblue',
        alpha=test_alpha,
        zorder=zorder_2test,
        path_effects=test_effects,
    )
    _plot_tiny_proportion_marker(
        ax, x['ω2_test'], x['ω2_test_P'], 'steelblue', alpha=1
    )
    ax.vlines(
        x['ω3_test'],
        0,
        _visible_proportion(x['ω3_test_P']),
        linewidth=20,
        color='goldenrod',
        alpha=test_alpha,
        zorder=zorder_3test,
        path_effects=test_effects,
    )
    _plot_tiny_proportion_marker(
        ax, x['ω3_test'], x['ω3_test_P'], 'goldenrod', alpha=1
    )

    if not has_zero:
        lowest_tick = 10 ** np.floor(np.log10(omega_values.min()))
        highest_tick = 10 ** np.ceil(np.log10(omega_values.max()))
        ax.set_xlim(
            left=xlim_min if xlim_min is not None else lowest_tick / TICK_EDGE_PADDING,
            right=xlim_max if xlim_max is not None else highest_tick * TICK_EDGE_PADDING,
        )
    else:
        highest_tick = 10 ** np.ceil(
            np.log10(omega_values[omega_values > 0].max())
        )
        ax.set_xlim(
            left=min(xlim_min, -thresh * 0.25) if xlim_min is not None else -thresh * 0.25,
            right=xlim_max if xlim_max is not None else highest_tick * TICK_EDGE_PADDING,
        )

    ax.axvline(1, linewidth=1, linestyle='dashed', color='k', alpha=0.5)
    if not has_zero:
        log_ticks = 10 ** np.arange(
            np.log10(lowest_tick), np.log10(highest_tick) + 1
        )
        ax.xaxis.set_major_locator(ticker.FixedLocator(log_ticks))
        ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
        ax.xaxis.set_minor_locator(ticker.NullLocator())
    if has_zero:
        positive_omegas = omega_values[omega_values > 0]
        lowest_tick = 10 ** np.ceil(
            np.log10(max(positive_omegas.min(), thresh))
        )
        log_ticks = 10 ** np.arange(
            np.log10(lowest_tick), np.log10(highest_tick) + 1
        )
        ax.xaxis.set_major_locator(ticker.FixedLocator(np.r_[0, log_ticks]))
        ax.xaxis.set_major_formatter(
            ticker.FuncFormatter(
                lambda value, _: '0' if value == 0 else f'$10^{{{int(np.log10(value))}}}$'
            )
        )
        ax.xaxis.set_minor_locator(ticker.NullLocator())
        zero_x = ax.transAxes.inverted().transform(
            ax.transData.transform((0, 0))
        )[0]
        first_tick_x = ax.transAxes.inverted().transform(
            ax.transData.transform((lowest_tick, 0))
        )[0]
        break_x = zero_x + (first_tick_x - zero_x) * 0.35
        for offset in (-0.012, 0.012):
            ax.plot(
                (break_x + offset - 0.007, break_x + offset + 0.007),
                (-0.015, 0.015),
                color='black', linewidth=1, transform=ax.transAxes,
                clip_on=False,
            )
    ax.tick_params(axis='x', which='both', zorder=3)

    for tick in ax.get_yticklabels():
        tick.set_fontweight('bold')
        tick.set_fontsize(14)

    for tick in ax.get_xticklabels():
        tick.set_fontsize(axis_fontsize)

    fig.supylabel('proportion of sites', x=0.01, weight='bold', fontsize=16)

    for ax in fig.axes:
        ax.tick_params(direction='in')

    if suptitle is not None:
        if subtitle is not None:
            fig.suptitle(f"{i}{suptitle}", y=0.96, fontsize=14, fontweight='bold')
        else:
            fig.suptitle(f"{i}{suptitle}", y=0.92, fontsize=14, fontweight='bold')

    if not build_in:
        if subtitle is not None:
            plt.title(subtitle, fontsize=12)

        if k:
            plt.text(
                0.7, 0.9,
                f'$\\it{{k}}$={round(x["k"], 2)}\n$\\it{{p}}$={x["p_value"]:.2e}',
                fontsize=13, ha='left', va='center', transform=ax.transAxes,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'),
            )
        else:
            plt.text(
                0.65, 0.9,
                f'$\\it{{p_1}}$={x["test_pval"]:.2e}\n'
                f'$\\it{{p_3}}$={x["shared_pval"]:.2e}',
                fontsize=13, ha='left', va='center', transform=ax.transAxes,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'),
            )

    if filename is not None:
        plt.savefig(filename, dpi=300, transparent=transparent)

    print(
        f"ω values for {gene}: ω1_test={x['ω1_test']:.4f}, "
        f"ω2_test={x['ω2_test']:.4f}, ω3_test={x['ω3_test']:.4f},\n"
        f"ω1_ref={x['ω1_ref']:.4f}, ω2_ref={x['ω2_ref']:.4f}, "
        f"ω3_ref={x['ω3_ref']:.4f};\n"
        f"Proportions: ω1_test_P={x['ω1_test_P']:.4f}, "
        f"ω2_test_P={x['ω2_test_P']:.4f}, "
        f"ω3_test_P={x['ω3_test_P']:.4f},\n"
        f"ω1_ref_P={x['ω1_ref_P']:.4f}, ω2_ref_P={x['ω2_ref_P']:.4f}, "
        f"ω3_ref_P={x['ω3_ref_P']:.4f}"
    )

    return fig, ax
