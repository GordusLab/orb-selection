import numpy as np
from matplotlib import patheffects as pe
from matplotlib import pyplot as plt
from matplotlib import ticker
from matplotlib.patches import Rectangle

OMEGA_COLORS = ('salmon', 'steelblue', 'goldenrod')
OMEGA_STROKE_COLORS = ('darkred', 'darkblue', 'brown')
MIN_VISIBLE_PROPORTION = 0.006
TICK_EDGE_PADDING = 10 ** 0.25


SAVEFIG_DPI = 300


def _min_top_gap(ax, dpi):
    """Data-unit height needed for a visible gap below the top spine at the given dpi.

    The top spine itself is drawn ~centered on ymax and has nonzero width, so a bar
    reaching to within one bare pixel of ymax still visually touches the spine. Clear
    half the spine's rendered width plus one pixel of actual white space.
    """
    pos = ax.get_position()
    pixel_height = pos.height * ax.get_figure().get_figheight() * dpi
    if pixel_height <= 0:
        return 0
    spine_px = ax.spines['top'].get_linewidth() * dpi / 72
    margin_px = spine_px / 2 + 1
    ymin, ymax = ax.get_ylim()
    return (ymax - ymin) * margin_px / pixel_height


def _visible_proportion(proportion, ax=None):
    """Keep nonzero bars visible without representing zero, or near-1, as exactly 1.

    Uses the lower of the figure's on-screen dpi and the eventual savefig dpi, so the
    gap stays visible whichever rendering (inline display or saved file) is viewed.
    """
    if proportion <= 0:
        return 0
    height = max(proportion, MIN_VISIBLE_PROPORTION)
    if ax is not None and proportion < 1.0:
        _, ymax = ax.get_ylim()
        dpi = min(ax.get_figure().dpi, SAVEFIG_DPI)
        height = min(height, ymax - _min_top_gap(ax, dpi))
    return height


def _hatched_bar_bounds(ax, x_data, linewidth):
    """Data-unit (left, right) bounds of a vlines bar of the given linewidth (pts).

    `ax.vlines` bars are fixed-width in points regardless of the x-axis scale;
    replicating that width for a Rectangle patch means converting through
    pixels at the bar's actual location rather than using a flat data-unit width.
    """
    dpi = ax.get_figure().dpi
    half_width_px = linewidth * dpi / 72 / 2
    x_px, y_px = ax.transData.transform((x_data, 0))
    x_left = ax.transData.inverted().transform((x_px - half_width_px, y_px))[0]
    x_right = ax.transData.inverted().transform((x_px + half_width_px, y_px))[0]
    return x_left, x_right


def _plot_combined_zero_bar(ax, height, alpha, zorder, path_effects=None):
    """Diagonally hatch an ω=0 bar in both ω1 (salmon) and ω2 (steelblue).

    When both rate classes collapse to the same ω=0 position, drawing them as
    two overlapping solid vlines leaves whichever is drawn second (ω2) fully
    occluding ω1, even when ω1 has a nonzero site proportion. A two-direction
    hatch makes both categories visible at once.

    The white halo (`path_effects`) is rendered once, on an invisible vlines
    at the same position/width -- applying it directly to each hatch
    Rectangle instead would stroke them independently and let one rectangle's
    halo paint over the other's hatch lines.
    """
    if path_effects:
        ax.vlines(0, 0, height, linewidth=20, color='none', zorder=zorder - 0.01,
                   path_effects=path_effects)

    x_left, x_right = _hatched_bar_bounds(ax, 0, linewidth=20)
    width = x_right - x_left
    for edgecolor, hatch in (('salmon', '//'), ('steelblue', '\\\\')):
        ax.add_patch(Rectangle(
            (x_left, 0), width, height,
            facecolor='none', edgecolor=edgecolor, hatch=hatch,
            linewidth=0, alpha=alpha, zorder=zorder,
        ))


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

    print(
        f"mean ω values for {group}: "
        f"ω1={means[0]:.4f}, ω2={means[1]:.4f}, ω3={means[2]:.4f}\n"
        f"mean site proportions for {group}: "
        f"ω1={weights[0]:.4f}, ω2={weights[1]:.4f}, ω3={weights[2]:.4f}"
    )

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


def _build_split_axes(
    x,
    test_group,
    ref_group,
    logbins,
    xlim_min,
    xlim_max,
    split_lower_max=None,
    split_upper_max=None,
):
    """Build a 6-row grid with a broken gap in both the test and ref y-axes.

    Each side (test/ref) gets a "near" panel (close to the break, showing the
    bulk of the ω2/ω3 counts) and a "far" panel (the ω3 tail, which otherwise
    dwarfs everything else on a linear count axis). Returns the figure, the
    6 data/spacer axes, and the two full-height overlay axes used for the
    mean-line twins (ax_avgs/ax_avgs2 in the caller).
    """
    # Derive one pair of y-axis intervals from both distributions. Using
    # separate intervals makes numerically similar breaks look misaligned
    # because their panels receive different physical heights.
    maxima = []
    for group in (test_group, ref_group):
        for rate_class in (2, 3):
            counts = np.histogram(x[f'ω{rate_class}_{group}'], bins=logbins)[0]
            maxima.extend((np.max(counts), np.partition(counts, -2)[-2]))

    ymax = np.max(maxima)
    lower_max = (
        np.partition(maxima, -2)[-2]
        if split_lower_max is None else split_lower_max
    )
    if not 0 < lower_max < ymax:
        raise ValueError(
            'split_lower_max must be greater than zero and below the upper '
            f'panel minimum ({ymax:.2f}).'
        )
    pad = ymax / 10
    round_scale = 10 ** np.floor(np.log10(ymax))
    upper_tick = (
        np.ceil(ymax / round_scale) * round_scale
        if split_upper_max is None else split_upper_max
    )
    upper_min = ymax - pad
    upper_span = upper_tick - upper_min
    upper_interval = (upper_min, upper_tick + 0.2 * upper_span)
    if upper_interval[1] < ymax:
        raise ValueError(
            'split_upper_max and its display padding must include the largest '
            f'histogram-bin count ({ymax:.2f}).'
        )
    lower_interval = (0, lower_max + pad)

    upper_height = upper_interval[1] - upper_interval[0]
    lower_height = lower_interval[1] - lower_interval[0]
    upper_ratio = 4 / (1 + lower_height / upper_height)
    lower_ratio = 4 - upper_ratio

    fig, axs = plt.subplots(
        6, 1, sharex=True,
        height_ratios=[upper_ratio, 0.2, lower_ratio, lower_ratio, 0.2, upper_ratio],
        figsize=(6, 5),
    )

    plt.subplots_adjust(hspace=0)
    plt.xlim(xlim_min, xlim_max)
    plt.xscale('log')
    plt.rcParams['font.family'] = 'Verdana'

    axs[1].set_axis_off()
    axs[4].set_axis_off()

    axs[0].set_ylim(*upper_interval)
    axs[2].set_ylim(*lower_interval)
    axs[3].set_ylim(*lower_interval)
    axs[3].invert_yaxis()
    axs[5].set_ylim(*upper_interval)
    axs[5].invert_yaxis()

    intervals = []
    for ax in (axs[0], axs[2]):
        _plot_distributions(ax, x, test_group, logbins)
        locator = ticker.AutoLocator()
        ax.yaxis.set_major_locator(locator)
        intervals.append(locator()[1] - locator()[0])

    for ax in (axs[3], axs[5]):
        _plot_distributions(ax, x, ref_group, logbins, labels=True)
        locator = ticker.AutoLocator()
        ax.yaxis.set_major_locator(locator)
        intervals.append(locator()[1] - locator()[0])

    tick_step = 100 if any(i >= 100 for i in intervals) else np.max(intervals)
    for j in (2, 3, 5):
        axs[j].yaxis.set_major_locator(ticker.MultipleLocator(tick_step))
    # Each upper panel is intentionally a single round-number outlier tick.
    for j in (0, 5):
        axs[j].yaxis.set_major_locator(ticker.FixedLocator([upper_tick]))

    axs[0].spines.bottom.set_visible(False)
    axs[2].spines.top.set_visible(False)
    axs[3].spines.bottom.set_visible(False)
    axs[5].spines.top.set_visible(False)

    # Overlay axes spanning each half of the grid, for the mean-line twins.
    gs = axs[0].get_gridspec()
    top_overlay = fig.add_subplot(gs[:3], sharex=axs[5])
    bottom_overlay = fig.add_subplot(gs[3:], sharex=axs[5])
    top_overlay.set_facecolor('none')
    top_overlay.set_axis_off()
    bottom_overlay.set_facecolor('none')
    bottom_overlay.set_axis_off()

    axs[0].tick_params(bottom=False, which='both')
    axs[3].tick_params(bottom=False, which='both')

    # Diagonal break markers at the right side of each gap.
    d = 0.1
    marker_kwargs = dict(
        marker=[(-1, -d), (1, d)], markersize=8, linestyle='none',
        color='k', mec='dimgrey', mew=1, clip_on=False,
    )
    for ax, y in ((axs[0], 0), (axs[2], 1), (axs[3], 0), (axs[5], 1)):
        ax.plot((None, 1), (y, y), transform=ax.transAxes, **marker_kwargs)

    hist_axes = [axs[0], axs[2], axs[3], axs[5]]
    return fig, axs, hist_axes, top_overlay, bottom_overlay


def plot_omega_distributions(
    x,
    result,
    top_title,
    bottom_title,
    numeral="",
    suptitle=None,
    filename=None,
    xlim_max=None,
    xlim_min=None,
    shift_top_title=False,
    transparent=True,
    groups=('test', 'ref'),
    split_axis=False,
    partitioned_desc=False,
    split_lower_max=None,
    split_upper_max=None,
):
    """Plot ω distributions and mean site proportions for test and reference.

    `groups` names the test/reference column suffixes to plot, e.g. the default
    ('test', 'ref') for the RELAX/BUSTED-PH alternative model, or
    ('test_pd', 'ref_pd') for the RELAX partitioned descriptive model.
    `partitioned_desc=True` is shorthand for `groups=('test_pd', 'ref_pd')` and
    takes precedence over `groups` if both are given.

    `split_axis` breaks each side's y-axis into a near/far pair of panels
    (with a diagonal-marked gap) so a handful of outlier ω3 counts don't
    flatten the rest of the distribution. Use it when the ω2/ω3 counts span
    a wide enough range that a single linear y-axis can't show both.
    `split_lower_max` optionally sets the maximum y-value displayed below the
    break for both mirrored distributions; otherwise it is inferred.
    `split_upper_max` optionally sets the single, round-number tick and maximum
    y-value shown above the break for both mirrored distributions.

    `xlim_min`/`xlim_max` default to the 0.1st/99.9th percentile of the
    plotted ω2/ω3 values, rounded outward to the nearest power of 10, rather
    than a fixed range or the true min/max -- RELAX in particular can fit a
    rate class to a vanishing proportion of sites with an essentially
    unconstrained ω estimate, producing underflow values (e.g. ~1e-98) that
    are numerically real but not a meaningful plot boundary; a fixed default
    would either clip real data or blow out the axis depending on the dataset,
    and true min/max stretches the whole axis to fit a single such outlier.
    """
    if partitioned_desc:
        groups = ('test_pd', 'ref_pd')
    test_group, ref_group = groups

    omega23 = np.concatenate([
        x[f'ω{index}_{group}'].values
        for index in (2, 3) for group in groups
    ])
    omega23 = omega23[omega23 > 0]
    if xlim_min is None:
        xlim_min = (
            10 ** np.floor(np.log10(np.percentile(omega23, 0.1)))
            if omega23.size else 0.001
        )
    if xlim_max is None:
        xlim_max = (
            10 ** np.ceil(np.log10(np.percentile(omega23, 99.9)))
            if omega23.size else 100000
        )

    logbins = np.geomspace(xlim_min, xlim_max, 100)

    if split_axis:
        fig, axs, hist_axes, top_overlay, bottom_overlay = _build_split_axes(
            x, test_group, ref_group, logbins, xlim_min, xlim_max,
            split_lower_max=split_lower_max,
            split_upper_max=split_upper_max,
        )
        bottom_tick_axis = axs[5]

        ax_avgs = top_overlay.twinx()
        ax_avgs2 = bottom_overlay.twinx()
    else:
        y_limits = []

        for group in groups:
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

        _plot_distributions(axs[0], x, test_group, logbins)
        _plot_distributions(axs[1], x, ref_group, logbins, labels=True)

        hist_axes = axs
        bottom_tick_axis = axs[1]

        ax_avgs = axs[0].twinx()
        ax_avgs2 = axs[1].twinx()

    ax_avgs.set_facecolor('none')
    _plot_mean_lines(
        ax_avgs, x, test_group, stroke_width=7.5 if result == 'busted' else 7
    )
    ax_avgs.spines.right.set_visible(False)

    ax_avgs2.set_facecolor('none')
    ax_avgs2.invert_yaxis()
    _plot_mean_lines(
        ax_avgs2,
        x,
        ref_group,
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
        ax_avgs.set_title(
            top_title, x=0.25, y=0.7, fontsize=13, color='white',
            weight='bold', backgroundcolor='lightgray',
        )
        ax_avgs2.set_title(
            bottom_title, x=0.25, y=0.05, fontsize=13, alpha=0.8,
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

    for axis in hist_axes:
        axis.yaxis.set_label_position('right')
        axis.yaxis.tick_right()
        for tick in axis.get_yticklabels():
            tick.set_fontsize(14)

    for tick in bottom_tick_axis.get_xticklabels():
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
    groups=('test', 'ref'),
    show_subticks=True,
):
    """Plot inferred ω values and site proportions for one gene.

    `groups` names the test/reference column suffixes to plot, e.g. the default
    ('test', 'ref') for the RELAX/BUSTED-PH alternative model, or
    ('test_pd', 'ref_pd') for the RELAX partitioned descriptive model.
    Set `show_subticks=True` to display unlabeled minor ticks at 2–9 within
    each logarithmic decade.
    """
    test_group, ref_group = groups

    x = df.loc[gene]

    fig, ax = plt.subplots(figsize=(5.5, 5))

    # Collapse numerical near-zero estimates so they share the zero position.
    if x[f'ω1_{ref_group}'] < 1e-9:
        x[f'ω1_{ref_group}'] = 0

    if x[f'ω1_{test_group}'] < 1e-9:
        x[f'ω1_{test_group}'] = 0

    omega_values = np.asarray([
        x[f'ω1_{ref_group}'], x[f'ω2_{ref_group}'], x[f'ω3_{ref_group}'],
        x[f'ω1_{test_group}'], x[f'ω2_{test_group}'], x[f'ω3_{test_group}'],
    ])
    has_zero = np.any(omega_values == 0)

    if np.min([x[f'ω2_{ref_group}'], x[f'ω2_{test_group}']]) > 0:
        if np.min([x[f'ω1_{ref_group}'], x[f'ω1_{test_group}']]) > 0:
            thresh = np.min([x[f'ω1_{ref_group}'], x[f'ω1_{test_group}']]) * 0.1
        elif np.max([x[f'ω1_{ref_group}'], x[f'ω1_{test_group}']]) > 0:
            thresh = np.max([x[f'ω1_{ref_group}'], x[f'ω1_{test_group}']]) * 0.1
        else:
            thresh = np.min([x[f'ω2_{ref_group}'], x[f'ω2_{test_group}']]) * 0.1
    elif np.max([x[f'ω2_{ref_group}'], x[f'ω2_{test_group}']]) > 0:
        thresh = np.max([x[f'ω2_{ref_group}'], x[f'ω2_{test_group}']]) * 0.1
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
        x[f'ω1_{ref_group}'],
        0,
        _visible_proportion(x[f'ω1_{ref_group}_P'], ax),
        linewidth=20,
        color='salmon',
        alpha=0.17,
        zorder=zorder_1ref,
        path_effects=(
            ref_path_effects if not offset_zero_w1 else ref_path_effects_w1_offset
        ),
    )
    _plot_tiny_proportion_marker(
        ax, x[f'ω1_{ref_group}'], x[f'ω1_{ref_group}_P'], 'salmon', alpha=0.17
    )
    ax.vlines(
        x[f'ω2_{ref_group}'],
        0,
        _visible_proportion(x[f'ω2_{ref_group}_P'], ax),
        linewidth=20,
        color='steelblue',
        alpha=0.17,
        zorder=zorder_2ref,
        path_effects=(
            ref_path_effects if not offset_zero_w2 else ref_path_effects_w2_offset
        ),
    )
    _plot_tiny_proportion_marker(
        ax, x[f'ω2_{ref_group}'], x[f'ω2_{ref_group}_P'], 'steelblue', alpha=0.17
    )
    ax.vlines(
        x[f'ω3_{ref_group}'],
        0,
        _visible_proportion(x[f'ω3_{ref_group}_P'], ax),
        linewidth=20,
        color='goldenrod',
        alpha=0.17,
        zorder=zorder_3ref,
        path_effects=ref_path_effects
    )
    _plot_tiny_proportion_marker(
        ax, x[f'ω3_{ref_group}'], x[f'ω3_{ref_group}_P'], 'goldenrod', alpha=0.3
    )

    test_alpha = 1
    test_effects = [] if build_in else [
        pe.withStroke(linewidth=22, foreground='white', offset=(-0.5, 1))
    ]

    ax.vlines(
        x[f'ω1_{test_group}'],
        0,
        _visible_proportion(x[f'ω1_{test_group}_P'], ax),
        linewidth=20,
        color='salmon',
        alpha=test_alpha,
        zorder=zorder_1test,
        path_effects=test_effects,
    )
    _plot_tiny_proportion_marker(
        ax, x[f'ω1_{test_group}'], x[f'ω1_{test_group}_P'], 'salmon', alpha=1
    )
    ax.vlines(
        x[f'ω2_{test_group}'],
        0,
        _visible_proportion(x[f'ω2_{test_group}_P'], ax),
        linewidth=20,
        color='steelblue',
        alpha=test_alpha,
        zorder=zorder_2test,
        path_effects=test_effects,
    )
    _plot_tiny_proportion_marker(
        ax, x[f'ω2_{test_group}'], x[f'ω2_{test_group}_P'], 'steelblue', alpha=1
    )
    ax.vlines(
        x[f'ω3_{test_group}'],
        0,
        _visible_proportion(x[f'ω3_{test_group}_P'], ax),
        linewidth=20,
        color='goldenrod',
        alpha=test_alpha,
        zorder=zorder_3test,
        path_effects=test_effects,
    )
    _plot_tiny_proportion_marker(
        ax, x[f'ω3_{test_group}'], x[f'ω3_{test_group}_P'], 'goldenrod', alpha=1
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
        ax.xaxis.set_minor_locator(
            ticker.LogLocator(base=10, subs=np.arange(2, 10))
            if show_subticks else ticker.NullLocator()
        )
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
        ax.xaxis.set_minor_locator(
            ticker.SymmetricalLogLocator(
                base=10, linthresh=thresh, subs=np.arange(2, 10)
            ) if show_subticks else ticker.NullLocator()
        )
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
                0.7, 0.85,
                f'$\\it{{p_1}}$={x["test_pval"]:.2e}\n'
                f'$\\it{{p_2}}$={x["background_pval"]}\n'
                f'$\\it{{p_3}}$={x["shared_pval"]:.2e}',
                fontsize=13, ha='left', va='center', transform=ax.transAxes,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'),
            )

    if filename is not None:
        plt.savefig(filename, dpi=300, transparent=transparent)

    print(
        f"ω values for {gene}: "
        f"ω1_{test_group}={x[f'ω1_{test_group}']:.4f}, "
        f"ω2_{test_group}={x[f'ω2_{test_group}']:.4f}, "
        f"ω3_{test_group}={x[f'ω3_{test_group}']:.4f},\n"
        f"ω1_{ref_group}={x[f'ω1_{ref_group}']:.4f}, "
        f"ω2_{ref_group}={x[f'ω2_{ref_group}']:.4f}, "
        f"ω3_{ref_group}={x[f'ω3_{ref_group}']:.4f};\n"
        f"Proportions: ω1_{test_group}_P={x[f'ω1_{test_group}_P']:.4f}, "
        f"ω2_{test_group}_P={x[f'ω2_{test_group}_P']:.4f}, "
        f"ω3_{test_group}_P={x[f'ω3_{test_group}_P']:.4f},\n"
        f"ω1_{ref_group}_P={x[f'ω1_{ref_group}_P']:.4f}, "
        f"ω2_{ref_group}_P={x[f'ω2_{ref_group}_P']:.4f}, "
        f"ω3_{ref_group}_P={x[f'ω3_{ref_group}_P']:.4f}"
    )

    return fig, ax
