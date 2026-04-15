"""Shared helpers for split-axis plotting layouts."""

import numpy as np
from matplotlib import ticker


def draw_y_axis_break_marks(ax_top, ax_bottom, dx=0.012, angle_deg=55.0):
    """Draw equal-angle diagonal break marks for stacked y-split axes."""

    def _dy_for_angle(ax, dx_local, angle_local):
        bbox = ax.get_position()
        width = max(bbox.width, 1e-9)
        height = max(bbox.height, 1e-9)
        return dx_local * (width / height) * np.tan(np.deg2rad(angle_local))

    dy_top = _dy_for_angle(ax_top, dx, angle_deg)
    dy_bottom = _dy_for_angle(ax_bottom, dx, angle_deg)
    kwargs = dict(color="k", clip_on=False, linewidth=1)

    ax_top.plot((-dx, +dx), (-dy_top, +dy_top), transform=ax_top.transAxes, **kwargs)
    ax_top.plot((1 - dx, 1 + dx), (-dy_top, +dy_top), transform=ax_top.transAxes, **kwargs)
    ax_bottom.plot(
        (-dx, +dx),
        (1 - dy_bottom, 1 + dy_bottom),
        transform=ax_bottom.transAxes,
        **kwargs,
    )
    ax_bottom.plot(
        (1 - dx, 1 + dx),
        (1 - dy_bottom, 1 + dy_bottom),
        transform=ax_bottom.transAxes,
        **kwargs,
    )


def compute_shared_split_tick_step(lower_ylim, upper_ylim, tick_count=4):
    """Compute a shared y-tick step for split axes."""
    intervals = []
    for ymin, ymax in (lower_ylim, upper_ylim):
        auto_locator = ticker.AutoLocator()
        ticks = auto_locator.tick_values(ymin, ymax)
        ticks = ticks[(ticks >= ymin) & (ticks <= ymax)]
        if ticks.size >= 2:
            diffs = np.diff(ticks)
            diffs = diffs[diffs > 0]
            if diffs.size:
                intervals.append(float(np.min(diffs)))

    if intervals:
        return max(intervals)

    span = max(lower_ylim[1] - lower_ylim[0], upper_ylim[1] - upper_ylim[0])
    if span <= 0:
        return 1.0
    return span / max(tick_count, 1)


def split_axis_ticks(lower_ylim, upper_ylim, tick_count=4):
    """Return split-axis y-ticks with a shared step and no ticks at the break."""
    lower_min, lower_max = lower_ylim
    upper_min, upper_max = upper_ylim
    step = compute_shared_split_tick_step(lower_ylim, upper_ylim, tick_count=tick_count)

    lower_start = np.ceil(lower_min / step) * step
    lower_ticks = np.arange(lower_start, lower_max + step * 0.5, step)

    upper_start = np.ceil(upper_min / step) * step
    upper_ticks = np.arange(upper_start, upper_max + step * 0.5, step)

    eps = step * 1e-6
    lower_ticks = lower_ticks[(lower_ticks >= lower_min - eps) & (lower_ticks < lower_max - eps)]
    upper_ticks = upper_ticks[(upper_ticks > upper_min + eps) & (upper_ticks <= upper_max + eps)]

    return lower_ticks, upper_ticks


def resolve_split_height_ratios(lower_ylim, upper_ylim, split_y_height_ratios):
    """Resolve split panel height ratios; 'auto' uses y-span proportions."""
    if split_y_height_ratios == "auto" or split_y_height_ratios is None:
        lower_span = float(lower_ylim[1] - lower_ylim[0])
        upper_span = float(upper_ylim[1] - upper_ylim[0])
        if lower_span <= 0 or upper_span <= 0:
            return (1, 1)
        return (upper_span, lower_span)

    if len(split_y_height_ratios) != 2:
        raise ValueError("split_y_height_ratios must contain exactly two values (top, bottom).")

    return split_y_height_ratios


__all__ = [
    "compute_shared_split_tick_step",
    "draw_y_axis_break_marks",
    "resolve_split_height_ratios",
    "split_axis_ticks",
]
