"""Shared helper utilities for the odds ratio permulation workflow."""

import os
import sys
from datetime import datetime
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

import id_converter

_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_SRC_DIR, os.pardir, os.pardir))
CONSOLE_PRINT_WIDTH = 96


def _fixed_width_lines(text: str, width: int = CONSOLE_PRINT_WIDTH) -> List[str]:
    """Return fixed-width wrapped lines for aligned console printing."""

    def _wrap_no_word_breaks(line: str) -> List[str]:
        if len(line) <= width:
            return [line]

        wrapped: List[str] = []
        remaining = line
        break_chars = (" ", "\t", "/", "_", "-", ",", ";")

        while len(remaining) > width:
            window = remaining[: width + 1]
            break_at = max(window.rfind(ch) for ch in break_chars)

            if break_at <= 0:
                wrapped.append(remaining[:width])
                remaining = remaining[width:]
                continue

            split_idx = break_at + 1
            wrapped.append(remaining[:split_idx].rstrip())

            if remaining[break_at] in {" ", "\t"}:
                remaining = remaining[split_idx:].lstrip()
            else:
                remaining = remaining[split_idx:]

        wrapped.append(remaining)
        return wrapped

    wrapped_lines: List[str] = []
    for raw_line in str(text).splitlines() or [""]:
        segments = _wrap_no_word_breaks(raw_line)
        if not segments:
            wrapped_lines.append("")
            continue
        wrapped_lines.extend(segments)
    return wrapped_lines


def _cprint(text: str = "", width: int = CONSOLE_PRINT_WIDTH) -> None:
    """Print fixed-width lines to stdout."""
    for line in _fixed_width_lines(text, width=width):
        print(line)


def _emit(text: str, file_obj, width: int = CONSOLE_PRINT_WIDTH) -> None:
    """Write output to file object; fixed-width only when writing to stdout."""
    if file_obj is sys.stdout:
        _cprint(text, width=width)
    else:
        print(text, file=file_obj)


def _resolve_repo_path(path: Optional[str]) -> Optional[str]:
    """Resolve user paths with repo-relative defaults and external path support.

    Resolution order:
    1. Absolute paths and '~' expanded paths are used as-is.
    2. Explicit './' or '../' paths are resolved from current working directory.
    3. Bare relative paths are resolved from repository root.
    """
    if path is None:
        return None

    expanded = os.path.expanduser(path)
    if os.path.isabs(expanded):
        return expanded

    if (
        expanded.startswith(f".{os.sep}")
        or expanded == "."
        or expanded.startswith(f"..{os.sep}")
        or expanded == ".."
    ):
        return os.path.abspath(expanded)

    return os.path.abspath(os.path.join(_REPO_ROOT, expanded))


def _unique_results_dir(
    parent_dir: str,
    time_obj: datetime,
    min_occ: int,
    max_occ: Optional[int],
    permulation_reps: int,
    dir_suffix: Optional[str] = None,
) -> str:
    """Create a dated parent folder and return a unique run subdirectory.
    
    Args:
        parent_dir: Parent directory path
        time_obj: DateTime object for directory dating
        min_occ: Minimum occupancy threshold
        max_occ: Maximum occupancy threshold
        permulation_reps: Number of permulation repetitions
        dir_suffix: Optional string to append to the end of the directory name
    """
    date_short = time_obj.strftime("%b%d")
    dated_parent = f"{parent_dir}/Results_{date_short}"
    os.makedirs(dated_parent, exist_ok=True)

    occ_range = f"{min_occ}-{max_occ}"
    base_name = f"occ_{occ_range}_{permulation_reps}x"

    existing_run_nums = []
    if os.path.exists(dated_parent):
        for entry in os.listdir(dated_parent):
            if entry.startswith("Run"):
                try:
                    run_num = int(entry.split("_")[0][3:])
                    existing_run_nums.append(run_num)
                except (ValueError, IndexError):
                    pass

    next_run_num = max(existing_run_nums) + 1 if existing_run_nums else 1
    subdir_name = f"Run{next_run_num}_{base_name}"
    if dir_suffix:
        subdir_name = f"{subdir_name}_{dir_suffix}"
    return f"{dated_parent}/{subdir_name}"


def drop_empty_cols(df, print_txt=True):
    """Drop columns where all entries are 0."""
    num_columns_before = df.shape[1]
    if print_txt:
        _cprint(
            f"Number of columns before dropping empty columns: {num_columns_before}"
        )

    df_cleaned = df.loc[:, (df.ne(0)).any(axis=0)]
    num_columns_after = df_cleaned.shape[1]

    if print_txt:
        _cprint(f"Number of columns after dropping empty columns: {num_columns_after}")
        _cprint("Species with no sequences in any orthogroup have been dropped.")

    return df_cleaned


def occupancy_filter(arr, min_occ, max_occ, total_occ_arr):
    """Filter an array of values by occupancy threshold."""
    if max_occ is None:
        max_occ = total_occ_arr.max()

    idx = np.asarray((total_occ_arr >= min_occ) & (total_occ_arr <= max_occ)).nonzero()[0]
    return arr[idx]


def filter_for_sp_of_interest(df, genecount_df, species_name):
    """Filter the DataFrame for HOGs that include a species of interest."""
    sp_of_interest_present = genecount_df[genecount_df[species_name] != 0]
    sp_of_interest_present_hogs = set(sp_of_interest_present.index.values)
    df_fltrd_sp_of_int = df[df.index.isin(sp_of_interest_present_hogs)]
    return len(df_fltrd_sp_of_int)


def calculate_odds(foreground_bool_arr, background_bool_arr, test_bool_mat, busco_arr):
    """Calculate the odds ratio and log odds ratio."""
    foreground_bool_arr = foreground_bool_arr.reshape(foreground_bool_arr.size, 1)
    background_bool_arr = background_bool_arr.reshape(background_bool_arr.size, 1)

    if busco_arr is not None:
        foreground_bool_arr = foreground_bool_arr.flatten() * busco_arr
        background_bool_arr = background_bool_arr.flatten() * busco_arr

        foreground_bool_arr = foreground_bool_arr.reshape(foreground_bool_arr.size, 1)
        background_bool_arr = background_bool_arr.reshape(background_bool_arr.size, 1)

    foreground_yes_arr = np.matmul(test_bool_mat, foreground_bool_arr)
    background_yes_arr = np.matmul(test_bool_mat, background_bool_arr)
    test_inv_bool_mat = 1 - test_bool_mat
    foreground_no_arr = np.matmul(test_inv_bool_mat, foreground_bool_arr)
    background_no_arr = np.matmul(test_inv_bool_mat, background_bool_arr)

    foreground_yes_arr += 0.5
    background_yes_arr += 0.5
    foreground_no_arr += 0.5
    background_no_arr += 0.5

    odds_foreground_arr = foreground_yes_arr / foreground_no_arr
    odds_background_arr = background_yes_arr / background_no_arr
    odds_ratio_arr = odds_foreground_arr / odds_background_arr
    return np.log(odds_ratio_arr)


def load_permulation_tip_values_from_csv(csv_path: str) -> List[Dict[str, float]]:
    """Load permulated tip values from a CSV file."""
    csv_path = _resolve_repo_path(csv_path)
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Permulation tip-values CSV not found: {csv_path}")

    tip_df = pd.read_csv(csv_path)
    if tip_df.empty:
        raise ValueError(f"Permulation tip-values CSV is empty: {csv_path}")

    species_cols = [c for c in tip_df.columns if c != "perm_id"]
    if not species_cols:
        raise ValueError(
            f"Permulation tip-values CSV has no species columns: {csv_path}"
        )

    tip_numeric = tip_df[species_cols].apply(pd.to_numeric, errors="raise")
    if tip_numeric.isnull().any().any():
        raise ValueError(f"Permulation tip-values CSV contains NaN values: {csv_path}")

    return [
        {species: float(value) for species, value in row.items()}
        for row in tip_numeric.to_dict(orient="records")
    ]


def save_loc_list(df, NX_path, output_path):
    """Convert HOGs to LOCs and save them one per line."""
    locs_df = id_converter.convert_hogs_to_locs(df, NX_path, show_progress=False)
    loc_ids = locs_df["LOC"].dropna().unique()
    with open(output_path, "w", encoding="utf-8") as file_obj:
        for loc in loc_ids:
            file_obj.write(f"{loc}\n")
    print(f"Wrote {len(loc_ids)} items to {output_path}")
    return loc_ids


__all__ = [
    "CONSOLE_PRINT_WIDTH",
    "_cprint",
    "_emit",
    "_fixed_width_lines",
    "_resolve_repo_path",
    "_unique_results_dir",
    "calculate_odds",
    "drop_empty_cols",
    "filter_for_sp_of_interest",
    "load_permulation_tip_values_from_csv",
    "occupancy_filter",
    "save_loc_list",
]
