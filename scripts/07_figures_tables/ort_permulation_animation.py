"""Click-through slides for the talk: how the odds ratio permulation test judges a gene.

Writes numbered PNGs (frame_01.png, frame_02.png, ...) to click through in PowerPoint, one per slide
(same canvas, so they can also be stacked in one slide and flipped through).

Each frame is one permulated labeling of the orb-weaver phenotype (a real permulation from the
1,000,000 set, not an illustration). The histogram is the log odds ratio of missing for every tested
HOG under that labeling, and the bin that holds each of two genes is highlighted:

  * a gene that is absent in every member of one family (HOG0029234, Linyphiids): because permulated
    labelings cluster on the tree like the real one, this gene is often just as extreme as it is under the
    real labeling;
  * the one candidate (HOG0071799): it is that extreme under the real labeling only.

A running tally counts the shown permulations in which each gene is at least as extreme as under the
real labeling, and the last frames show the real labeling and the p-values from all 1,000,000 permulations.

The shown permulations are drawn at random from the first 100,000, stratified on whether the clade gene is as
extreme as the real labeling. By default the share equals its p-value; --n-extreme sets the count directly
(the talk version shows 8 of 20, more than the p-value of 0.14 implies, to make the point visible; the last slide
still gives the true p-value). Suspect samples (strong + moderate flags) are treated as missing in
their HOG, as in the final volcano plots.

Run from anywhere:  python scripts/07_figures_tables/ort_permulation_animation.py [--out DIR] [--frames N] [--seed S]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path[:0] = [os.path.join(REPO, "src"), os.path.join(REPO, "scripts", "04_permulation_loss_dup")]
import odds_ratio_test as ort  # noqa: E402
from odds_ratio_test_helpers import (  # noqa: E402
    calculate_odds_batch, calculate_odds_batch_weighted, load_flagged_cells, permulated_foreground_matrices,
)

HIST_FILL, HIST_EDGE = "#EAE1F0", "#C4A9D7"
COL_FG, COL_BG = "salmon", "steelblue"
GENES = [  # (HOG, short description, colour)
    ("N5.HOG0029234", "absent in every Linyphiid", "goldenrod"),
    ("N5.HOG0071799", "real candidate", COL_BG),
]
FULL_RUN = os.path.join(REPO, "results/odds_ratio_test/Results_Nov02/Run5_occ_30-88_1000000x_1M_flags_strong_moderate/results_all.csv")
PERMS_CSV = os.path.join(REPO, "data", "perms_tip_values_100k.csv")      # identical to the first 100,000 rows of the 1M file
BINS = np.arange(-8, 8.01, 0.25)


def load_setup():
    flagged = load_flagged_cells(os.path.join(REPO, "data/flag_table_genes.csv"), strengths=("strong", "moderate"))
    odds = ort.OddsRatioResults(
        genecount_csv=os.path.join(REPO, "data/N5.GeneCount.tsv"), hog_node_genes_tsv=os.path.join(REPO, "data/N5.tsv"),
        time="animation", foreground_list_filename=os.path.join(REPO, "data/orbweavers-list.txt"),
        buscos_filename=os.path.join(REPO, "data/buscos.csv"), flagged_cells=flagged)
    incl = (odds.foreground_bool_arr + odds.background_bool_arr).astype(bool).nonzero()[0]
    tips = pd.read_csv(PERMS_CSV)
    tip_values = [dict(zip(tips.columns[1:], row)) for row in tips[tips.columns[1:]].to_numpy()]
    fg, bg = permulated_foreground_matrices(odds.all_species_arr, tip_values, incl)
    occ = odds.total_occ_arr.ravel()
    loss_idx = np.where((occ >= 30) & (occ <= 88))[0]
    return odds, fg, bg, loss_idx


def loss_lor(odds, fg, bg, rows):
    """Log odds ratio of loss for HOG rows (indices into the full HOG list) under each column of fg/bg."""
    out = calculate_odds_batch(fg, bg, odds.loss_bool_mat[rows], busco_arr=odds.loss_busco_arr)
    flagged = np.where(np.isin(rows, odds.flagged_hog_idx))[0]
    if len(flagged):
        out[flagged] = calculate_odds_batch_weighted(fg, bg, odds.loss_bool_mat[rows[flagged]], odds.weight_matrix("loss", rows[flagged]))
    return out


def build_frames(n_frames, seed, n_extreme=None):
    odds, fg, bg, loss_idx = load_setup()
    hog_row = {h: list(odds.hog_list).index(h) for h, _, _ in GENES}
    gene_rows = np.array([hog_row[h] for h, _, _ in GENES])
    true = np.array([float(odds.loss_lor_arr.ravel()[r]) for r in gene_rows])
    # the two genes under each of the 100,000 permulations, to pick the frames
    gene_perm = np.hstack([loss_lor(odds, fg[:, i:i + 5000], bg[:, i:i + 5000], gene_rows) for i in range(0, fg.shape[1], 5000)])
    extreme0 = gene_perm[0] <= true[0] + ort.TIE_TOLERANCE                      # clade gene as extreme as the real labeling
    p_full = pd.read_csv(FULL_RUN, index_col="HOG").loc[[h for h, _, _ in GENES], "P-value loss more likely in bg"].to_numpy()
    rng = np.random.default_rng(seed)
    k_ext = int(round(p_full[0] * n_frames)) if n_extreme is None else n_extreme     # how many shown permulations have the clade gene as extreme as the real labeling
    pick = np.concatenate([rng.choice(np.where(extreme0)[0], k_ext, replace=False),
                           rng.choice(np.where(~extreme0)[0], n_frames - k_ext, replace=False)])
    rng.shuffle(pick)
    hists = []
    for i in pick:                                                                # whole distribution of the tested HOGs under this labeling
        lor = loss_lor(odds, fg[:, i:i + 1], bg[:, i:i + 1], loss_idx).ravel()
        hists.append(np.histogram(lor, bins=BINS)[0])
    true_all = odds.loss_lor_arr.ravel()[loss_idx]
    frames = {"perm": pick, "hist": np.array(hists), "gene": gene_perm[:, pick].T,       # (frames, genes)
              "true_hist": np.histogram(true_all, bins=BINS)[0], "true_gene": true, "p_full": p_full}
    return frames


LABEL_HALF_W = 1.45         # half the width of a gene label, in log-odds-ratio units (checked against the real text extents by check_frame)


def draw(ax_h, ax_t, F, k, n_frames):
    """Frame k < n_frames: a permulation; k == n_frames: the real labeling; k > n_frames: the summary."""
    ax_h.clear(); ax_t.clear(); ax_t.axis("off")
    real = k >= n_frames
    hist = F["true_hist"] if real else F["hist"][k]
    vals = F["true_gene"] if real else F["gene"][k]
    centers = BINS[:-1] + np.diff(BINS) / 2
    w = np.diff(BINS)[0]
    ax_h.bar(centers, hist, width=w, color=HIST_FILL, edgecolor=HIST_EDGE, linewidth=0.8)
    ymax = F["ymax"]
    shown = min(k + 1, n_frames)                                                   # permulations counted so far
    bins_of = [min(max(np.searchsorted(BINS, v, side="right") - 1, 0), len(centers) - 1) for v in vals]
    # label positions: above the highest bar under the label, over its own bar; if the two labels would collide, push them apart sideways (same height)
    lab_x = [float(np.clip(centers[b], -8 + LABEL_HALF_W + 0.2, 8 - LABEL_HALF_W - 0.2)) for b in bins_of]
    gap = 2 * LABEL_HALF_W + 0.4
    if abs(lab_x[0] - lab_x[1]) < gap:
        mid = (lab_x[0] + lab_x[1]) / 2
        left, right = (0, 1) if lab_x[0] <= lab_x[1] else (1, 0)
        lab_x[left], lab_x[right] = mid - gap / 2, mid + gap / 2
        lo, hi = -8 + LABEL_HALF_W + 0.2, 8 - LABEL_HALF_W - 0.2
        shift = max(lo - lab_x[left], 0) - max(lab_x[right] - hi, 0)                 # keep the pair inside the plot, same distance apart
        lab_x[left] += shift; lab_x[right] += shift
    lab_y = []
    for x_l in lab_x:
        under = hist[(centers > x_l - LABEL_HALF_W - w) & (centers < x_l + LABEL_HALF_W + w)]
        lab_y.append(max(ymax * 0.22, under.max() + ymax * 0.12))
    for j, ((hog, desc, col), v, b) in enumerate(zip(GENES, vals, bins_of)):
        narrow = j == 1 and bins_of[0] == bins_of[1]
        ax_h.bar(centers[b], max(hist[b], ymax * 0.015), width=w * (0.45 if narrow else 1.0), color=col, edgecolor="black", linewidth=0.8, zorder=3)
        if not real:                                                                # under the real labeling the bar is at the dotted line
            ax_h.plot([F["true_gene"][j]] * 2, [0, ymax * 0.9], color=col, ls=":", lw=1.8, zorder=1)
        ax_h.annotate(hog.replace("N5.", ""), xy=(centers[b], max(hist[b], ymax * 0.015) + ymax * 0.012), xytext=(lab_x[j], lab_y[j]), ha="center", va="bottom",
                      fontsize=17, fontweight="bold", color=col, arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8, shrinkA=3, shrinkB=0), zorder=5,
                      bbox=dict(fc="white", ec="none", pad=2.0))
    ax_h.set_xlim(BINS[0], BINS[-1]); ax_h.set_ylim(0, ymax)
    ax_h.set_xlabel("Log odds ratio of missing", fontsize=18, fontweight="bold")
    ax_h.set_ylabel("HOGs", fontsize=18, fontweight="bold")
    ax_h.tick_params(labelsize=16)
    for s in ("top", "right"): ax_h.spines[s].set_visible(False)
    # captions below the x axis, at the two ends (the same row as the axis label, which is centred)
    ax_h.text(0.0, -0.115, "missing more likely in\nnon-orb-weavers", transform=ax_h.transAxes, ha="left", va="top", fontsize=15, color=COL_BG)
    ax_h.text(1.0, -0.115, "missing more likely in\norb-weavers", transform=ax_h.transAxes, ha="right", va="top", fontsize=15, color=COL_FG)
    if not real:
        ax_h.text(0.5, 1.012, "dotted line: the gene under the real labeling", transform=ax_h.transAxes, ha="center", va="bottom", fontsize=15, color="0.35")
    title = "All 1,000,000 permulations" if k > n_frames else "Real labeling" if real else f"Permulation {k + 1} of {n_frames}"
    ax_h.set_title(title, fontsize=22, fontweight="bold", color=("black" if real else "0.25"), pad=40)       # room for the dotted-line note under the title
    # tallies under the plot: gene name on the left, its description above and the result below, with room between the lines
    for j, (hog, desc, col) in enumerate(GENES):
        y = 0.74 - 0.52 * j
        if k <= n_frames - 1:
            n_ext = int((F["gene"][:shown, j] <= F["true_gene"][j] + ort.TIE_TOLERANCE).sum())
            txt = f"as extreme as under the real labeling in {n_ext} of {shown} permulation{'s' if shown != 1 else ''}"
        elif k == n_frames:
            txt = f"real log odds ratio {F['true_gene'][j]:.1f}"
        else:
            p = F["p_full"][j]
            n_as = int(round(p * 1_000_001)) - 1                                    # p = (count + 1) / (permulations + 1)
            txt = f"as extreme in {n_as:,} of 1,000,000 permulations (p = {p:.2g})"
        ax_t.text(0.0, y, f"{hog.replace('N5.', '')}", color=col, fontsize=19, fontweight="bold", va="center", transform=ax_t.transAxes)
        ax_t.text(0.215, y + 0.15, desc, color="0.3", fontsize=15, va="center", transform=ax_t.transAxes)
        ax_t.text(0.215, y - 0.10, txt, color="black", fontsize=17, va="center", transform=ax_t.transAxes)


def _boxes_overlap(a, b, pad=0):
    return not (a.x1 + pad <= b.x0 or b.x1 + pad <= a.x0 or a.y1 + pad <= b.y0 or b.y1 + pad <= a.y0)


def _seg_hits_box(p, q, box, n=60):
    return any(box.x0 <= p[0] + (q[0] - p[0]) * t <= box.x1 and box.y0 <= p[1] + (q[1] - p[1]) * t <= box.y1 for t in np.linspace(0, 1, n))


def check_frame(fig, ax_h, ax_t, F, k, n_frames):
    """Problems in the drawn frame: any text touching another text, an arrow running through another label, a label sitting on a bar."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    anns = [t for t in ax_h.texts if hasattr(t, "arrow_patch") and t.arrow_patch is not None]
    texts = [(t.get_text().replace("\n", " ")[:30], matplotlib.text.Text.get_window_extent(t, r)) for t in list(ax_h.texts) + list(ax_t.texts) + [ax_h.title, ax_h.xaxis.label]]
    probs = []
    pad = 3                                                                       # pixels of clear space required between texts
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            if _boxes_overlap(texts[i][1], texts[j][1], pad):
                probs.append(f"text '{texts[i][0]}' touches '{texts[j][0]}'")
    boxes = {id(t): matplotlib.text.Text.get_window_extent(t, r) for t in anns}
    segs = {}
    for t in anns:
        bx = boxes[id(t)]
        tip = ax_h.transData.transform(t.xy)
        segs[id(t)] = ((bx.x0 + bx.x1) / 2, bx.y0), tuple(tip)
    for t in anns:
        for name, bx in texts:
            if bx.x0 == boxes[id(t)].x0 and bx.y0 == boxes[id(t)].y0:
                continue
            if _seg_hits_box(*segs[id(t)], bx):
                probs.append(f"arrow of '{t.get_text()}' runs through '{name}'")
    if len(anns) == 2:
        (a0, a1), (b0, b1) = segs[id(anns[0])], segs[id(anns[1])]
        if abs(a0[0] - b0[0]) < 12 and abs(a1[0] - b1[0]) < 12:
            probs.append("the two arrows run on top of each other")
    hist = F["true_hist"] if k >= n_frames else F["hist"][k]
    centers = BINS[:-1] + np.diff(BINS) / 2
    for t in anns:
        bx = boxes[id(t)]
        (x0, y0), (x1, _) = ax_h.transData.inverted().transform([(bx.x0, bx.y0), (bx.x1, bx.y0)])
        under = hist[(centers + 0.125 > x0) & (centers - 0.125 < x1)]
        if len(under) and y0 < under.max() + 1:
            probs.append(f"label '{t.get_text()}' sits on a bar")
        if bx.x0 < ax_h.bbox.x0 or bx.x1 > ax_h.bbox.x1 or bx.y1 > ax_h.bbox.y1:
            probs.append(f"label '{t.get_text()}' leaves the plot area")
    return probs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/Users/calvin/Library/CloudStorage/OneDrive-JohnsHopkins/0 gordus lab/0 Selection paper final/talks/ort_final_nov17")
    ap.add_argument("--frames", type=int, default=20, help="permulations shown (plus the real labeling and the summary)")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--n-extreme", type=int, default=None, help="shown permulations in which the clade gene is as extreme as the real labeling (default: its p-value x frames). The talk uses 8 of 20, more than the p-value of 0.14 implies, to make the point visible; the last slide still gives the true p-value.")
    ap.add_argument("--folder", default="permulation_clickthrough", help="subfolder of --out for the PNGs")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    F = build_frames(args.frames, args.seed, args.n_extreme)
    F["ymax"] = float(max(F["hist"].max(), F["true_hist"].max()) * 1.28)       # headroom for the labels above the tallest bars
    n = args.frames
    order = list(range(n)) + [n, n + 1]                                            # permulations, the real labeling, the summary
    fig = plt.figure(figsize=(12, 8.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.1, 1.45], hspace=0.30)
    ax_h, ax_t = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    fig.subplots_adjust(left=0.08, right=0.98, top=0.87, bottom=0.07)
    out = os.path.join(args.out, args.folder)
    os.makedirs(out, exist_ok=True)
    for i, k in enumerate(order, 1):
        draw(ax_h, ax_t, F, k, n)
        for pr in check_frame(fig, ax_h, ax_t, F, k, n):
            print(f"frame {i:02d}: {pr}")
        fig.savefig(os.path.join(out, f"frame_{i:02d}.png"), dpi=200, facecolor="white")
    print(f"{len(order)} frames written to {out}")


if __name__ == "__main__":
    main()
