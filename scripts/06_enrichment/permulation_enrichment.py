"""Permulation-based GO enrichment of gene loss / duplication log odds ratios.

topGO's Fisher test on a list of significant genes treats genes as independent draws. Here every
HOG is scored against the same species labeling, so a single atypical labeling shifts many HOGs
(and many GO terms) the same way. This module instead scores each GO term directly and builds the
null from the permulated labelings, which keeps that dependence.

For each test (loss, duplication) and each labeling (the true one and every permulation):
  1. each HOG gets a signed score z = (LOR - mean_perm) / sd_perm, where mean_perm and sd_perm are
     that HOG's log odds ratio mean and standard deviation across the permulations
     (positive = event more likely in the foreground, i.e. orb-weavers);
  2. HOGs map to U. diversus genes through their best BLAST hit (a gene's score is the mean over
     its HOGs), and a GO term's score is the mean score of its genes (annotations are propagated up
     the GO hierarchy, as topGO does).

The term p-values compare the true labeling's score with the permulated scores (add-one
estimator, ties counted as extreme). Both directions are reported: "fg" = higher in the
foreground, "bg" = higher in the background. Corrections are applied within each ontology:
Benjamini-Hochberg over both directions stacked, and a max-T family-wise adjusted p-value that uses
the same permulations (so it accounts for overlap between terms exactly).

Usage (from the repo root):
    python scripts/06_enrichment/permulation_enrichment.py \
        --permulations-csv data/perms_tip_values_100k.csv --out-dir results/permulation_enrichment/run1
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
from tqdm.auto import tqdm

_HERE = os.path.dirname(os.path.abspath(__file__))
_STAGE04 = os.path.join(_HERE, os.pardir, "04_permulation_loss_dup")
_SRC = os.path.join(_HERE, os.pardir, os.pardir, "src")
sys.path[:0] = [_STAGE04, _SRC]

import odds_ratio_test as ort  # noqa: E402
from odds_ratio_test_helpers import (  # noqa: E402
    _resolve_repo_path,
    benjamini_hochberg,
    calculate_odds_batch,
    load_permulation_matrices,
)

TIE_TOLERANCE = ort.TIE_TOLERANCE
ONTOLOGIES = ("BP", "MF", "CC")


# ----------------------------------------------------------------------------- inputs


def load_hog_gene_map(blast_tsv):
    """HOG -> reference gene id (numeric string, no 'LOC' prefix); drops HOGs without a hit."""
    blast = pd.read_csv(blast_tsv, sep="\t", index_col=0)
    gene = blast["Gene"].astype(str)
    gene = gene[gene.str.startswith("LOC")].str[3:]
    return gene


def load_term_sets(term_file, universe_genes, min_genes, max_genes):
    """GO terms with their genes restricted to the universe, filtered by size."""
    terms = pd.read_csv(term_file, sep="\t", dtype=str)
    universe = set(universe_genes)
    kept = []
    for row in terms.itertuples(index=False):
        genes = [g for g in row.genes.split(",") if g in universe]
        if min_genes <= len(genes) <= max_genes:
            kept.append((row.GOID, row.ontology, row.term, genes))
    return pd.DataFrame(kept, columns=["GOID", "ontology", "term", "genes"])


# ----------------------------------------------------------------------- score matrix


def build_term_hog_matrix(terms, hog_ids, hog_gene):
    """Sparse (terms x HOGs) averaging matrix: term score = mean over its genes of the gene score,
    gene score = mean over its HOGs of the HOG score."""
    genes = pd.Index(sorted({g for gs in terms["genes"] for g in gs}))
    hog_gene_index = genes.get_indexer(hog_gene.reindex(hog_ids).fillna("").to_numpy())
    has_gene = hog_gene_index >= 0

    # gene x HOG, each gene averages its HOGs
    hogs_per_gene = np.bincount(hog_gene_index[has_gene], minlength=len(genes))
    cols = np.where(has_gene)[0]
    rows = hog_gene_index[has_gene]
    gene_hog = sp.csr_matrix((1.0 / hogs_per_gene[rows], (rows, cols)), shape=(len(genes), len(hog_ids)))

    # term x gene, each term averages its genes
    t_rows, t_cols = [], []
    for t, gs in enumerate(terms["genes"]):
        t_rows.extend([t] * len(gs))
        t_cols.extend(genes.get_indexer(gs))
    sizes = terms["genes"].map(len).to_numpy()
    term_gene = sp.csr_matrix(
        (1.0 / sizes[t_rows], (t_rows, t_cols)), shape=(len(terms), len(genes))
    )
    return (term_gene @ gene_hog).tocsr()


def lor_chunks(mat_rows, busco, foreground, background, chunk_size):
    """Yield (start, stop, log odds ratios of shape (n_rows, chunk)) over permulation chunks."""
    n_perm = foreground.shape[1]
    for start in range(0, n_perm, chunk_size):
        stop = min(start + chunk_size, n_perm)
        yield start, stop, calculate_odds_batch(
            foreground[:, start:stop], background[:, start:stop], mat_rows, busco
        )


def hog_null_moments(mat_rows, busco, foreground, background, chunk_size, desc):
    """Mean and standard deviation of every HOG's log odds ratio across the permulations."""
    total = np.zeros(mat_rows.shape[0])
    total_sq = np.zeros(mat_rows.shape[0])
    n_perm = foreground.shape[1]
    n_chunks = -(-n_perm // chunk_size)
    for _, _, lor in tqdm(
        lor_chunks(mat_rows, busco, foreground, background, chunk_size), total=n_chunks, desc=desc
    ):
        total += lor.sum(axis=1)
        total_sq += (lor**2).sum(axis=1)
    mean = total / n_perm
    sd = np.sqrt(np.maximum(total_sq / n_perm - mean**2, 0.0))
    return mean, sd


# --------------------------------------------------------------------------- the test


def run_enrichment(
    test,
    true_odds,
    hog_gene,
    term_file,
    foreground,
    background,
    min_occ,
    max_occ,
    min_genes,
    max_genes,
    chunk_size=500,
):
    """Permulation enrichment for one test ("loss" or "dup"); returns a DataFrame of GO terms."""
    n_perm = foreground.shape[1]
    occ = true_odds.total_occ_arr.ravel()
    in_range = (occ >= min_occ) & ((occ <= max_occ) if test == "loss" else True)
    universe_rows = np.where(in_range)[0]
    universe_hogs = pd.Index(true_odds.hog_list[universe_rows])

    terms = load_term_sets(
        term_file, hog_gene.reindex(universe_hogs).dropna().unique(), min_genes, max_genes
    )
    a_matrix = build_term_hog_matrix(terms, universe_hogs, hog_gene)

    mat = getattr(true_odds, f"{test}_bool_mat")
    busco = getattr(true_odds, f"{test}_busco_arr", None)
    observed_lor = getattr(true_odds, f"{test}_lor_arr").ravel()

    # Only HOGs that feed at least one term need scoring
    used = np.unique(a_matrix.nonzero()[1])
    a_used = a_matrix[:, used].tocsr()
    mat_used = mat[universe_rows[used]]
    obs_used = observed_lor[universe_rows[used]]
    print(
        f"[{test}] {len(universe_hogs)} HOGs in range, {len(used)} map to a gene with GO terms, "
        f"{len(terms)} terms (BP/MF/CC: "
        f"{', '.join(str(int((terms.ontology == o).sum())) for o in ONTOLOGIES)})"
    )

    # Pass 1: per-HOG null mean / sd across permulations
    mean, sd = hog_null_moments(mat_used, busco, foreground, background, chunk_size, f"{test}: HOG moments")
    usable = sd > 1e-12
    scale = np.where(usable, sd, np.inf)
    z_obs = np.where(usable, (obs_used - mean) / scale, 0.0)
    score_obs = a_used @ z_obs

    # Pass 2: term score null moments and exceedance counts
    n_terms = len(terms)
    t_sum, t_sq = np.zeros(n_terms), np.zeros(n_terms)
    count_ge, count_le = np.zeros(n_terms), np.zeros(n_terms)
    n_chunks = -(-n_perm // chunk_size)
    for _, _, lor in tqdm(
        lor_chunks(mat_used, busco, foreground, background, chunk_size),
        total=n_chunks,
        desc=f"{test}: term scores",
    ):
        scores = a_used @ ((lor - mean[:, None]) / scale[:, None])
        t_sum += scores.sum(axis=1)
        t_sq += (scores**2).sum(axis=1)
        count_ge += (scores >= score_obs[:, None] - TIE_TOLERANCE).sum(axis=1)
        count_le += (scores <= score_obs[:, None] + TIE_TOLERANCE).sum(axis=1)
    null_mean = t_sum / n_perm
    null_sd = np.sqrt(np.maximum(t_sq / n_perm - null_mean**2, 0.0))
    term_ok = null_sd > 1e-12
    safe_sd = np.where(term_ok, null_sd, np.inf)
    std_obs = (score_obs - null_mean) / safe_sd

    # Pass 3: max-T null per ontology (largest / smallest standardized term score per permulation)
    max_null = {o: np.full(n_perm, -np.inf) for o in ONTOLOGIES}
    min_null = {o: np.full(n_perm, np.inf) for o in ONTOLOGIES}
    ont_rows = {o: np.where((terms.ontology == o).to_numpy() & term_ok)[0] for o in ONTOLOGIES}
    for start, stop, lor in tqdm(
        lor_chunks(mat_used, busco, foreground, background, chunk_size),
        total=n_chunks,
        desc=f"{test}: max-T null",
    ):
        std_scores = (a_used @ ((lor - mean[:, None]) / scale[:, None]) - null_mean[:, None]) / safe_sd[:, None]
        for o, rows in ont_rows.items():
            if len(rows):
                max_null[o][start:stop] = std_scores[rows].max(axis=0)
                min_null[o][start:stop] = std_scores[rows].min(axis=0)

    result = terms[["GOID", "ontology", "term"]].copy()
    result["n_genes"] = terms["genes"].map(len)
    result["observed_score"] = score_obs
    result["null_mean"] = null_mean
    result["null_sd"] = null_sd
    result["standardized_score"] = std_obs
    result["p_fg"] = (count_ge + 1) / (n_perm + 1)
    result["p_bg"] = (count_le + 1) / (n_perm + 1)
    result.loc[~term_ok, ["p_fg", "p_bg"]] = np.nan

    # Corrections within each ontology: BH with both directions stacked; max-T adjusted p-values
    for col in ("q_fg", "q_bg", "p_maxT_fg", "p_maxT_bg"):
        result[col] = np.nan
    for o, rows in ont_rows.items():
        if not len(rows):
            continue
        stacked = np.concatenate([result["p_fg"].to_numpy()[rows], result["p_bg"].to_numpy()[rows]])
        q = benjamini_hochberg(stacked)
        result.loc[result.index[rows], "q_fg"] = q[: len(rows)]
        result.loc[result.index[rows], "q_bg"] = q[len(rows):]
        for col, null_vals, sign in (("p_maxT_fg", max_null[o], 1), ("p_maxT_bg", min_null[o], -1)):
            sorted_null = np.sort(sign * null_vals)
            obs = sign * std_obs[rows]
            n_ge = n_perm - np.searchsorted(sorted_null, obs - TIE_TOLERANCE, side="left")
            result.loc[result.index[rows], col] = (n_ge + 1) / (n_perm + 1)
    return result.sort_values("p_fg").reset_index(drop=True)


# ----------------------------------------------------------------- calibration check


def null_on_null(
    test,
    true_odds,
    hog_gene,
    term_file,
    foreground,
    background,
    min_occ,
    max_occ,
    min_genes,
    max_genes,
    n_null,
    n_pseudo,
    chunk_size=500,
):
    """Treat n_pseudo permulations as if they were the observed labeling and test them against
    n_null other permulations. With a valid null the term p-values are uniform and the
    family-wise adjusted p-values reach 0.05 for about 5% of the pseudo-labelings."""
    null_fg, null_bg = foreground[:, :n_null], background[:, :n_null]
    pseudo_fg = foreground[:, n_null : n_null + n_pseudo]
    pseudo_bg = background[:, n_null : n_null + n_pseudo]

    occ = true_odds.total_occ_arr.ravel()
    in_range = (occ >= min_occ) & ((occ <= max_occ) if test == "loss" else True)
    universe_rows = np.where(in_range)[0]
    universe_hogs = pd.Index(true_odds.hog_list[universe_rows])
    terms = load_term_sets(
        term_file, hog_gene.reindex(universe_hogs).dropna().unique(), min_genes, max_genes
    )
    a_matrix = build_term_hog_matrix(terms, universe_hogs, hog_gene)
    used = np.unique(a_matrix.nonzero()[1])
    a_used = a_matrix[:, used].tocsr()
    mat_used = getattr(true_odds, f"{test}_bool_mat")[universe_rows[used]]
    busco = getattr(true_odds, f"{test}_busco_arr", None)

    mean, sd = hog_null_moments(mat_used, busco, null_fg, null_bg, chunk_size, f"{test}: null moments")
    scale = np.where(sd > 1e-12, sd, np.inf)

    def term_scores(fg, bg):
        out = []
        for _, _, lor in lor_chunks(mat_used, busco, fg, bg, chunk_size):
            out.append(a_used @ ((lor - mean[:, None]) / scale[:, None]))
        return np.concatenate(out, axis=1)

    null_scores = term_scores(null_fg, null_bg)  # (terms, n_null)
    pseudo_scores = term_scores(pseudo_fg, pseudo_bg)  # (terms, n_pseudo)

    null_mean = null_scores.mean(axis=1)
    null_sd = null_scores.std(axis=1)
    ok = null_sd > 1e-12
    sorted_null = np.sort(null_scores, axis=1)

    p_fg = np.ones((len(terms), n_pseudo))
    for t in np.where(ok)[0]:
        n_ge = n_null - np.searchsorted(sorted_null[t], pseudo_scores[t] - TIE_TOLERANCE, side="left")
        p_fg[t] = (n_ge + 1) / (n_null + 1)
    p_fg = p_fg[ok]

    # family-wise: max standardized score over terms (within ontology) vs its null
    fwer_hits = {}
    for o in ONTOLOGIES:
        rows = np.where((terms.ontology == o).to_numpy() & ok)[0]
        if not len(rows):
            continue
        std_null = (null_scores[rows] - null_mean[rows, None]) / null_sd[rows, None]
        std_pseudo = (pseudo_scores[rows] - null_mean[rows, None]) / null_sd[rows, None]
        sorted_max = np.sort(std_null.max(axis=0))
        n_ge = n_null - np.searchsorted(sorted_max, std_pseudo.max(axis=0) - TIE_TOLERANCE, side="left")
        fwer_hits[o] = float((((n_ge + 1) / (n_null + 1)) <= 0.05).mean())

    return {
        "test": test,
        "n_terms": int(ok.sum()),
        "n_null": n_null,
        "n_pseudo": n_pseudo,
        "frac_p_le_0.05": float((p_fg <= 0.05).mean()),
        "frac_p_le_0.01": float((p_fg <= 0.01).mean()),
        "frac_p_le_0.001": float((p_fg <= 0.001).mean()),
        "fwer_fraction_maxT_le_0.05": fwer_hits,
    }


# ------------------------------------------------------------------------------ main


def write_summary(path, results, args, n_perm):
    lines = [
        "Permulation-based GO enrichment",
        f"permulations: {n_perm} ({args.permulations_csv})",
        f"term size filter: {args.min_genes}-{args.max_genes} genes; occupancy {args.min_occ}-{args.max_occ} (loss), >= {args.min_occ} (duplication)",
        "fg = higher log odds ratio in orb-weavers (event more likely in orb-weavers); bg = more likely in non-orb-weavers",
        "",
    ]
    for test, df in results.items():
        lines.append(f"== {test} ==")
        for o in ONTOLOGIES:
            sub = df[df.ontology == o]
            if sub.empty:
                continue
            lines.append(
                f"{o}: {len(sub)} terms | p<=0.025: fg {int((sub.p_fg <= .025).sum())}, bg {int((sub.p_bg <= .025).sum())}"
                f" | BH q<0.05: fg {int((sub.q_fg < .05).sum())}, bg {int((sub.q_bg < .05).sum())}"
                f" | max-T p<0.05: fg {int((sub.p_maxT_fg < .05).sum())}, bg {int((sub.p_maxT_bg < .05).sum())}"
                f" | smallest q: fg {sub.q_fg.min():.3f}, bg {sub.q_bg.min():.3f}"
            )
        lines.append("")
    with open(path, "w", encoding="utf-8") as file_obj:
        file_obj.write("\n".join(lines))
    print("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--permulations-csv", default="data/perms_tip_values_100k.csv")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--tests", nargs="+", default=["loss", "dup"], choices=["loss", "dup"])
    parser.add_argument("--blast-tsv", default="data/N5_occ_30_blasted_udiv.tsv")
    parser.add_argument("--term-file", default="data/udiv_go_term_gene_sets.tsv")
    parser.add_argument("--min-occ", type=int, default=30)
    parser.add_argument("--max-occ", type=int, default=88, help="Upper occupancy bound for loss")
    parser.add_argument("--min-genes", type=int, default=5)
    parser.add_argument("--max-genes", type=int, default=500)
    parser.add_argument("--reps", type=int, default=None, help="Use only the first N permulations")
    parser.add_argument("--chunk-size", type=int, default=500)
    parser.add_argument(
        "--diagnostic",
        nargs=2,
        type=int,
        metavar=("N_NULL", "N_PSEUDO"),
        help="Run the null-on-null calibration check instead of the real test",
    )
    args = parser.parse_args()

    repo = _resolve_repo_path(".")
    true_odds = ort.OddsRatioResults(
        genecount_csv=_resolve_repo_path("data/N5.GeneCount.tsv"),
        hog_node_genes_tsv=_resolve_repo_path("data/N5.tsv"),
        time=time.strftime("%Y-%m-%d at %H:%M:%S"),
        foreground_list_filename=_resolve_repo_path("data/orbweavers-list.txt"),
        background_list_filename=None,
        buscos_filename=_resolve_repo_path("data/buscos.csv"),
    )
    included = (true_odds.foreground_bool_arr + true_odds.background_bool_arr).astype(bool).nonzero()[0]
    foreground, background = load_permulation_matrices(
        _resolve_repo_path(args.permulations_csv), true_odds.all_species_arr, included
    )
    if args.reps is not None:
        foreground, background = foreground[:, : args.reps], background[:, : args.reps]
    hog_gene = load_hog_gene_map(_resolve_repo_path(args.blast_tsv))
    term_file = _resolve_repo_path(args.term_file)
    shared = dict(
        true_odds=true_odds,
        hog_gene=hog_gene,
        term_file=term_file,
        foreground=foreground,
        background=background,
        min_occ=args.min_occ,
        max_occ=args.max_occ,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        chunk_size=args.chunk_size,
    )

    out_dir = _resolve_repo_path(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    if args.diagnostic:
        n_null, n_pseudo = args.diagnostic
        if n_null + n_pseudo > foreground.shape[1]:
            raise ValueError("n_null + n_pseudo exceeds the number of permulations")
        lines = []
        for test in args.tests:
            res = null_on_null(test, n_null=n_null, n_pseudo=n_pseudo, **shared)
            lines.append(str(res))
            print(res)
        with open(os.path.join(out_dir, "null_on_null.txt"), "w", encoding="utf-8") as file_obj:
            file_obj.write("\n".join(lines) + "\n")
        return

    results = {}
    for test in args.tests:
        results[test] = run_enrichment(test, **shared)
        results[test].to_csv(os.path.join(out_dir, f"{test}_terms.csv"), index=False)
    write_summary(os.path.join(out_dir, "summary.txt"), results, args, foreground.shape[1])
    print(f"Wrote results to {os.path.relpath(out_dir, repo)}")


if __name__ == "__main__":
    main()
