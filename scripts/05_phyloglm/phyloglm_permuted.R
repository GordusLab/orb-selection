#!/usr/bin/env Rscript
# Permuted-label null calibration for the phyloglm test (Option A).
#
# Fits the same per-gene model as phyloglm.R (gene_count ~ orb_weaving, Poisson GEE
# on the species tree) for a fixed random subset of genes under many permuted
# phenotype labelings, and compares genome-wide summary statistics from the
# real labeling to their permuted distribution. Every labeling is applied to the
# SAME genes, so labeling-level correlation between genes is kept (a labeling
# that puts orb-weavers in a gene-rich clade shifts many genes at once).
#
# Per labeling: fraction of genes with p <= 0.05, fraction with a positive
# coefficient, median |z|, and Storey's pi0 (lambda = 0.5). Also saves the full
# z/p matrices (genes x labelings) for any further analysis.
#
# Usage:
#   Rscript phyloglm_permuted.R --outdir DIR [--n-genes 300] [--n-perms 200]
#                               [--seed 1] [--cores 4]
#                               [--perms data/perms_tip_values_100k.csv]
#
# Resumable: finished labelings (DIR/perm_<id>.rds) are skipped on re-run.
# perms file: perm_id + one 0/1 column per tip, 1 = orb-weaver (same layout as
# data/perms_tip_values_100k.csv).

suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(phylolm)
  library(ape)
  library(parallel)
  library(here)
})

# --- args ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
outdir <- get_arg("--outdir")
if (is.null(outdir)) stop("--outdir is required")
n_genes <- as.integer(get_arg("--n-genes", 300))
n_perms <- as.integer(get_arg("--n-perms", 200))
seed <- as.integer(get_arg("--seed", 1))
cores <- as.integer(get_arg("--cores", 4))
perms_file <- get_arg("--perms", here("data", "perms_tip_values_100k.csv"))
# --busco: add standardized BUSCO fraction complete + fraction duplicated as covariates
use_busco <- "--busco" %in% args
# --genes FILE: reuse a previous genes.txt so two runs compare on the same genes
genes_file <- get_arg("--genes")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- data (same preprocessing as phyloglm.R) -----------------------------
speciesTree <- read.tree(here("data", "SpeciesTree_full_brlen.nwk"))
orb_weavers <- read_lines(here("data", "orbweavers-list.txt"))

gc <- read_tsv(here("data", "N5.GeneCount.tsv"), show_col_types = FALSE)
species_cols <- gc %>%
  select(-HOG) %>%
  select(where(~ any(. != 0))) %>%
  colnames()
gc <- gc %>%
  select(HOG, all_of(species_cols)) %>%
  filter(Occupancy >= 30) %>%
  select(-OG, -`Gene Tree Parent Clade`, -Occupancy)

tips <- intersect(speciesTree$tip.label, setdiff(colnames(gc), "HOG"))
tree <- drop.tip(speciesTree, setdiff(speciesTree$tip.label, tips))
counts <- as.matrix(gc[, tips])
rownames(counts) <- gc$HOG

# per-species BUSCO covariates (columns as in odds_ratio_test.py), standardized over tips
busco <- read.csv(here("data", "buscos.csv"), header = FALSE,
                  col.names = c("species", "single", "dup", "total", "f_sc", "f_total"))
rownames(busco) <- busco$species
busco <- busco[tips, ]
busco_cov <- data.frame(
  f_total = as.numeric(scale(busco$f_total)),
  f_dup = as.numeric(scale(busco$dup / (busco$single + busco$dup))),
  row.names = tips
)
model_formula <- if (use_busco) count ~ orb_weaving + f_total + f_dup else count ~ orb_weaving

# --- single fit ----------------------------------------------------------
# returns c(z, p, converged); NA if the fit errors
fit_gene <- function(g, orb) {
  y <- as.integer(counts[g, ])
  if (length(unique(y)) < 2) return(c(NA, NA, NA))
  df <- data.frame(count = y, orb_weaving = orb[tips], busco_cov[tips, ], row.names = tips)
  converged <- TRUE
  fit <- tryCatch(
    withCallingHandlers(
      phyloglm(model_formula, data = df, phy = tree, method = "poisson_GEE"),
      warning = function(w) {
        if (grepl("converge", conditionMessage(w))) converged <<- FALSE
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(NA, NA, NA))
  co <- summary(fit)$coefficients["orb_weavingTRUE", ]
  c(co[["z.value"]], co[["p.value"]], as.numeric(converged))
}

fit_all <- function(genes, orb) {
  out <- t(vapply(genes, fit_gene, numeric(3), orb = orb))
  colnames(out) <- c("z", "p", "converged")
  out
}

# --- choose genes: random, restricted to genes that fit under the real labeling
set.seed(seed)
real_orb <- setNames(tips %in% orb_weavers, tips)
candidates <- sample(rownames(counts)[apply(counts, 1, function(r) length(unique(r)) > 1)])
if (!is.null(genes_file)) {
  # sample() above still runs so the labeling draw below is identical across runs
  genes <- read_lines(genes_file)
  real_fit <- fit_all(genes, real_orb)
} else {
  real_fit <- fit_all(candidates[seq_len(min(length(candidates), ceiling(n_genes * 1.3)))], real_orb)
  ok <- !is.na(real_fit[, "p"])
  genes <- rownames(real_fit)[ok][seq_len(n_genes)]
  real_fit <- real_fit[genes, ]
}
write_lines(genes, file.path(outdir, "genes.txt"))
saveRDS(real_fit, file.path(outdir, "real_fit.rds"))
cat(length(genes), "genes selected\n")

# --- permuted labelings --------------------------------------------------
perms <- read.csv(perms_file, check.names = FALSE)
stopifnot(all(tips %in% colnames(perms)))
perms <- perms[sample(nrow(perms), min(n_perms, nrow(perms))), ]
stopifnot(!anyNA(perms[, tips]))

run_perm <- function(i) {
  id <- perms$perm_id[i]
  f <- file.path(outdir, paste0("perm_", id, ".rds"))
  if (file.exists(f)) return(invisible(NULL))
  orb <- setNames(as.logical(unlist(perms[i, tips])), tips)
  saveRDS(fit_all(genes, orb), f)
}
cat("Running", nrow(perms), "labelings on", cores, "cores\n")
invisible(mclapply(seq_len(nrow(perms)), run_perm, mc.cores = cores, mc.preschedule = FALSE))

# --- summarize -----------------------------------------------------------
pi0 <- function(p, lambda = 0.5) min(1, mean(p > lambda, na.rm = TRUE) / (1 - lambda))
stats <- function(m, coef_sign = NULL) {
  p <- m[, "p"]
  c(frac_sig = mean(p <= 0.05, na.rm = TRUE),
    frac_pos = mean(m[, "z"] > 0, na.rm = TRUE),
    median_abs_z = median(abs(m[, "z"]), na.rm = TRUE),
    pi0 = pi0(p),
    frac_nonconverged = mean(m[, "converged"] == 0, na.rm = TRUE))
}

files <- list.files(outdir, pattern = "^perm_.*\\.rds$", full.names = TRUE)
null_fits <- lapply(files, readRDS)
names(null_fits) <- sub("^perm_|\\.rds$", "", basename(files))
null_stats <- as.data.frame(t(sapply(null_fits, stats)))
null_stats$perm_id <- names(null_fits)
write_csv(null_stats, file.path(outdir, "null_stats.csv"))
saveRDS(list(z = sapply(null_fits, function(m) m[, "z"]), p = sapply(null_fits, function(m) m[, "p"])),
        file.path(outdir, "null_matrices.rds"))

real_stats <- stats(real_fit)
cat(sprintf("\n%d genes, %d labelings\n", length(genes), nrow(null_stats)))
for (s in names(real_stats)) {
  v <- null_stats[[s]]
  cat(sprintf("%-18s real=%7.3f  null mean=%7.3f [%.3f, %.3f]  p(null>=real)=%.3f  p(null<=real)=%.3f\n",
              s, real_stats[[s]], mean(v), min(v), max(v),
              (1 + sum(v >= real_stats[[s]])) / (1 + length(v)),
              (1 + sum(v <= real_stats[[s]])) / (1 + length(v))))
}
