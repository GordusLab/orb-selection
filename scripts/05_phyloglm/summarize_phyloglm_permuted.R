#!/usr/bin/env Rscript
# Combine phyloglm_permuted_chunk.R outputs into per-gene empirical p-values and
# genome-wide null statistics.
#
# Per gene: p_emp = (1 + #labelings with |z_perm| >= |z_real|) / (1 + #valid
# labelings), then Benjamini-Hochberg q_emp. Labelings where the permuted fit
# failed (NA z) are not counted for that gene.
# Per labeling (genome-wide): fraction of genes with model p <= alpha, fraction
# of positive z, median |z|, pi0 (lambda = 0.5), fraction non-converged -- for the
# real labeling and for each permuted one, plus an empirical p for each.
#
# Usage:
#   Rscript summarize_phyloglm_permuted.R --dir DIR [--alpha 0.05]
# Reads DIR/chunk_00000.rds (real) and DIR/chunk_*.rds (permuted); writes
# DIR/gene_empirical_p.csv, DIR/null_stats.csv and DIR/real_stats.csv.
# Chunks are read one at a time, so memory use stays at one chunk.

suppressMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
dir <- get_arg("--dir")
if (is.null(dir)) stop("--dir is required")
alpha <- as.numeric(get_arg("--alpha", 0.05))

# p-values are recomputed from z (the GEE summary uses a normal reference)
z_to_p <- function(z) 2 * pnorm(-abs(z))
pi0 <- function(p, lambda = 0.5) min(1, mean(p > lambda, na.rm = TRUE) / (1 - lambda))
stats <- function(z, conv) {
  p <- z_to_p(z)
  c(frac_sig = mean(p <= alpha, na.rm = TRUE),
    frac_pos = mean(z > 0, na.rm = TRUE),
    median_abs_z = median(abs(z), na.rm = TRUE),
    pi0 = pi0(p),
    frac_nonconverged = mean(conv == 0, na.rm = TRUE))
}

real <- readRDS(file.path(dir, "chunk_00000.rds"))
genes <- rownames(real$z)
z_real <- real$z[, 1]
real_stats <- stats(z_real, real$converged[, 1])

exceed <- n_valid <- setNames(integer(length(genes)), genes)
null_rows <- list()
files <- setdiff(list.files(dir, pattern = "^chunk_[0-9]+\\.rds$", full.names = TRUE),
                 file.path(dir, "chunk_00000.rds"))
for (f in sort(files)) {
  ch <- readRDS(f)
  stopifnot(identical(rownames(ch$z), genes))
  valid <- !is.na(ch$z) & !is.na(z_real)
  exceed <- exceed + rowSums(valid & (abs(ch$z) >= abs(z_real)))
  n_valid <- n_valid + rowSums(valid)
  for (j in seq_len(ncol(ch$z))) {
    null_rows[[length(null_rows) + 1]] <-
      c(perm_id = colnames(ch$z)[j], stats(ch$z[, j], ch$converged[, j]))
  }
}
null_stats <- as.data.frame(do.call(rbind, null_rows), stringsAsFactors = FALSE)
for (col in setdiff(names(null_stats), "perm_id")) null_stats[[col]] <- as.numeric(null_stats[[col]])

gene_tab <- data.frame(
  HOG = genes,
  z_real = z_real,
  p_model = z_to_p(z_real),
  exceed = exceed,
  n_valid = n_valid,
  p_emp = ifelse(is.na(z_real), NA, (1 + exceed) / (1 + n_valid))
)
gene_tab$q_emp <- p.adjust(gene_tab$p_emp, method = "BH")
write_csv(gene_tab, file.path(dir, "gene_empirical_p.csv"))
write_csv(null_stats, file.path(dir, "null_stats.csv"))
write_csv(as.data.frame(t(real_stats)), file.path(dir, "real_stats.csv"))

cat(sprintf("%d genes, %d permuted labelings\n", length(genes), nrow(null_stats)))
for (s in names(real_stats)) {
  v <- null_stats[[s]]
  cat(sprintf("%-18s real=%7.3f  null mean=%7.3f [%.3f, %.3f]  p(null>=real)=%.4f  p(null<=real)=%.4f\n",
              s, real_stats[[s]], mean(v), min(v), max(v),
              (1 + sum(v >= real_stats[[s]])) / (1 + length(v)),
              (1 + sum(v <= real_stats[[s]])) / (1 + length(v))))
}
ok <- !is.na(gene_tab$p_emp)
cat(sprintf("\nmodel p <= %.2f: %d genes;  p_emp <= %.2f: %d;  q_emp <= 0.05: %d;  q_emp <= 0.1: %d;  min p_emp = %.2g\n",
            alpha, sum(gene_tab$p_model <= alpha, na.rm = TRUE), alpha,
            sum(gene_tab$p_emp <= alpha, na.rm = TRUE),
            sum(gene_tab$q_emp <= 0.05, na.rm = TRUE), sum(gene_tab$q_emp <= 0.1, na.rm = TRUE),
            min(gene_tab$p_emp, na.rm = TRUE)))
