#!/usr/bin/env Rscript
# Merge extra permuted labelings for a gene subset into the per-gene empirical
# p-values from summarize_phyloglm_permuted.R.
#
# The subset was re-run with phyloglm_permuted_chunk.R --genes FILE on labelings
# beyond the first batch (e.g. rows 10,001+ of data/perms_tip_values_500k.csv).
# For each subset gene, the exceedances (#labelings with |z_perm| >= |z_real|) and
# valid-labeling counts from the original run and the extension are summed, and
# p_emp = (1 + exceed) / (1 + n_valid). All other genes keep their original p_emp;
# BH is then redone across every gene with a p_emp, so q_emp is comparable to the
# original table.
#
# Genes still looking strong (fewer than --stop-exceed exceedances) are written to
# DIR/genes_next.txt for a further stage; the rest cannot reach FDR significance
# however many labelings are added.
#
# Usage:
#   Rscript summarize_phyloglm_extend.R --orig-dir DIR --ext-dir DIR [--stop-exceed 20]
# Reads ORIG/gene_empirical_p.csv, ORIG/chunk_00000.rds, EXT/genes.txt, EXT/chunk_*.rds;
# writes EXT/gene_empirical_p_extended.csv and EXT/genes_next.txt.

suppressMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
orig_dir <- get_arg("--orig-dir")
ext_dir <- get_arg("--ext-dir")
if (is.null(orig_dir) || is.null(ext_dir)) stop("--orig-dir and --ext-dir are required")
stop_exceed <- as.integer(get_arg("--stop-exceed", 20))

tab <- as.data.frame(read_csv(file.path(orig_dir, "gene_empirical_p.csv"), show_col_types = FALSE))
genes <- read_lines(file.path(ext_dir, "genes.txt"))
stopifnot(all(genes %in% tab$HOG))

z_real <- readRDS(file.path(orig_dir, "chunk_00000.rds"))$z[genes, 1]
stopifnot(isTRUE(all.equal(unname(z_real), tab$z_real[match(genes, tab$HOG)])))

exceed <- tab$exceed[match(genes, tab$HOG)]
n_valid <- tab$n_valid[match(genes, tab$HOG)]
names(exceed) <- names(n_valid) <- genes

files <- sort(list.files(ext_dir, pattern = "^chunk_[0-9]+\\.rds$", full.names = TRUE))
n_labelings <- 0
for (f in files) {
  z <- readRDS(f)$z
  stopifnot(identical(rownames(z), genes))
  valid <- !is.na(z) & !is.na(z_real)
  exceed <- exceed + rowSums(valid & (abs(z) >= abs(z_real)))
  n_valid <- n_valid + rowSums(valid)
  n_labelings <- n_labelings + ncol(z)
}
cat(sprintf("%d genes, %d extra labelings from %d chunks\n", length(genes), n_labelings, length(files)))

idx <- match(genes, tab$HOG)
tab$p_emp_10k <- tab$p_emp
tab$exceed[idx] <- exceed
tab$n_valid[idx] <- n_valid
tab$p_emp[idx] <- (1 + exceed) / (1 + n_valid)
tab$q_emp <- p.adjust(tab$p_emp, method = "BH")
write_csv(tab, file.path(ext_dir, "gene_empirical_p_extended.csv"))

nxt <- genes[exceed < stop_exceed]
write_lines(nxt, file.path(ext_dir, "genes_next.txt"))

ext <- tab[idx, ]
ext <- ext[order(ext$p_emp), ]
cat(sprintf("p_emp <= 0.05: %d of %d tested genes; q_emp <= 0.05: %d; q_emp <= 0.1: %d; min p_emp = %.2g\n",
            sum(tab$p_emp <= 0.05, na.rm = TRUE), sum(!is.na(tab$p_emp)),
            sum(tab$q_emp <= 0.05, na.rm = TRUE), sum(tab$q_emp <= 0.1, na.rm = TRUE), min(tab$p_emp, na.rm = TRUE)))
cat(sprintf("extended genes with < %d exceedances (continue to next stage): %d of %d\n", stop_exceed, length(nxt), length(genes)))
cat("labelings per gene (median):", median(n_valid), "\n")
print(head(ext[, c("HOG", "z_real", "p_model", "p_emp_10k", "exceed", "n_valid", "p_emp", "q_emp")], 12), row.names = FALSE)
