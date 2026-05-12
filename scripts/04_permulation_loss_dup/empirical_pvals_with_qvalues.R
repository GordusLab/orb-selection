# Calculate one-tailed empirical p-values and q-values for continuous phyloglm permulations
# Uses parallel processing over genes

library(foreach)
library(doParallel)

# --- User settings ---
perm_dirs <- sprintf("/scratch4/agordus1/crunnel2/tmp_permulate_cont_%d", 1:10)
obs_csv <- "/home/crunnel2/orb-selection/results/phyloglm_continuous.csv"
out_csv <- "/home/crunnel2/orb-selection/results/phyloglm_continuous_perm_pvals.csv"
tmp_results_dir <- "/scratch4/agordus1/crunnel2/pval_results_tmp"
n_cores <- as.numeric(commandArgs(trailingOnly = TRUE)[1])  # Adjust as needed

# Use whichever z-value column exists in each file format.
z_col_candidates <- c("coef_orb_weavingTRUE_z.value", "z_value")

extract_z <- function(x, candidates) {
  for (nm in candidates) {
    if (!is.null(x[[nm]])) return(as.numeric(x[[nm]]))
  }
  NA_real_
}

# --- Load observed z-values ---
obs <- read.csv(obs_csv, stringsAsFactors = FALSE)
obs_z <- sapply(seq_len(nrow(obs)), function(i) extract_z(obs[i, , drop = FALSE], z_col_candidates))
genes <- obs$HOG
names(obs_z) <- genes

# Create temp directory for results
if (!dir.exists(tmp_results_dir)) dir.create(tmp_results_dir, recursive = TRUE)

# --- Set up parallel backend ---
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# --- Parallelize over genes ---
foreach(gene = genes, .packages = "base") %dopar% {
  # For each gene, loop through all permutation directories and accumulate counts
  orb_cnt <- 0L
  non_orb_cnt <- 0L
  valid_cnt <- 0L
  
  z_obs <- obs_z[gene]
  
  for (perm_dir in perm_dirs) {
    rds_file <- file.path(perm_dir, paste0(gene, ".rds"))
    if (!file.exists(rds_file)) next
    
    out_perm <- readRDS(rds_file)
    z_perm <- extract_z(out_perm, z_col_candidates)
    
    if (!is.na(z_perm) && !is.na(z_obs)) {
      valid_cnt <- valid_cnt + 1L
      if (z_perm >= z_obs) orb_cnt <- orb_cnt + 1L
      if (z_perm <= z_obs) non_orb_cnt <- non_orb_cnt + 1L
    }
  }
  
  # Calculate p-values with pseudocount
  orb_p <- (orb_cnt + 1) / (valid_cnt + 1)
  non_orb_p <- (non_orb_cnt + 1) / (valid_cnt + 1)
  
  if (valid_cnt == 0) {
    orb_p <- NA_real_
    non_orb_p <- NA_real_
  }
  
  # Save results to temp file
  gene_result <- list(
    HOG = gene,
    z_obs = as.numeric(z_obs),
    orb_count = as.integer(orb_cnt),
    non_orb_count = as.integer(non_orb_cnt),
    n_perm_valid = as.integer(valid_cnt),
    orb_p = as.numeric(orb_p),
    non_orb_p = as.numeric(non_orb_p)
  )
  
  tmpfile <- file.path(tmp_results_dir, paste0(gene, ".rds"))
  saveRDS(gene_result, tmpfile)
  NULL
}

stopCluster(cl)

cat(sprintf("Finished parallelized empirical p-value calculations for %d genes.\n", length(genes)))

# --- Load and combine all temp results ---
result_list <- list()
for (gene in genes) {
  tmpfile <- file.path(tmp_results_dir, paste0(gene, ".rds"))
  if (file.exists(tmpfile)) {
    result_list[[gene]] <- readRDS(tmpfile)
  }
}

# Convert list to data frame
result <- data.frame(
  HOG = sapply(result_list, function(x) x$HOG),
  z_obs = sapply(result_list, function(x) x$z_obs),
  orb_count = sapply(result_list, function(x) x$orb_count),
  non_orb_count = sapply(result_list, function(x) x$non_orb_count),
  n_perm_valid = sapply(result_list, function(x) x$n_perm_valid),
  orb_p = sapply(result_list, function(x) x$orb_p),
  non_orb_p = sapply(result_list, function(x) x$non_orb_p),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# --- Calculate q-values ---
if (!requireNamespace("qvalue", quietly = TRUE)) {
  install.packages("qvalue")
}
library(qvalue)

result$orb_q <- NA_real_
result$non_orb_q <- NA_real_

idx_orb <- which(!is.na(result$orb_p))
idx_non_orb <- which(!is.na(result$non_orb_p))

if (length(idx_orb) > 0) {
  result$orb_q[idx_orb] <- qvalue(result$orb_p[idx_orb])$qvalues
}
if (length(idx_non_orb) > 0) {
  result$non_orb_q[idx_non_orb] <- qvalue(result$non_orb_p[idx_non_orb])$qvalues
}

# --- Write final results ---
write.csv(result, out_csv, row.names = FALSE)

cat(sprintf("Wrote results to %s\n", out_csv))

# Optional: clean up temp files
unlink(tmp_results_dir, recursive = TRUE)
