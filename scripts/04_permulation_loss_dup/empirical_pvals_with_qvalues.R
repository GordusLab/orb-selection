# Calculate empirical p-values and q-values for phyloglm permulations (loss and duplication)
# Adjust file paths and perm_dirs as needed

# --- User settings ---
perm_dirs <- sprintf("/scratch4/agordus1/crunnel2/tmp_permulate_loss_%d", 1:1000) # adjust as needed
obs_loss_csv <- "/home/crunnel2/orb-selection/results/phyloglm_loss.csv"
obs_dup_csv  <- "/home/crunnel2/orb-selection/results/phyloglm_dup.csv"
z_col <- "coef_orb_weavingTRUE_z.value"

# --- Load observed z-values ---
obs_loss <- read.csv(obs_loss_csv, row.names = 1)
obs_dup  <- read.csv(obs_dup_csv, row.names = 1)
genes <- intersect(obs_loss$HOG, obs_dup$HOG)
obs_loss_z <- setNames(obs_loss[[z_col]], obs_loss$HOG)
obs_dup_z  <- setNames(obs_dup[[z_col]],  obs_dup$HOG)

# --- Initialize counters ---
loss_orb     <- setNames(rep(0, length(genes)), genes)
loss_non_orb <- setNames(rep(0, length(genes)), genes)
dup_orb      <- setNames(rep(0, length(genes)), genes)
dup_non_orb  <- setNames(rep(0, length(genes)), genes)

# --- Loop over permulation directories ---
for (perm_dir in perm_dirs) {
  for (gene in genes) {
    # Loss
    rds_file_loss <- file.path(perm_dir, paste0(gene, ".rds"))
    if (file.exists(rds_file_loss)) {
      out_loss <- readRDS(rds_file_loss)
      z_perm_loss <- out_loss[[z_col]]
      z_obs_loss  <- obs_loss_z[gene]
      if (!is.na(z_perm_loss) && !is.na(z_obs_loss)) {
        if (z_perm_loss >= z_obs_loss) loss_orb[gene]     <- loss_orb[gene] + 1
        if (z_perm_loss <= z_obs_loss) loss_non_orb[gene] <- loss_non_orb[gene] + 1
      }
    }
    # Duplication
    rds_file_dup <- file.path(gsub("loss", "dup", perm_dir), paste0(gene, ".rds"))
    if (file.exists(rds_file_dup)) {
      out_dup <- readRDS(rds_file_dup)
      z_perm_dup <- out_dup[[z_col]]
      z_obs_dup  <- obs_dup_z[gene]
      if (!is.na(z_perm_dup) && !is.na(z_obs_dup)) {
        if (z_perm_dup >= z_obs_dup) dup_orb[gene]     <- dup_orb[gene] + 1
        if (z_perm_dup <= z_obs_dup) dup_non_orb[gene] <- dup_non_orb[gene] + 1
      }
    }
  }
}

n_perm <- length(perm_dirs)
result <- data.frame(
  HOG = genes,
  loss_orb_count     = loss_orb,
  loss_non_orb_count = loss_non_orb,
  dup_orb_count      = dup_orb,
  dup_non_orb_count  = dup_non_orb,
  n_perm = n_perm,
  loss_orb_p     = loss_orb / n_perm,
  loss_non_orb_p = loss_non_orb / n_perm,
  dup_orb_p      = dup_orb / n_perm,
  dup_non_orb_p  = dup_non_orb / n_perm
)

# --- Calculate q-values ---
if (!requireNamespace("qvalue", quietly = TRUE)) {
  install.packages("qvalue")
}
library(qvalue)

result$loss_orb_q     <- qvalue(result$loss_orb_p)$qvalues
result$loss_non_orb_q <- qvalue(result$loss_non_orb_p)$qvalues
result$dup_orb_q      <- qvalue(result$dup_orb_p)$qvalues
result$dup_non_orb_q  <- qvalue(result$dup_non_orb_p)$qvalues

write.csv(result, "/home/crunnel2/orb-selection/results/phyloglm_perm_pvals.csv", row.names = FALSE)
