library(phylolm)
library(foreach)
library(doParallel)

# Command line arguments: 1) loss or duplication, 2) permulation array number, 3) number of cpus to use in parallel
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript olm_permulate.R <loss|dup> <permulation_array_number> <CPUs>")
}

test <- args[1]
array_id <- as.integer(args[2])
n_cores <- as.integer(args[3])

# Read in species tree
speciesTree_pruned <- readRDS(file = sprintf("/scratch4/agordus1/crunnel2/tmp_permulate_%s_%d/speciesTree_pruned.rds", test, array_id))
# Read in long format data frame
long_df <- readRDS(file = sprintf("/scratch4/agordus1/crunnel2/tmp_permulate_%s_%d/long_df.rds", test, array_id))
# Set up temporary directory for results
tmp_dir <- sprintf("/scratch4/agordus1/crunnel2/tmp_permulate_%s_%d", test, array_id)
if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)

unique_genes <- unique(long_df$HOG)

# Set up parallel backend
cl <- makeCluster(n_cores)
registerDoParallel(cl)

foreach(g = unique_genes, .packages = c("phylolm", "phytools", "ape")) %dopar% {
  df_gene <- subset(long_df, HOG == g)
  tab <- table(df_gene$gene_copy_var)
  out <- list(HOG = g)
  if (length(tab) < 2 || any(tab < 2)) {
    out$error <- "Insufficient variation"
  } else {
    rownames(df_gene) <- df_gene$species
    fit <- tryCatch(
      phyloglm(gene_copy_var ~ orb_weaving, data = df_gene, phy = speciesTree_pruned, method = "poisson_GEE"),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      out$error <- fit$message
    } else {
      out$error <- NA
      s <- summary(fit)
      coef_table <- as.data.frame(s$coefficients)
      out$z_value <- coef_table["orb_weavingTRUE", "z.value"]
      out$p_value <- coef_table["orb_weavingTRUE", "p.value"]
    }
  }
  tmpfile <- file.path(tmp_dir, paste0(g, ".rds"))
  saveRDS(out, tmpfile)
  cat(g, file = progress_file, append = TRUE, sep = "\n")
  NULL
}

stopCluster(cl)

cat(sprintf("Finished parallelized phyloglm regressions for gene loss, permulation number %d.\n", array_id))

# Note: For a binary predictor like orb_weaving, the coefficient for orb_weavingTRUE is the effect size (log-odds or log-rate ratio) for TRUE vs FALSE.
# The effect for orb_weaving == FALSE is the negative of this coefficient, and the p-value is the same.