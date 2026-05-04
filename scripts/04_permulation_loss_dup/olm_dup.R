library(readr)
library(dplyr)
library(tidyr)
library(phylolm)
library(here)
library(ape)
library(foreach)
library(doParallel)

# Read in species tree
treefile = "/home/crunnel2/orb-selection/data/SpeciesTree_full_brlen.nwk"
speciesTree <- read.tree(treefile)

# Read gene count table
gene_counts <- read_tsv("/home/crunnel2/orb-selection/data/N5.GeneCount.tsv")

# Read orb-weaver species list
orb_weavers <- read_lines("/home/crunnel2/orb-selection/data/orbweavers-list.txt")

# Filter rows with occupancy >= 30
gene_counts <- gene_counts %>%
  filter(Occupancy >= 30)

# Remove unwanted columns
gene_counts <- gene_counts %>%
  select(-OG, -`Gene Tree Parent Clade`, -Occupancy)

# Drop species columns with all zeros
species_cols <- gene_counts %>%
  select(-HOG) %>%
  select(where(~ any(. != 0))) %>%
  colnames()

gene_counts <- gene_counts %>%
  select(HOG, all_of(species_cols))


# Pivot to long format and convert to data.frame
long_df <- gene_counts %>%
  pivot_longer(
    cols = -HOG,
    names_to = "species",
    values_to = "gene_count"
  ) %>%
  mutate(
    orb_weaving = species %in% orb_weavers
  ) %>%
  mutate(
    HOG = as.character(HOG),
    species = as.character(species),
    gene_count = as.integer(gene_count),
    orb_weaving = as.logical(orb_weaving)
  )
long_df <- as.data.frame(long_df)

# Add binary variable to long_df
long_df$gene_duplicated <- ifelse(long_df$gene_count > 1, 1, 0)

common_species <- intersect(speciesTree$tip.label, long_df$species)
speciesTree_pruned <- ape::drop.tip(speciesTree, setdiff(speciesTree$tip.label, common_species))

# Parallelized phyloglm regressions for gene duplication

# Create temp directory for results
tmp_dup_dir <- "/scratch4/agordus1/crunnel2/tmp_dup"
if (!dir.exists(tmp_dup_dir)) dir.create(tmp_dup_dir, recursive = TRUE)

# Helper to combine temp files into a CSV
combine_tmp_results <- function(tmp_dir, out_csv) {
  files <- list.files(tmp_dir, full.names = TRUE, pattern = "\\.rds$")
  results <- lapply(files, readRDS)
  if (length(results) == 0) return()
  all_fields <- unique(unlist(lapply(results, names)))
  df <- as.data.frame(do.call(rbind, lapply(results, function(x) {
    row <- setNames(rep(NA, length(all_fields)), all_fields)
    for (n in names(x)) row[[n]] <- x[[n]]
    row
  })))
  write.csv(df, file = out_csv, row.names = FALSE)
}

# Set up parallel backend
n_cores <- parallel::detectCores()
cl <- makeCluster(n_cores)
registerDoParallel(cl)

unique_genes <- unique(long_df$HOG)

# Progress file setup
progress_file <- "/scratch4/agordus1/crunnel2/phyloglm_progress_dup.txt"
if (file.exists(progress_file)) file.remove(progress_file)
file.create(progress_file)

# Progress monitor function (runs in main process)
progress_monitor <- function(total, interval = 300) {
  repeat {
    Sys.sleep(interval)
    done <- 0
    if (file.exists(progress_file)) {
      done <- length(readLines(progress_file))
    }
    cat(sprintf("[Progress] %d of %d genes completed at %s\n", done, total, format(Sys.time(), "%H:%M:%S")))
    if (done >= total) break
    # If a stop file exists, exit early (for test runs or manual stop)
    if (file.exists(paste0(progress_file, ".stop"))) break
  }
}

cat("Timing gene duplication regressions for all genes...\n")
progress_pid2 <- parallel:::mcparallel(progress_monitor(length(unique_genes), interval = 300))
dup_time <- system.time({
  foreach(g = unique_genes, .packages = c("phylolm", "phytools", "ape")) %dopar% {
    df_gene <- subset(long_df, HOG == g)
    tab <- table(df_gene$gene_duplicated)
    out <- list(HOG = g)
    if (length(tab) < 2 || any(tab < 2)) {
      out$error <- "Insufficient variation"
      out$fit <- NA
    } else {
      rownames(df_gene) <- df_gene$species
      fit <- tryCatch(
        phyloglm(gene_duplicated ~ orb_weaving, data = df_gene, phy = speciesTree_pruned, method = "poisson_GEE"),
        error = function(e) e
      )
      if (inherits(fit, "error")) {
        out$error <- fit$message
        out$fit <- NA
      } else {
        out$error <- NA
        out$fit <- NA
        s <- summary(fit)
        coef_table <- as.data.frame(s$coefficients)
        for (rn in rownames(coef_table)) {
          for (cn in colnames(coef_table)) {
            out[[paste0("coef_", rn, "_", cn)]] <- coef_table[rn, cn]
          }
        }
        out$logLik <- fit$logLik
        out$alpha <- fit$alpha
        out$converged <- fit$converged
      }
    }
    tmpfile <- file.path(tmp_dup_dir, paste0(g, ".rds"))
    saveRDS(out, tmpfile)
    cat(g, file = progress_file, append = TRUE, sep = "\n")
    NULL
  }
})

# Wait for progress monitor to finish, but always kill if main loop is done
try({
  parallel:::mccollect(progress_pid2, wait = FALSE)
  Sys.sleep(1)
  if (parallel:::selectChildren(list(progress_pid2), timeout = 0)[[1]] == 0) {
    file.create(paste0(progress_file, ".stop"))
    parallel:::mckill(progress_pid2, signal = 9)
    parallel:::mccollect(progress_pid2)
  }
}, silent = TRUE)

cat("Elapsed time for gene duplication regressions (user, system, elapsed):\n")
print(dup_time)

stopCluster(cl)

cat("Finished parallelized phyloglm regressions for gene loss and duplication.\n")

# Combine temp files into final CSV
combine_tmp_results(tmp_dup_dir, here("results/phyloglm_dup.csv"))

# Note: For a binary predictor like orb_weaving, the coefficient for orb_weavingTRUE is the effect size (log-odds or log-rate ratio) for TRUE vs FALSE.
# The effect for orb_weaving == FALSE is the negative of this coefficient, and the p-value is the same.