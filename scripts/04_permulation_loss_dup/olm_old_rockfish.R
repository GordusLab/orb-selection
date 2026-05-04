library(readr)
library(dplyr)
library(tidyr)
library(phylolm)
# library(here)
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

# Add binary variables to long_df
long_df$gene_lost <- ifelse(long_df$gene_count == 0, 1, 0)
long_df$gene_duplicated <- ifelse(long_df$gene_count > 1, 1, 0)

common_species <- intersect(speciesTree$tip.label, long_df$species)
speciesTree_pruned <- ape::drop.tip(speciesTree, setdiff(speciesTree$tip.label, common_species))


# Parallelized phyloglm regressions for gene loss and duplication

# Set up parallel backend
n_cores <- parallel::detectCores()
cl <- makeCluster(n_cores)
registerDoParallel(cl)

unique_genes <- unique(long_df$HOG)

# Progress file setup
progress_file <- "./phyloglm_progress.txt"
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


cat("Timing gene loss regressions for all genes...\n")
# Start progress monitor in background
progress_pid <- parallel:::mcparallel(progress_monitor(length(unique_genes), interval = 300))

loss_time <- system.time({
  results_loss <- foreach(g = unique_genes, .packages = c("phylolm", "phytools", "ape")) %dopar% {
    df_gene <- subset(long_df, HOG == g)
    tab <- table(df_gene$gene_lost)
    if (length(tab) < 2 || any(tab < 2)) {
      # Mark as done in progress file
      cat(g, file = progress_file, append = TRUE, sep = "\n")
      return(list(HOG = g, fit = NA, error = "Insufficient variation"))
    }
    # Set rownames for phyloglm
    rownames(df_gene) <- df_gene$species
    fit <- tryCatch(
      phyloglm(gene_lost ~ orb_weaving, data = df_gene, phy = speciesTree_pruned, method = "poisson_GEE"),
      error = function(e) e
    )
    # Mark as done in progress file
    cat(g, file = progress_file, append = TRUE, sep = "\n")
    if (inherits(fit, "error")) {
      return(list(HOG = g, fit = NA, error = fit$message))
    }
    list(HOG = g, fit = fit, error = NA)
  }
})


# Wait for progress monitor to finish, but always kill if main loop is done
try({
  parallel:::mccollect(progress_pid, wait = FALSE)
  Sys.sleep(1) # Give monitor a chance to print final update
  if (parallel:::selectChildren(list(progress_pid), timeout = 0)[[1]] == 0) {
    # Still running, create stop file and kill
    file.create(paste0(progress_file, ".stop"))
    parallel:::mckill(progress_pid, signal = 9)
    parallel:::mccollect(progress_pid)
  }
}, silent = TRUE)

cat("Elapsed time for gene loss regressions (user, system, elapsed):\n")
print(loss_time)


cat("Timing gene duplication regressions for all genes...\n")
# Reset progress file
if (file.exists(progress_file)) file.remove(progress_file)
file.create(progress_file)
progress_pid2 <- parallel:::mcparallel(progress_monitor(length(unique_genes), interval = 300))

dup_time <- system.time({
  results_dup <- foreach(g = unique_genes, .packages = c("phylolm", "phytools", "ape")) %dopar% {
    df_gene <- subset(long_df, HOG == g)
    tab <- table(df_gene$gene_duplicated)
    if (length(tab) < 2 || any(tab < 2)) {
      # Mark as done in progress file
      cat(g, file = progress_file, append = TRUE, sep = "\n")
      return(list(HOG = g, fit = NA, error = "Insufficient variation"))
    }
    # Set rownames for phyloglm
    rownames(df_gene) <- df_gene$species
    fit <- tryCatch(
      phyloglm(gene_duplicated ~ orb_weaving, data = df_gene, phy = speciesTree_pruned, method = "poisson_GEE"),
      error = function(e) e
    )
    # Mark as done in progress file
    cat(g, file = progress_file, append = TRUE, sep = "\n")
    if (inherits(fit, "error")) {
      return(list(HOG = g, fit = NA, error = fit$message))
    }
    list(HOG = g, fit = fit, error = NA)
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

# Extract coefficients and p-values and save to CSV
loss_coef <- sapply(results_loss, function(res) {
  fit <- res$fit
  if (inherits(fit, "phyloglm")) {
    coef_table <- summary(fit)$coefficients
    if ("orb_weavingTRUE" %in% rownames(coef_table)) {
      return(coef_table["orb_weavingTRUE", "Estimate"])
    }
  }
  return(NA)
})

loss_pvals <- sapply(results_loss, function(res) {
  fit <- res$fit
  if (inherits(fit, "phyloglm")) {
    coef_table <- summary(fit)$coefficients
    if ("orb_weavingTRUE" %in% rownames(coef_table)) {
      return(coef_table["orb_weavingTRUE", "p.value"])
    }
  }
  return(NA)
})

dup_coef <- sapply(results_dup, function(res) {
  fit <- res$fit
  if (inherits(fit, "phyloglm")) {
    coef_table <- summary(fit)$coefficients
    if ("orb_weavingTRUE" %in% rownames(coef_table)) {
      return(coef_table["orb_weavingTRUE", "Estimate"])
    }
  }
  return(NA)
})

dup_pvals <- sapply(results_dup, function(res) {
  fit <- res$fit
  if (inherits(fit, "phyloglm")) {
    coef_table <- summary(fit)$coefficients
    if ("orb_weavingTRUE" %in% rownames(coef_table)) {
      return(coef_table["orb_weavingTRUE", "p.value"])
    }
  }
  return(NA)
})

# Combine results into a data.frame
results_df <- data.frame(
  HOG = unique_genes,
  loss_coef = as.numeric(loss_coef),
  loss_pval = as.numeric(loss_pvals),
  dup_coef = as.numeric(dup_coef),
  dup_pval = as.numeric(dup_pvals)
)

# Save to CSV
write.csv(results_df, file = "/home/crunnel2/orb-selection/results/phyloglm_pvalues.csv", row.names = FALSE)

# Note: For a binary predictor like orb_weaving, the coefficient for orb_weavingTRUE is the effect size (log-odds or log-rate ratio) for TRUE vs FALSE.
# The effect for orb_weaving == FALSE is the negative of this coefficient, and the p-value is the same.