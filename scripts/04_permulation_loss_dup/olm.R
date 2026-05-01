library(readr)
library(dplyr)
library(tidyr)
library(phylolm)
library(here)
library(ape)

# Read in species tree
treefile = here("data/SpeciesTree_full_brlen.nwk")
speciesTree <- read.tree(treefile)

# Read gene count table
gene_counts <- read_tsv(here("data/N5.GeneCount.tsv"))

# Read orb-weaver species list
orb_weavers <- read_lines(here("data/orbweavers-list.txt"))

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

# Pivot to long format
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


# Add binary variables to long_df
long_df$gene_lost <- ifelse(long_df$gene_count == 0, 1, 0)
long_df$gene_duplicated <- ifelse(long_df$gene_count > 1, 1, 0)

common_species <- intersect(speciesTree$tip.label, long_df$species)
speciesTree_pruned <- ape::drop.tip(speciesTree, setdiff(speciesTree$tip.label, common_species))

# Parallelized phyloglm regressions for gene loss and duplication
library(foreach)
library(doParallel)

# Set up parallel backend
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

unique_genes <- unique(long_df$HOG)[1:11]

cat("Timing gene loss regressions for first 10 genes...\n")
loss_time <- system.time({
  results_loss <- foreach(g = unique_genes, .packages = c("phylolm", "phytools", "ape")) %dopar% {
    df_gene <- subset(long_df, HOG == g)
    tab <- table(df_gene$gene_lost)
    if (length(tab) < 2 || any(tab < 2)) {
      return(list(HOG = g, fit = NA, error = "Insufficient variation"))
    }
    # Convert to standard data.frame
    df_gene <- as.data.frame(df_gene)
    rownames(df_gene) <- df_gene$species
    fit <- tryCatch(
      phyloglm(gene_lost ~ orb_weaving, data = df_gene, phy = speciesTree_pruned, method = "poisson_GEE"),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      return(list(HOG = g, fit = NA, error = fit$message))
    }
    list(HOG = g, fit = fit, error = NA)
  }
})
cat("Elapsed time for gene loss regressions (user, system, elapsed):\n")
print(loss_time)

cat("Timing gene duplication regressions for first 10 genes...\n")
dup_time <- system.time({
  results_dup <- foreach(g = unique_genes, .packages = c("phylolm", "phytools", "ape")) %dopar% {
    df_gene <- subset(long_df, HOG == g)
    tab <- table(df_gene$gene_duplicated)
    if (length(tab) < 2 || any(tab < 2)) {
      return(list(HOG = g, fit = NA, error = "Insufficient variation"))
    }
    # Convert to standard data.frame
    df_gene <- as.data.frame(df_gene)
    rownames(df_gene) <- df_gene$species
    fit <- tryCatch(
      phyloglm(gene_duplicated ~ orb_weaving, data = df_gene, phy = speciesTree_pruned, method = "poisson_GEE"),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      return(list(HOG = g, fit = NA, error = fit$message))
    }
    list(HOG = g, fit = fit, error = NA)
  }
})
cat("Elapsed time for gene duplication regressions (user, system, elapsed):\n")
print(dup_time)

stopCluster(cl)

cat("Finished parallelized phyloglm regressions for gene loss and duplication.\n")
