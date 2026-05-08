library(readr)
library(dplyr)
library(tidyr)
library(ape)

# Command line arguments: 1) loss or duplication, 2) permulation array number
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript olm_permulate_prep.R <loss|dup> <permulation_array_number>")
}

test <- args[1]
array_id <- as.integer(args[2])

# Read in species tree
treefile = "/home/crunnel2/orb-selection/data/SpeciesTree_full_brlen.nwk"
speciesTree <- read.tree(treefile)

# Read gene count table
gene_counts <- read_tsv("/home/crunnel2/orb-selection/data/N5.GeneCount.tsv")

# Read permulation tip values csv
perm_df <- read_csv("/home/crunnel2/orb-selection/data/perms_tip_values.csv")
current_perm <- perm_df[array_id, , drop = FALSE]

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
    orb_weaving = as.integer(current_perm[1, species])
  ) %>%
  mutate(
    HOG = as.character(HOG),
    species = as.character(species),
    gene_count = as.integer(gene_count),
    orb_weaving = as.logical(orb_weaving)
  )
long_df <- as.data.frame(long_df)

# Add binary variable to long_df
if (test == "loss") {
  long_df$gene_copy_var <- ifelse(long_df$gene_count == 0, 1, 0)
} else if (test == "dup") {
  long_df$gene_copy_var <- ifelse(long_df$gene_count > 1, 1, 0)
} else {
  stop("First argument must be 'loss' or 'dup'")
}

common_species <- intersect(speciesTree$tip.label, long_df$species)
speciesTree_pruned <- ape::drop.tip(speciesTree, setdiff(speciesTree$tip.label, common_species))

# Parallelized phyloglm regressions for gene loss

# Create temp directory for results
tmp_dir <- sprintf("/scratch4/agordus1/crunnel2/tmp_permulate_%s_%d", test, array_id)
if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)

saveRDS(long_df, file = sprintf("/scratch4/agordus1/crunnel2/tmp_permulate_%s_%d/long_df.rds", test, array_id))
saveRDS(speciesTree_pruned, file = sprintf("/scratch4/agordus1/crunnel2/tmp_permulate_%s_%d/speciesTree_pruned.rds", test, array_id))
