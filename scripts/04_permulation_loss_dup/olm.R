library(readr)
library(dplyr)
library(tidyr)
library(phylolm)
library(here)
library(ape)

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

treefile = here("data/SpeciesTree_full_brlen.nwk")
speciesTree <- read.tree(treefile)

## Subset to the first gene in long_df
first_gene <- unique(long_df$HOG)[1]
df_gene <- subset(long_df, HOG == first_gene)

# Filter for complete cases
complete_cases <- complete.cases(df_gene$gene_count, df_gene$orb_weaving)
df_gene_complete <- df_gene[complete_cases, ]


# Debug: check for mismatches between tree and data
cat("Species in tree but not in data:\n")
print(setdiff(speciesTree$tip.label, df_gene_complete$species))

cat("Species in data but not in tree:\n")
print(setdiff(df_gene_complete$species, speciesTree$tip.label))

# Match tree and data
common_species <- intersect(speciesTree$tip.label, df_gene_complete$species)
speciesTree_pruned <- ape::drop.tip(speciesTree, setdiff(speciesTree$tip.label, common_species))
df_gene_matched <- df_gene_complete[df_gene_complete$species %in% common_species, ]

# Add binary variables to long_df
long_df$gene_lost <- ifelse(long_df$gene_count == 0, 1, 0)
long_df$gene_duplicated <- ifelse(long_df$gene_count > 1, 1, 0)

# Convert to standard data.frame
df_gene_matched <- as.data.frame(df_gene_matched)
rownames(df_gene_matched) <- df_gene_matched$species

# Run phyloglm for this gene
fit <- phylolm(gene_count ~ orb_weaving, data = df_gene_matched, phy = speciesTree_pruned)
summary(fit)

# Run phyloglm for gene_lost
cat("\nphyloglm for gene_lost (gene_count == 0):\n")
fit_lost <- phyloglm(gene_lost ~ orb_weaving, data = df_gene_matched, phy = speciesTree_pruned, method = "poisson_GEE")
print(summary(fit_lost))

# Run phyloglm for gene_duplicated
cat("\nphyloglm for gene_duplicated (gene_count > 1):\n")
fit_dup <- phyloglm(gene_duplicated ~ orb_weaving, data = df_gene_matched, phy = speciesTree_pruned, method = "poisson_GEE")
print(summary(fit_dup))