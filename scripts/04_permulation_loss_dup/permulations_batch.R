# Generate one batch of categorical permulations and write its tip values to CSV.
# Several batches can run as separate processes at once (see run_permulation_batches.sh)
# and be concatenated afterwards with combine_permulation_batches.py.
#
# Usage: Rscript permulations_batch.R <batch_id> [ntrees=10000] [outdir=data/perms_batches]

library(RERconverge)
library(here)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: permulations_batch.R <batch_id> [ntrees] [outdir]")
batch_id <- as.integer(args[1])
ntrees <- if (length(args) >= 2) as.integer(args[2]) else 10000L
outdir <- if (length(args) >= 3) args[3] else here("data/perms_batches")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Defines run_categorical_permulations() and save_tip_values(); does not run the 10,000-tree default
source(here("scripts/04_permulation_loss_dup/permulations.R"))

# Distinct, reproducible seed per batch
set.seed(20260000 + batch_id)

testCatPerms <- run_categorical_permulations(
    foreground_list_filename = here("data/orbweavers-list.txt"),
    ntrees = ntrees
)

# Write to a temp name first so an interrupted run never leaves a partial CSV that looks finished
final_csv <- file.path(outdir, sprintf("batch_%02d_tip_values.csv", batch_id))
tmp_csv <- paste0(final_csv, ".tmp")
save_tip_values(testCatPerms, tmp_csv)
file.rename(tmp_csv, final_csv)
message("Wrote ", final_csv)
