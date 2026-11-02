# Generate one batch of categorical permulations and write its tip values to CSV.
# Several batches can run as separate processes at once (see run_permulation_batches.sh)
# and be concatenated afterwards with combine_permulation_batches.py.
#
# Usage: Rscript permulations_batch.R <batch_id> [ntrees=10000] [outdir=data/perms_batches] \
#            [exclude_file=NONE] [seed_base=20260000]
#   exclude_file: text file with one species per line to drop from the tree (in addition to the default
#                 excluded tips), e.g. for a robustness analysis without some species. Foreground species
#                 that are dropped simply no longer count.
#   seed_base:    the seed is seed_base + batch_id; use a different base for a different tree/analysis.

library(RERconverge)
library(here)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: permulations_batch.R <batch_id> [ntrees] [outdir]")
batch_id <- as.integer(args[1])
ntrees <- if (length(args) >= 2) as.integer(args[2]) else 10000L
outdir <- if (length(args) >= 3) args[3] else here("data/perms_batches")
exclude_file <- if (length(args) >= 4) args[4] else "NONE"
seed_base <- if (length(args) >= 5) as.integer(args[5]) else 20260000L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Defines run_categorical_permulations() and save_tip_values(); does not run the 10,000-tree default
source(here("scripts/04_permulation_loss_dup/permulations.R"))

# Distinct, reproducible seed per batch
set.seed(seed_base + batch_id)

# Tips dropped from the tree: the defaults of run_categorical_permulations plus any extra species
excluded_tips <- eval(formals(run_categorical_permulations)$excluded_tips)
if (exclude_file != "NONE") {
    extra <- trimws(readLines(exclude_file))
    excluded_tips <- unique(c(excluded_tips, extra[nzchar(extra)]))
    message("Excluding ", length(excluded_tips), " tips in total (", length(extra), " from ", exclude_file, ")")
}

testCatPerms <- run_categorical_permulations(
    foreground_list_filename = here("data/orbweavers-list.txt"),
    ntrees = ntrees,
    excluded_tips = excluded_tips
)

# Write to a temp name first so an interrupted run never leaves a partial CSV that looks finished
final_csv <- file.path(outdir, sprintf("batch_%02d_tip_values.csv", batch_id))
tmp_csv <- paste0(final_csv, ".tmp")
save_tip_values(testCatPerms, tmp_csv)
file.rename(tmp_csv, final_csv)
message("Wrote ", final_csv)
