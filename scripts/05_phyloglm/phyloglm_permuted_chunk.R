#!/usr/bin/env Rscript
# One chunk of the full gene x labeling phyloglm permulation null.
#
# Fits gene_count ~ orb_weaving (+ BUSCO fraction complete + fraction duplicated,
# unless --no-busco) with the same Poisson GEE as phyloglm.R, for EVERY gene
# under a block of permuted phenotype labelings (or under the real labeling),
# and saves the orb_weaving z-values and convergence flags. Summarize across
# chunks with summarize_phyloglm_permuted.R.
#
# Usage:
#   Rscript phyloglm_permuted_chunk.R --chunk K --outdir DIR
#       [--chunk-size 100] [--cores N] [--perms data/perms_tip_values_100k.csv]
#       [--genes FILE] [--perm-ids FILE] [--no-busco] [--model poisson|logistic]
#       [--flagged-cells FILE [--flag-strengths strong] [--flagged-only]]
#
#   --chunk 0     fit the real labeling (data/orbweavers-list.txt)
#   --chunk K>=1  fit labelings (K-1)*chunk-size + 1 .. K*chunk-size of the perms
#                 file, in file order (or in --perm-ids order, if given)
#   --genes       restrict to the HOGs in FILE (testing); default: every HOG with
#                 occupancy >= 30 and more than one distinct count
#   --model       poisson: count ~ orb_weaving (Poisson GEE, as phyloglm.R); logistic: presence (count > 0)
#                 ~ orb_weaving (phylogenetic logistic regression, comparable to the odds ratio loss test);
#                 for logistic the default gene set is the HOGs whose presence varies
#   --flagged-cells  gene-level contamination flag table (data/flag_table_genes.csv: species, HOG, strength).
#                 Same semantics as the odds ratio test: flagged genes are counted per HOG x species cell and
#                 subtracted from the gene count; a cell left with no copies is missing (NA) and that species
#                 is dropped from that HOG's fit (the tree is pruned for that gene only); occupancy is
#                 recomputed without the missing cells before the occupancy >= 30 filter. Unflagged genes
#                 are fit exactly as without the option.
#   --flag-strengths  comma-separated strengths to apply (default strong; e.g. strong,moderate)
#   --flagged-only   only fit HOGs with at least one missing cell
#   --cores       parallel over genes within each labeling (default:
#                 SLURM_CPUS_PER_TASK, else 1)
#
# Resumable: an existing DIR/chunk_<K>.rds is left alone.
# perms file: perm_id + one 0/1 column per tip, 1 = orb-weaver.

suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(phylolm)
  library(ape)
  library(parallel)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
chunk <- as.integer(get_arg("--chunk"))
outdir <- get_arg("--outdir")
if (is.na(chunk) || is.null(outdir)) stop("--chunk and --outdir are required")
chunk_size <- as.integer(get_arg("--chunk-size", 100))
cores <- as.integer(get_arg("--cores", Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
perms_file <- get_arg("--perms", here("data", "perms_tip_values_100k.csv"))
genes_file <- get_arg("--genes")
perm_ids_file <- get_arg("--perm-ids")
use_busco <- !("--no-busco" %in% args)
# --model poisson (gene count, Poisson GEE; default) or logistic (gene presence/absence, logistic_MPLE)
model <- get_arg("--model", "poisson")
if (!model %in% c("poisson", "logistic")) stop("--model must be poisson or logistic")
flag_file <- get_arg("--flagged-cells")
flag_strengths <- strsplit(get_arg("--flag-strengths", "strong"), ",")[[1]]
flagged_only <- "--flagged-only" %in% args
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(outdir, sprintf("chunk_%05d.rds", chunk))
if (file.exists(out_file)) {
  cat("Chunk", chunk, "already complete.\n")
  quit(save = "no")
}

# --- data (same preprocessing as phyloglm.R) -----------------------------
speciesTree <- read.tree(here("data", "SpeciesTree_full_brlen.nwk"))
orb_weavers <- read_lines(here("data", "orbweavers-list.txt"))

gc <- read_tsv(here("data", "N5.GeneCount.tsv"), show_col_types = FALSE)
species_cols <- gc %>%
  select(-HOG) %>%
  select(where(~ any(. != 0))) %>%
  colnames()
gc <- gc %>% select(HOG, all_of(species_cols))

tips <- intersect(speciesTree$tip.label, setdiff(colnames(gc), "HOG"))
tree <- drop.tip(speciesTree, setdiff(speciesTree$tip.label, tips))
counts <- as.matrix(gc[, tips])
rownames(counts) <- gc$HOG
occupancy <- gc$Occupancy

# contamination flags: subtract flagged copies; cells left with none become NA (missing, not zero)
if (!is.null(flag_file)) {
  fl <- read_csv(flag_file, show_col_types = FALSE) %>%
    filter(strength %in% flag_strengths) %>%
    count(HOG, species, name = "n_flagged")
  i <- match(fl$HOG, rownames(counts)); j <- match(fl$species, tips)
  ok <- !is.na(i) & !is.na(j)
  ok[ok] <- counts[cbind(i[ok], j[ok])] > 0
  fl <- fl[ok, ]; i <- i[ok]; j <- j[ok]
  new <- pmax(counts[cbind(i, j)] - fl$n_flagged, 0)
  counts[cbind(i, j)] <- new
  gone <- new == 0
  counts[cbind(i[gone], j[gone])] <- NA
  occupancy <- occupancy - tabulate(i[gone], nbins = length(occupancy))
  cat(sprintf("flags (%s): %d cells applied, %d left missing, %d HOGs touched\n",
              paste(flag_strengths, collapse = "+"), length(i), sum(gone), length(unique(i))))
}
counts <- counts[occupancy >= 30, , drop = FALSE]

genes <- if (!is.null(genes_file)) {
  read_lines(genes_file)
} else {
  rownames(counts)[apply(counts, 1, function(r) {
    r <- r[!is.na(r)]
    length(unique(if (model == "logistic") r > 0 else r)) > 1
  })]
}
# HOGs that fall below occupancy 30 once flagged cells are removed are no longer fit
genes <- intersect(genes, rownames(counts))
if (flagged_only) genes <- intersect(genes, rownames(counts)[rowSums(is.na(counts)) > 0])

# per-species BUSCO covariates (columns as in odds_ratio_test.py), standardized
busco <- read.csv(here("data", "buscos.csv"), header = FALSE,
                  col.names = c("species", "single", "dup", "total", "f_sc", "f_total"))
rownames(busco) <- busco$species
busco <- busco[tips, ]
busco_cov <- data.frame(
  f_total = as.numeric(scale(busco$f_total)),
  f_dup = as.numeric(scale(busco$dup / (busco$single + busco$dup))),
  row.names = tips
)
model_formula <- if (use_busco) count ~ orb_weaving + f_total + f_dup else count ~ orb_weaving

# returns c(z, converged); NA if the fit errors
fit_gene <- function(g, orb) {
  y <- as.integer(counts[g, ])
  if (model == "logistic") y <- as.integer(y > 0)
  gene_tips <- tips
  gene_tree <- tree
  if (anyNA(y)) {   # flagged cells: fit without those species
    gene_tips <- tips[!is.na(y)]
    y <- y[!is.na(y)]
    gene_tree <- drop.tip(tree, setdiff(tips, gene_tips))
  }
  if (length(unique(y)) < 2) return(c(NA, NA))
  df <- data.frame(count = y, orb_weaving = orb[gene_tips], busco_cov[gene_tips, ], row.names = gene_tips)
  converged <- TRUE
  fit <- tryCatch(
    withCallingHandlers(
      phyloglm(model_formula, data = df, phy = gene_tree,
               method = if (model == "logistic") "logistic_MPLE" else "poisson_GEE"),
      warning = function(w) {
        if (grepl("converge", conditionMessage(w))) converged <<- FALSE
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(NA, NA))
  c(summary(fit)$coefficients["orb_weavingTRUE", "z.value"], as.numeric(converged))
}

fit_labeling <- function(orb, gene_cores = 1) {
  res <- mclapply(genes, fit_gene, orb = orb, mc.cores = gene_cores)
  res <- do.call(rbind, res)
  rownames(res) <- genes
  res
}

# --- labelings for this chunk --------------------------------------------
if (chunk == 0) {
  labelings <- list(real = setNames(tips %in% orb_weavers, tips))
} else {
  perms <- read.csv(perms_file, check.names = FALSE)
  stopifnot(all(tips %in% colnames(perms)), !anyNA(perms[, tips]))
  if (!is.null(perm_ids_file)) {
    perms <- perms[match(read_lines(perm_ids_file), perms$perm_id), ]
  }
  rows <- ((chunk - 1) * chunk_size + 1):min(chunk * chunk_size, nrow(perms))
  if (length(rows) == 0 || rows[1] > nrow(perms)) stop("chunk ", chunk, " is past the end of the perms file")
  labelings <- lapply(rows, function(i) setNames(as.logical(unlist(perms[i, tips])), tips))
  names(labelings) <- perms$perm_id[rows]
}

cat(sprintf("chunk %d: %d genes x %d labelings on %d cores (model=%s, busco=%s)\n",
            chunk, length(genes), length(labelings), cores, model, use_busco))
t0 <- Sys.time()
# Each finished labeling is saved to DIR/partial/ as it completes, so a chunk killed
# by the time limit or a node failure resumes from there, and progress can be read off
# `ls DIR/partial | wc -l`. The partial files are removed once the chunk file is written.
partial_dir <- file.path(outdir, "partial")
dir.create(partial_dir, showWarnings = FALSE)
fit_labeling_cached <- function(i) {
  f <- file.path(partial_dir, sprintf("%05d_%s.rds", chunk, names(labelings)[i]))
  if (file.exists(f)) return(readRDS(f))
  res <- fit_labeling(labelings[[i]])
  saveRDS(res, paste0(f, ".tmp"))
  file.rename(paste0(f, ".tmp"), f)
  res
}

# Parallelize over labelings (each worker fits all genes serially) so cores never
# wait on the slowest gene; a single labeling (the real one) parallelizes over genes.
fits <- if (length(labelings) == 1) {
  lapply(labelings, fit_labeling, gene_cores = cores)
} else {
  mclapply(seq_along(labelings), fit_labeling_cached,
           mc.cores = min(cores, length(labelings)), mc.preschedule = FALSE)
}
names(fits) <- names(labelings)

z <- sapply(fits, function(m) m[, 1])
converged <- sapply(fits, function(m) m[, 2])
storage.mode(converged) <- "integer"
rownames(z) <- rownames(converged) <- genes

tmp <- paste0(out_file, ".tmp")
saveRDS(list(z = z, converged = converged, busco = use_busco, model = model), tmp)
file.rename(tmp, out_file)
unlink(list.files(partial_dir, pattern = sprintf("^%05d_", chunk), full.names = TRUE))
cat(sprintf("chunk %d done in %.1f min\n", chunk, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
