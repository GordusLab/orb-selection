library(RERconverge)
library(here)

rerpath <- find.package("RERconverge")

run_categorical_permulations <- function(
  foreground_list_filename,
  ntrees = 10000,
  treefile = here("data/SpeciesTree_full_brlen.nwk"),
  excluded_tips = c(
      "Drosophila_melanogaster",
      "Antrodiaetus_roretzi",
      "Orchestina_okitsui",
      "Falcileptoneta_japonica",
      "Masirana_silvicola"
  ),
  rm = "ER",
  rp = "auto",
  save_rdata_path = NULL,
  save_object_name = "testCatPerms"
) {
    # Read species tree and wrap in a list as required by categoricalPermulations.
    speciesTree <- read.tree(treefile)
    testTrees <- list(masterTree = speciesTree)

    # Build phenotype vector: foreground tips as 2, background as 1.
    foreground <- readLines(foreground_list_filename)
    allspecs <- testTrees$masterTree$tip.label
    included_tips <- allspecs[!allspecs %in% excluded_tips]
    phenvec <- ifelse(included_tips %in% foreground, 2, 1)
    names(phenvec) <- included_tips

    testCatPerms <- categoricalPermulations(
        testTrees,
        phenvec,
        rm = rm,
        rp = rp,
        ntrees = ntrees
    )

    # Optionally save permulations object to an RData file for downstream inspection.
    if (!is.null(save_rdata_path)) {
        dir.create(dirname(save_rdata_path), recursive = TRUE, showWarnings = FALSE)
        assign(save_object_name, testCatPerms)
        save(list = save_object_name, file = save_rdata_path)
        message("Saved permulations RData to: ", save_rdata_path)
    }

    return(testCatPerms)
}

save_tip_values <- function(testCatPerms, output_csv_path) {
    # Canonical species order from the first permulation
    tip_order <- names(testCatPerms$trees[[1]]$tips)

    # Optional safety checks
    stopifnot(all(vapply(
        testCatPerms$trees,
        function(tr) setequal(names(tr$tips), tip_order),
        logical(1)
    )))

    # Build matrix: 10,000 rows (permulations) x N species columns
    tip_mat <- do.call(
        rbind,
        lapply(testCatPerms$trees, function(tr) {
            # reorder each named tip vector to match canonical order
            tr$tips[tip_order]
        })
    )

    # Convert to data frame
    tip_df <- as.data.frame(tip_mat, check.names = FALSE)
    tip_df$perm_id <- seq_len(nrow(tip_df)) # optional row id
    tip_df <- tip_df[, c("perm_id", tip_order)] # put perm_id first

    # tip_df now has species as columns and ntrees rows of 1 and 2 values
    dim(tip_df)
    head(tip_df[, 1:6])

    tip_df[tip_order] <- tip_df[tip_order] - 1

    write.csv(tip_df, output_csv_path, row.names = FALSE)
}

# Initial run
testCatPerms <- run_categorical_permulations(
    foreground_list_filename = here("data/orbweavers-list.txt"),
    ntrees = 10000,
    save_rdata_path = here("data/perms10000.RData"),
)

save_tip_values(testCatPerms, here("data/perms_tip_values.csv"))
