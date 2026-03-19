library(RERconverge)

rerpath <- find.package("RERconverge")

run_categorical_permulations <- function(
    foreground_list_filename,
    ntrees = 10000,
    treefile = "~/orb-selection/assets/SpeciesTree_full_brlen.nwk",
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
    speciesTree <- read.tree(path.expand(treefile))
    testTrees <- list(masterTree = speciesTree)

    # Build phenotype vector: foreground tips as 2, background as 1.
    foreground <- readLines(path.expand(foreground_list_filename))
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
        save_rdata_path <- path.expand(save_rdata_path)
        dir.create(dirname(save_rdata_path), recursive = TRUE, showWarnings = FALSE)
        assign(save_object_name, testCatPerms)
        save(list = save_object_name, file = save_rdata_path)
        message("Saved permulations RData to: ", save_rdata_path)
    }

    return(testCatPerms)
}


testCatPerms <- run_categorical_permulations(
        foreground_list_filename = "~/orb-selection/assets/orbweavers-list.txt",
        ntrees = 1000,
        save_rdata_path = "~/orb-selection/assets/perms.RData"
)
