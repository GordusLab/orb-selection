# Pipeline adapted from https://github.com/lyijin/topGO_pipeline/
library(here)
library(topGO)


folders <- dir(here("results/go_enrichment/loc_lists"), full.names = TRUE)
annot_filename <- here("data/udiv_go_annots.all.tsv")

for (folder in folders) {
  # hit files are all files without "universe" in the name
  hit_files <- grep(
    list.files(folder),
    pattern = "*universe*",
    inv = TRUE,
    value = TRUE
  )

  # shrink list of all GO terms down to the correct universe
  universe_file <- grep(
    list.files(folder),
    pattern = "*universe*",
    value = TRUE
  )

  universe_genes <- scan(
    paste0(folder, "/", universe_file),
    character(0),
    sep = "\n"
  )

  for (hits in hit_files) {
    gene_id_to_go <- readMappings(file = annot_filename)
    gene_id_to_go <- gene_id_to_go[universe_genes]
    gene_id_to_go <- gene_id_to_go[gene_id_to_go != "no_hit"]
    gene_names <- names(gene_id_to_go)

    for (go_category in c("bp", "cc", "mf")) {
      print(paste("Current file:", hits))
      genes_of_interest_filename <- paste0(folder, "/", hits)
      genes_of_interest <- scan(
        genes_of_interest_filename,
        character(0),
        sep = "\n"
      )

      genelist <- factor(as.integer(gene_names %in% genes_of_interest))
      names(genelist) <- gene_names

      go_data <- try(new(
        "topGOdata",
        ontology = toupper(go_category),
        allGenes = genelist,
        gene2GO = gene_id_to_go,
        annotationFun = annFUN.gene2GO
      ))

      # handle error
      if (class(go_data) == "try-error") {
        print(paste0("Error for file", hits, "!"))
        next
      }

      # weight01 is the default algorithm used in Alexa et al. (2006)
      weight01_fisher <- runTest(go_data, statistic = "fisher")

      # generate a results table (for only the top 1000 GO terms)
      #   topNodes: highest 1000 GO terms shown
      #   numChar: truncates GO term descriptions at 1000 chars
      #   (basically, disables truncation)
      if (length(genes_of_interest) < 500) {
        results_table <- GenTable(
          go_data,
          P_value = weight01_fisher,
          orderBy = "P_value",
          topNodes = 100,
          numChar = 1000
        )
      } else {
        results_table <- GenTable(
          go_data,
          P_value = weight01_fisher,
          orderBy = "P_value",
          topNodes = 300,
          numChar = 1000
        )
      }

      # write it out into a file for python post-processing
      # Remove _hits.txt from the filename to get clean base name
      base_name <- sub("_hits\\.txt$", "", basename(hits))

      output_dir <- here(
        "results/go_enrichment/topgo_results", base_name
      )
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }

      output_filename <- file.path(output_dir, paste0(go_category, "_", hits))

      write.table(
        results_table,
        file = output_filename,
        quote = FALSE,
        sep = "\t"
      )
    }
  }
}
