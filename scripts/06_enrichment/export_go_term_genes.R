# Export GO term -> gene sets for a reference genome, with annotations propagated up the
# GO hierarchy (a gene annotated to a term is also annotated to all of its ancestors), the same
# way topGO does it. Used by permulation_enrichment.py.
#
# Usage (from the repo root):
#   Rscript scripts/06_enrichment/export_go_term_genes.R [annots=data/udiv_go_annots.all.tsv] \
#       [out=data/udiv_go_term_gene_sets.tsv]
#
# Input: topGO gene-to-GO mapping file (gene <tab> comma-separated GO ids).
# Output columns: GOID, ontology (BP/MF/CC), term, n_genes, genes (comma-separated).

suppressMessages({
  library(here)
  library(topGO)
  library(GO.db)
})

args <- commandArgs(trailingOnly = TRUE)
annot_file <- if (length(args) >= 1) args[1] else here("data/udiv_go_annots.all.tsv")
out_file <- if (length(args) >= 2) args[2] else here("data/udiv_go_term_gene_sets.tsv")

gene_to_go <- readMappings(file = annot_file)
gene_to_go <- gene_to_go[names(gene_to_go) != "no_hit"]

direct <- unique(data.frame(
  gene = rep(names(gene_to_go), lengths(gene_to_go)),
  go = unlist(gene_to_go, use.names = FALSE),
  stringsAsFactors = FALSE
))
direct <- direct[direct$go != "no_hit" & grepl("^GO:", direct$go), ]

term_info <- AnnotationDbi::select(
  GO.db,
  keys = unique(direct$go),
  columns = c("TERM", "ONTOLOGY"),
  keytype = "GOID"
)
n_dropped <- length(setdiff(unique(direct$go), term_info$GOID))
message(n_dropped, " annotated GO ids not found in GO.db (obsolete) were dropped")
direct <- direct[direct$go %in% term_info$GOID, ]
ontology_of <- setNames(term_info$ONTOLOGY, term_info$GOID)

ancestors <- list(
  BP = as.list(GOBPANCESTOR),
  MF = as.list(GOMFANCESTOR),
  CC = as.list(GOCCANCESTOR)
)

expanded <- list()
for (ont in names(ancestors)) {
  d <- direct[ontology_of[direct$go] == ont, ]
  anc <- lapply(d$go, function(g) setdiff(c(g, ancestors[[ont]][[g]]), "all"))
  expanded[[ont]] <- unique(data.frame(
    gene = rep(d$gene, lengths(anc)),
    go = unlist(anc, use.names = FALSE),
    ontology = ont,
    stringsAsFactors = FALSE
  ))
}
expanded <- do.call(rbind, expanded)

genes_by_term <- split(expanded$gene, expanded$go)
term_names <- AnnotationDbi::select(
  GO.db,
  keys = names(genes_by_term),
  columns = c("TERM", "ONTOLOGY"),
  keytype = "GOID"
)

out <- data.frame(
  GOID = term_names$GOID,
  ontology = term_names$ONTOLOGY,
  term = term_names$TERM,
  n_genes = lengths(genes_by_term[term_names$GOID]),
  genes = vapply(
    genes_by_term[term_names$GOID],
    function(g) paste(sort(unique(g)), collapse = ","),
    character(1)
  ),
  stringsAsFactors = FALSE
)
out <- out[order(out$ontology, -out$n_genes), ]

write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(
  "Wrote ", nrow(out), " terms (",
  paste(names(table(out$ontology)), table(out$ontology), collapse = ", "), ") to ", out_file
)
