library(RERconverge)
library(phangorn)

rerpath <- find.package('RERconverge')

treefile <- "~/orb-selection/assets/SpeciesTree_full_brlen.nwk"

#define the root species
root <- "Drosophila_melanogaster"

#read in the species tree
speciesTree <- read.tree(treefile)
# speciesTree <- drop.tip(
#   speciesTree, 
#   c("Antrodiaetus_roretzi", "Orchestina_okitsui", "Falcileptoneta_japonica", "Masirana_silvicola")
# )

#make empty list - don't need gene trees for permulation
testTrees <- list()

#add species tree to "list" of trees
testTrees$masterTree<-speciesTree 

#vector of foreground tips (orb-weavers)
orbweavers <- readLines("~/orb-selection/assets/orbweavers-list.txt")
allspecs <- testTrees$masterTree$tip.label
notN5 <- c("Drosophila_melanogaster", "Antrodiaetus_roretzi", "Orchestina_okitsui", "Falcileptoneta_japonica", "Masirana_silvicola")
N5 <- allspecs[!allspecs %in% notN5]
phenvec <- ifelse(N5 %in% orbweavers, 2, 1)
names(phenvec) <- N5

testCatPerm10000<-categoricalPermulations(
  testTrees, 
  phenvec, 
  rm="ER", 
  rp="auto", 
  ntrees=10000
  )
