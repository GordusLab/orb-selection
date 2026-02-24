library(RERconverge)
rerpath = find.package('RERconverge')

treefile = "~/orb-selection/assets/SpeciesTree_rooted.txt"

#define the root species
root_sp = "Antrodiaetus_roretzi"
masterTree = readTrees(treefile, max.read = 200)



# #perform binary CC permulation
# permCC = getPermsBinary(100, marineFg, sisters_marine, root_sp, mamRERw, toyTrees,
#                         masterTree, permmode="cc")
# #calculate permulation p-values
# permpvalCC = permpvalcor(res,permCC)