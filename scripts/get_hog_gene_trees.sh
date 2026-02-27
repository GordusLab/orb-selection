results_dir_path="/media/will/wyrm/kelvin_data/orthofinder/run5/OrthoFinder/Results_Jan06_1"
gene_trees_path="${results_dir_path}/Resolved_Gene_Trees"
out_path=${results_dir_path}/N5_HOG_Gene_Trees
HOG_tsv_dir_path="${results_dir_path}/Phylogenetic_Hierarchical_Orthogroups"

mkdir -p ${out_path}

tail -n +2 ${HOG_tsv_dir_path}/N5.tsv | while IFS=$'\t' read -r HOG OG NODE rest; do
  gotree subtree -i ${gene_trees_path}/${OG}_tree.txt -n "^${NODE}$" > ${out_path}/${HOG}_tree.txt
done