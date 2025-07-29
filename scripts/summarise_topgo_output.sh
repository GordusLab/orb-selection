#!/bin/bash

# Get the root directory of the git repository
repo_root=$(git rev-parse --show-toplevel)

# generate new summary files
# two criteria to pass:
# 1. At least 5 terms in the universe ($4)
# 2. P value < 0.05 OR if it contains the character '<' (because R outputs stuff like '< 1e-30') ($7)
for test in "$repo_root"/results/go_enrichment/topgo_results/*/; do
  # Remove trailing slash and get just the directory name
  test_name=$(basename "$test")
  echo "Processing directory: $test_name"
  
  cd "$test"
  # remove previous summary files
  rm -f summary_*.txt
  rm -f network_*.txt

  for a in bp*.txt; do
      echo -e 'Description\tCount\tp\tOntology' >> ${a/bp_/summary_}
      
      for b in $a ${a/bp/cc} ${a/bp/mf}; do

          # Create summary that can be used for network plots
          echo -e 'GO.ID\tDescription\tp.Val' >> network_${b}
    
          awk -v ont=$(basename $b | cut -d'_' -f1) -F $'\t' '
            BEGIN{OFS="\t"}
            {
              if ($1 ~ /[0-9]/ && $4 >= 5 && ($7 < 0.05 || $7 ~ /^</)) 
              print $3, $5, $7, ont
          }
          ' $b >> ${a/bp_/summary_}

          awk -F $'\t' '
            BEGIN{OFS="\t"}
            {
              if ($1 ~ /[0-9]/)
              print $2, $3, $7
          }
          ' $b >> network_${b}

      done
  done
done
