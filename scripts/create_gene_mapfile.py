#!/usr/bin/env python3
"""
Create a mapfile from N5.tsv where:
- First column: species_name + '_' + gene_name (with periods replaced by underscores)
- Second column: species_name

Example output:
Poecilopachys_australasia_ICKB01022636_1_p1	Poecilopachys_australasia
"""

import pandas as pd
import re

def create_mapfile(input_file, output_file):
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')
    
    # Get species names (column headers starting from column 3)
    species_columns = df.columns[3:]  # Skip HOG, OG, Gene Tree Parent Clade
    
    mapfile_entries = []
    total_rows = len(df)
    rows_with_75plus_species = 0
    
    # Process each row
    for _, row in df.iterrows():
        # Count how many species have at least one gene in this row
        species_with_genes = 0
        for species_name in species_columns:
            genes_str = str(row[species_name])
            # Check if this species has any genes (not empty or NaN)
            if genes_str != 'nan' and genes_str.strip() != '':
                species_with_genes += 1
        
        # Only process this row if at least 75 species have genes
        if species_with_genes < 75:
            continue
            
        rows_with_75plus_species += 1
            
        # Process each species column
        for species_name in species_columns:
            genes_str = str(row[species_name])
            
            # Skip empty cells (NaN becomes 'nan')
            if genes_str == 'nan' or genes_str.strip() == '':
                continue
            
            # Split by comma and space if multiple genes
            genes = [gene.strip() for gene in genes_str.split(',')]
            
            for gene in genes:
                if gene and gene != 'nan':
                    # Replace periods with underscores in gene name
                    gene_clean = gene.replace('.', '_')
                    
                    # Create the combined name: species_gene
                    combined_name = f"{species_name}_{gene_clean}"
                    
                    # Add to mapfile entries
                    mapfile_entries.append([combined_name, species_name])
    
    # Create DataFrame and save
    mapfile_df = pd.DataFrame(mapfile_entries, columns=['gene_id', 'species'])
    mapfile_df.to_csv(output_file, sep='\t', header=False, index=False)
    
    print(f"Total rows in input: {total_rows}")
    print(f"Rows with ≥75 species having genes: {rows_with_75plus_species}")
    print(f"Rows filtered out: {total_rows - rows_with_75plus_species}")
    print(f"Created mapfile with {len(mapfile_entries)} entries")
    print(f"Output saved to: {output_file}")
    
    # Show first few entries as example
    print("\nFirst 10 entries:")
    for i in range(min(10, len(mapfile_entries))):
        print(f"{mapfile_entries[i][0]}\t{mapfile_entries[i][1]}")

if __name__ == "__main__":
    input_file = "/Users/calvin/orb-selection/assets/N5.tsv"
    output_file = "/Users/calvin/orb-selection/assets/gene_species_mapfile.tsv"
    
    create_mapfile(input_file, output_file)
