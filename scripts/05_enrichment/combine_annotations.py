#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine individual HOG annotation TSV files into a single master file."
    )
    parser.add_argument("input_dir", help="Directory containing the individual *_blasted.tsv files.")
    parser.add_argument("output_file", help="Path for the final combined TSV file.")
    return parser.parse_args()

def main():
    args = parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory not found at '{args.input_dir}'")
        return

    all_files = [os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir) if f.endswith('_blasted.tsv')]
    
    if not all_files:
        print(f"No '*_blasted.tsv' files found in '{args.input_dir}'")
        return

    df_list = []
    print(f"Found {len(all_files)} files to combine. Reading...")
    for f in tqdm(all_files, desc="Combining files"):
        try:
            df = pd.read_csv(f, sep='\t')
            df_list.append(df)
        except pd.errors.EmptyDataError:
            print(f"Warning: Skipping empty file {f}")
        except Exception as e:
            print(f"Warning: Could not read {f}. Error: {e}")

    if not df_list:
        print("No valid dataframes to combine. Exiting.")
        return

    print("Concatenating dataframes...")
    combined_df = pd.concat(df_list, ignore_index=True)
    
    # Ensure output directory exists
    output_dir = os.path.dirname(args.output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    print(f"Saving combined file to {args.output_file}...")
    combined_df.to_csv(args.output_file, sep='\t', index=False)
    print("Done.")

if __name__ == "__main__":
    main()
