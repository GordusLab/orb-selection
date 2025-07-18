"""
Functions to filter orthogroups based on occupancy,
presence of Uloborus diversus, and single copy representation."""

import argparse
import pandas as pd


def drop_empty_cols(df, print_txt=False):
    """Drops columns where all entries (ignoring headers) are 0,
    to get rid of species not included in the current node/hierarchy"""

    # Print the number of columns before cleaning
    num_columns_before = df.shape[1]
    print(f"Number of columns before dropping empty columns: {num_columns_before}")

    # Drop columns where all entries (ignoring headers) are 0
    df_cleaned = df.loc[:, (df.ne(0)).any(axis=0)]

    # Print the number of columns after cleaning
    num_columns_after = df_cleaned.shape[1]

    if print_txt:
        print(f"Number of columns after dropping empty columns: {num_columns_after}")
        print("Species with no sequences in any orthogroup have been dropped.")

    return df_cleaned


def udiv_present(df):
    """Filters the DataFrame to keep only orthogroups with Uloborus diversus present."""

    num_rows_before = df.shape[0]
    print(f"Number of orthogroups before filtering for U.div: {num_rows_before}")

    df_udiv_present = df[(df["Uloborus_diversus"] > 0)]
    num_rows_after = df_udiv_present.shape[0]

    print(f"Number of orthogroups after filtering for U.div: {num_rows_after}")
    print("Orthogroups with no U. div representation have been dropped.")

    return df_udiv_present


def set_occupancy_threshold(df, occ):
    """Filters the DataFrame to keep only orthogroups with occupancy above a specified threshold."""

    num_rows_before = df.shape[0]
    print(f"Number of orthogroups before filtering for occupancy: {num_rows_before}")
    numeric_cols = df.select_dtypes(include="number").columns

    # #add occupancy column
    # df['Occupancy'] = (df[numeric_cols] > 0).sum(axis=1)

    df_occupancy = df[(df[numeric_cols] > 0).sum(axis=1) >= occ]
    num_rows_after = df_occupancy.shape[0]

    print(f"Number of orthogroups after filtering for occupancy: {num_rows_after}")
    print(f"Orthogroups with occupancy less than {occ} have been dropped.")

    return df_occupancy


def single_copy_filter(df, threshold):
    """Set a minimum percentage of single copy species"""

    num_rows_before = df.shape[0]
    print(
        f"Number of orthogroups before filtering for single copy reps: {num_rows_before}"
    )

    numeric_cols = df.select_dtypes(include="number").columns

    # #add single copy reps column
    # df['No. single copy reps'] = (df[numeric_cols] == 1).sum(axis=1)
    df_single = df[(df[numeric_cols] == 1).sum(axis=1) >= threshold]
    num_rows_after = df_single.shape[0]

    print(
        f"Number of orthogroups after filtering for single copy reps: {num_rows_after}"
    )
    print(
        f"Orthogroups with less than {threshold}\
           species present in single copy have been dropped."
    )

    return df_single


def main():
    """Main function to parse arguments and filter the gene count DataFrame."""
    parser = argparse.ArgumentParser()
    parser.add_argument("genecount_tsv_file")
    parser.add_argument("output_filename")
    parser.add_argument("-o", "--occupancy", type=int)
    parser.add_argument("-u", "--udiv_present", action="store_true")
    parser.add_argument("-s", "--single_copy_threshold", type=int)
    parser.add_argument("-l", "--save_hog_list", action="store_true")
    args = parser.parse_args()

    # Read the CSV file into a DataFrame
    df = pd.read_csv(args.genecount_tsv_file, sep="\t", index_col="HOG")

    df = drop_empty_cols(df)

    if args.udiv_present is True:
        df = udiv_present(df)
    else:
        pass

    if args.single_copy_threshold is not None:
        df = single_copy_filter(df, args.single_copy_threshold)
    else:
        pass

    if args.occupancy is not None:
        df = set_occupancy_threshold(df, args.occupancy)
    else:
        pass

    # Save the cleaned DataFrame to a new TSV file
    df.to_csv(f"{args.output_filename}.tsv", sep="\t")
    print(f"Filtered GeneCount file saved as {args.output_filename}.tsv.")

    if args.save_hog_list is True:
        hogs = list(df.index.values)
        with open(f"{args.output_filename}_list.txt", "w", encoding="utf-8") as f:
            for line in hogs:
                f.write(f"{line}\n")
        print(f"Filtered HOG list saved as {args.output_filename}_list.txt.")

    else:
        pass


if __name__ == "__main__":
    main()
