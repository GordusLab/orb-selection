"""Adapted from David Emms OrthoFinder tool orthogroup_gene_count.py"""

import os
import sys
import csv

def main(inFN):

    outFN = os.path.splitext(inFN)[0] + ".GeneCount.tsv"
    with open(inFN, 'r') as infile, open(outFN, 'w') as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")
        header = next(reader)
        n_col_skip = 3 if header[0] == "HOG" else 1
        writer.writerow(header)
        for line in reader:
            writer.writerow(line[:n_col_skip] + [0 if "" == cell else len(cell.split(", ")) for cell in line[n_col_skip:]])
    print("Orthogroup gene count table has been written to %s" % outFN)

    return outFN

if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] == "-h" or sys.argv[1] == "-help" or sys.argv[1] == "--help":
        print("Usage: orthogroup_gene_count.py N#.tsv")
        sys.exit()

    hog_file = sys.argv[1]

    main(hog_file)
