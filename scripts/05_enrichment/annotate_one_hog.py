#!/usr/bin/env python3

import argparse
import os
import re
import tempfile
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbitblastnCommandline

DEFAULT_NO_HIT = {
    "Best_hit": "NO_HIT",
    "Subject_id": "NO_HIT",
    "Gene": "NO_HIT",
    "Db_xref": "NO_HIT",
    "Protein": "NO_HIT",
    "Protein_id": "NO_HIT",
    "Location": "NO_HIT",
    "Gbkey": "NO_HIT",
    "E_val": "NO_HIT",
    "iden": "NO_HIT",
    "len": "NO_HIT",
    "Coords": "NO_HIT",
}

def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate a single Orthofinder group by blasting its genes against a reference BLAST DB."
    )
    parser.add_argument("hog", help="The HOG identifier to process.")
    parser.add_argument("n0_tsv", help="Orthofinder N0 TSV file.")
    parser.add_argument("ortho_dir", help="Directory containing per-species FASTA files.")
    parser.add_argument("ref_db", help="Reference BLAST database path.")
    parser.add_argument("output_dir", help="Directory to save the output TSV file.")
    return parser.parse_args()

def load_species_sequences(ortho_dir, species):
    infasta = os.path.join(ortho_dir, f"{species}.faa")
    if not os.path.exists(infasta):
        return {}
    return SeqIO.to_dict(SeqIO.parse(infasta, "fasta"))

def blast_gene(seq_record, ref_db):
    with tempfile.TemporaryDirectory(prefix="hog_blast_") as tmpdir:
        query_path = os.path.join(tmpdir, "query.fasta")
        xml_path = os.path.join(tmpdir, "blast.xml")

        SeqIO.write(seq_record, query_path, "fasta")

        blast_cline = NcbitblastnCommandline(
            query=query_path,
            db=ref_db,
            evalue=1e-5,
            outfmt=5,
            out=xml_path,
            num_threads=1,
        )
        _, stderr = blast_cline()

        if stderr:
            print(f"BLAST ERROR for {seq_record.id}: {stderr}")

        if not os.path.exists(xml_path) or os.path.getsize(xml_path) == 0:
            return DEFAULT_NO_HIT.copy()

        with open(xml_path, encoding="utf-8", errors="replace") as handle:
            for record in NCBIXML.parse(handle):
                if not record.alignments:
                    continue

                align = record.alignments[0]
                hsps = align.hsps[0]
                parsed = parse_subject_header(align.title)
                return {
                    "Best_hit": align.title,
                    "Subject_id": parsed["subject_id"],
                    "Gene": parsed["gene"],
                    "Db_xref": parsed["db_xref"],
                    "Protein": parsed["protein"],
                    "Protein_id": parsed["protein_id"],
                    "Location": parsed["location"],
                    "Gbkey": parsed["gbkey"],
                    "E_val": hsps.expect,
                    "iden": hsps.identities * 3,
                    "len": hsps.align_length * 3,
                    "Coords": str([hsps.sbjct_start, hsps.sbjct_end]),
                }
    return DEFAULT_NO_HIT.copy()

def parse_subject_header(title):
    subject_id = str(title).split(" [", 1)[0].strip()
    bracket_fields = {
        key: value
        for key, value in re.findall(r"\[(\w+)=([^\]]+)\]", str(title))
    }
    return {
        "subject_id": subject_id or "NO_HIT",
        "gene": bracket_fields.get("gene", "NO_HIT"),
        "db_xref": bracket_fields.get("db_xref", "NO_HIT"),
        "protein": bracket_fields.get("protein", "NO_HIT"),
        "protein_id": bracket_fields.get("protein_id", "NO_HIT"),
        "location": bracket_fields.get("location", "NO_HIT"),
        "gbkey": bracket_fields.get("gbkey", "NO_HIT"),
    }

def main():
    args = parse_args()
    
    df = pd.read_csv(args.n0_tsv, sep="\t", index_col="HOG")
    
    if args.hog not in df.index:
        print(f"HOG {args.hog} not found in {args.n0_tsv}")
        return

    hog_row = df.loc[[args.hog]].reset_index()

    result_data = DEFAULT_NO_HIT.copy()
    result_data["HOG"] = args.hog
    result_data["OG"] = hog_row.at[0, "OG"]

    species_cols = [col for col in hog_row.columns if col not in ["HOG", "OG", "Gene Tree Parent Clade"]]

    for species in species_cols:
        if pd.notna(hog_row.at[0, species]):
            genes_str = str(hog_row.at[0, species])
            gene = genes_str.split(",")[0].strip()
            if not gene:
                continue

            seq_dict = load_species_sequences(args.ortho_dir, species)
            seq_record = seq_dict.get(gene)
            
            if seq_record:
                blast_result = blast_gene(seq_record, args.ref_db)
                if blast_result["Best_hit"] != "NO_HIT":
                    result_data.update(blast_result)
                    break 
    
    os.makedirs(args.output_dir, exist_ok=True)
    out_df = pd.DataFrame([result_data])
    outname = os.path.join(args.output_dir, f"{args.hog}_blasted.tsv")
    out_df.to_csv(outname, sep="\t", index=False)
    print(f"Annotation for {args.hog} saved to {outname}")

if __name__ == "__main__":
    main()
