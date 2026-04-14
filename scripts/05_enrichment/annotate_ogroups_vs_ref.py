#!/usr/bin/env python3

import argparse
import os
import re
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor

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


# Worker config set by _init_worker.
WORKER_CONFIG = {"ortho_dir": None, "ref_db": None}
WORKER_SEQ_CACHE = {}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Annotate Orthofinder groups by blasting the first available gene in each HOG "
            "against a reference BLAST DB."
        )
    )
    parser.add_argument("n0_tsv", help="Orthofinder N0 TSV file")
    parser.add_argument("ortho_dir", help="Directory containing per-species FASTA files")
    parser.add_argument("ref_db", help="Reference BLAST database path")
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 1) - 1),
        help="Number of worker processes (default: CPU count - 1)",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=25,
        help="Task chunksize for process pool mapping",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=500,
        help="Print progress every N processed HOG rows (set <=0 to disable)",
    )
    return parser.parse_args()


def _init_worker(ortho_dir, ref_db):
    WORKER_CONFIG["ortho_dir"] = ortho_dir
    WORKER_CONFIG["ref_db"] = ref_db
    WORKER_SEQ_CACHE.clear()


def _load_species_sequences(species):
    if species in WORKER_SEQ_CACHE:
        return WORKER_SEQ_CACHE[species]

    infasta = os.path.join(WORKER_CONFIG["ortho_dir"], f"{species}.faa")
    if not os.path.exists(infasta):
        WORKER_SEQ_CACHE[species] = {}
        return WORKER_SEQ_CACHE[species]

    seq_dict = SeqIO.to_dict(SeqIO.parse(infasta, "fasta"))
    WORKER_SEQ_CACHE[species] = seq_dict
    return seq_dict


def _blast_gene(seq_record):
    with tempfile.TemporaryDirectory(prefix="hog_blast_") as tmpdir:
        query_path = os.path.join(tmpdir, "query.fasta")
        xml_path = os.path.join(tmpdir, "blast.xml")

        SeqIO.write(seq_record, query_path, "fasta")

        blast_cline = NcbitblastnCommandline(
            query=query_path,
            db=WORKER_CONFIG["ref_db"],
            evalue=1e-5,
            outfmt=5,
            out=xml_path,
            num_threads=1,
        )
        _, stderr = blast_cline()

        if stderr:
            # Keep stderr visible but continue processing.
            print(f"BLAST ERROR: {stderr}")

        with open(xml_path, encoding="utf-8", errors="replace") as handle:
            for record in NCBIXML.parse(handle):
                if not record.alignments:
                    continue

                align = record.alignments[0]
                hsps = align.hsps[0]
                parsed = _parse_subject_header(align.title)
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
                    # Convert AA-based alignment metrics to nucleotide units.
                    "iden": hsps.identities * 3,
                    "len": hsps.align_length * 3,
                    "Coords": str([hsps.sbjct_start, hsps.sbjct_end]),
                }

    return DEFAULT_NO_HIT.copy()


def _process_row(task):
    row_idx, species_to_genes = task

    for species, genes_str in species_to_genes:
        gene = str(genes_str).split(",")[0].strip()
        if not gene:
            continue

        seq_dict = _load_species_sequences(species)
        seq_record = seq_dict.get(gene)
        if seq_record is None:
            continue

        row_result = _blast_gene(seq_record)
        if row_result["Best_hit"] != "NO_HIT":
            return row_idx, row_result

    return row_idx, DEFAULT_NO_HIT.copy()


def _parse_subject_header(title):
    subject_id = str(title).split(" [", 1)[0].strip()
    bracket_fields = {
        key: value
        for key, value in re.findall(r"\[(\w+)=([^\]]+)\]", str(title))
    }

    return {
        "subject_id": subject_id if subject_id else "NO_HIT",
        "gene": bracket_fields.get("gene", "NO_HIT"),
        "db_xref": bracket_fields.get("db_xref", "NO_HIT"),
        "protein": bracket_fields.get("protein", "NO_HIT"),
        "protein_id": bracket_fields.get("protein_id", "NO_HIT"),
        "location": bracket_fields.get("location", "NO_HIT"),
        "gbkey": bracket_fields.get("gbkey", "NO_HIT"),
    }


def _build_tasks(df):
    tasks = []
    species_cols = list(df.columns[3:])

    for row_idx in range(len(df)):
        present_species = []
        for species in species_cols:
            val = df.at[row_idx, species]
            if pd.notna(val):
                present_species.append((species, val))
        tasks.append((row_idx, present_species))

    return tasks


def _consume_results(results_iter, out, total_tasks, progress_every):
    start = time.time()
    processed = 0

    for row_idx, row_result in results_iter:
        for key, value in row_result.items():
            out.at[row_idx, key] = value

        processed += 1
        if progress_every > 0 and (processed % progress_every == 0 or processed == total_tasks):
            elapsed = max(time.time() - start, 1e-9)
            rate = processed / elapsed
            pct = (processed / total_tasks) * 100 if total_tasks else 100.0
            print(
                f"[progress] {processed}/{total_tasks} rows ({pct:.1f}%) | "
                f"elapsed: {elapsed:.1f}s | rate: {rate:.2f} rows/s"
            )


def main():
    args = parse_args()

    df = pd.read_csv(args.n0_tsv, sep="\t")
    out = pd.read_csv(args.n0_tsv, sep="\t", usecols=["HOG", "OG"])

    for col, default in DEFAULT_NO_HIT.items():
        out[col] = default

    tasks = _build_tasks(df)
    total_tasks = len(tasks)

    print(
        f"Starting annotation: {total_tasks} rows | workers={args.workers} | "
        f"progress_every={args.progress_every}"
    )

    if args.workers <= 1:
        _init_worker(args.ortho_dir, args.ref_db)
        results_iter = map(_process_row, tasks)
        _consume_results(results_iter, out, total_tasks, args.progress_every)
    else:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=_init_worker,
            initargs=(args.ortho_dir, args.ref_db),
        ) as ex:
            results_iter = ex.map(_process_row, tasks, chunksize=args.chunksize)
            _consume_results(results_iter, out, total_tasks, args.progress_every)

            outname = args.n0_tsv[:-4] + "_blasted.tsv"
            out.to_csv(outname, sep="\t", index=False)
            return

    outname = args.n0_tsv[:-4] + "_blasted.tsv"
    out.to_csv(outname, sep="\t", index=False)


if __name__ == "__main__":
    main()
