# Flag table scripts

These scripts build `data/flag_table_genes.csv` and `data/flag_table_cells.csv`: suspected cross-sample contamination in the nine suspect samples of BioProject PRJDB7399 (the species in `data/hub_species_to_exclude.txt`). What the tables contain, the strength rules and the caveats are in [`data/flag_table_README.txt`](../../../data/flag_table_README.txt). The tables are used by the log odds ratio test (`--flagged-cells`, see [the stage 04 README](../README.md#cell-level-flags-suspected-contamination)) and the phyloGLM permutation script.

The tables were first built from interactive commands (2026-11-02). These scripts are those commands, cleaned up and given arguments instead of fixed paths.

## Steps

| Step | Script | What it does |
|---|---|---|
| 1 | `01_scan_suspects_vs_all.sh` | DIAMOND scan of the suspect samples' proteins against all 102 transcriptomes' proteins; a second scan with low-complexity masking off for queries with no hit. Run on a server. |
| 2 | `02_select_candidates.py` | Keeps matches to a different family at >= 98% protein identity over >= 100 aa; counts the species sharing each gene (`k`); keeps the best partner per gene. Writes `candidates.tsv` and the ID lists for step 3. |
| 3 | `03_extract_sequences.sh` | `seqkit` extraction of the proteins and coding sequences of the candidate genes and their partners. |
| 4 | `04_build_flag_tables.py` | Coding-sequence identity from a protein-guided codon alignment, strength labels (`--kmax 20`), gene-level table, cell-level table. |

Inputs not in the repo: all transcriptome proteins (`all_pep.faa`, IDs `<species>__<gene>`) and the TransDecoder coding sequences from stage 01.

```bash
S=scripts/04_permulation_loss_dup/flag_table
bash $S/01_scan_suspects_vs_all.sh data/hub_species_to_exclude.txt PATH_TO_ALL_PEP_FAA PATH_TO_FLAG_WORKDIR 3
python $S/02_select_candidates.py --scan PATH_TO_FLAG_WORKDIR/suspects_vs_all.tsv PATH_TO_FLAG_WORKDIR/unaligned_vs_all.tsv \
  --accessions data/txptome_accessions.csv --out-dir PATH_TO_FLAG_WORKDIR/candidates
bash $S/03_extract_sequences.sh PATH_TO_FLAG_WORKDIR/candidates PATH_TO_ALL_PEP_FAA PATH_TO_TRANSDECODER_CDS PATH_TO_FLAG_WORKDIR/seqs
python $S/04_build_flag_tables.py --candidates PATH_TO_FLAG_WORKDIR/candidates/candidates.tsv \
  --proteins PATH_TO_FLAG_WORKDIR/seqs/flag_proteins.faa --cds PATH_TO_FLAG_WORKDIR/seqs/flag_cds.fa \
  --n5-tsv data/N5.tsv --suspects data/hub_species_to_exclude.txt --accessions data/txptome_accessions.csv \
  --out-dir PATH_TO_OUTPUT
```

Steps 2 and 4 need `pandas` and `biopython`; steps 1 and 3 need `diamond` and `seqkit`.

## What was checked

- Steps 2 and 4 were re-run on the saved DIAMOND tables and sequences: `flag_table_genes.csv` (16,063 rows) and `flag_table_cells.csv` (12,350 rows) came out identical to the committed files (apart from the renamed `partner_is_suspect` column).
- Steps 1 and 3 were written from the commands run on the `will` server and were not re-run, so they are untested as scripts. The original run keyed the sequence extraction by gene ID (not the full FASTA header), as `03_extract_sequences.sh` does.
- `data/txptome_accessions.csv` lists all five *Phrurolithus* species under Corinnidae, but the original run treated *P. coreanus* as Phrurolithidae and the scripts keep that as the default (`--family-override`). The other four stay Corinnidae, so a congener would count as a different family from *P. coreanus*. No row of the committed table pairs *P. coreanus* with a congener, so this has no effect on it; check it if the scan is re-run.
- The thresholds (98% protein, 99%/97% nucleotide, 20 species) are first choices and have not been calibrated per species pair with BUSCO orthologs.
