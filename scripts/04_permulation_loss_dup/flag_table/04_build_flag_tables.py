"""Score candidates by coding-sequence identity, label their strength, and write the flag tables.

For each suspect gene's best partner (candidates.tsv from 02_select_candidates.py) the nucleotide identity is
measured over the aligned codons: the two proteins are aligned globally (BLOSUM62, gap open -11, extend -1) and
the alignment is mapped back onto the coding sequences, comparing codon by codon.

Strength: strong = nt identity >= 99 and shared by <= KMAX species; moderate = 97 <= identity < 99 and
<= KMAX species, or identity >= 99 and > KMAX species; weak = everything else.
(Above about 20 species the medians fall to ordinary conserved-gene levels, so wide sharing is not
treated as strong evidence.)

A flagged gene copy is subtracted from its species' copy count in the HOG; a species left with none is 'missing'
(weight 0 in the odds ratio test), not 'absent'.

Writes flag_table_genes.csv and flag_table_cells.csv in --out-dir.
"""
import argparse
import collections
import csv
import os

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--candidates", required=True, help="candidates.tsv from 02_select_candidates.py")
ap.add_argument("--proteins", nargs="+", required=True, help="protein FASTA(s), IDs '<species>__<gene>'")
ap.add_argument("--cds", nargs="+", required=True, help="coding-sequence FASTA(s), IDs '<species>__<gene>'")
ap.add_argument("--n5-tsv", required=True, help="data/N5.tsv (HOG per suspect gene)")
ap.add_argument("--suspects", required=True, help="text file: the suspect species, one per line (or whitespace separated)")
ap.add_argument("--accessions", required=True, help="data/txptome_accessions.csv")
ap.add_argument("--family-override", nargs="*", default=["Phrurolithus_coreanus=Phrurolithidae"])
ap.add_argument("--kmax", type=int, default=20)
ap.add_argument("--out-dir", required=True)
args = ap.parse_args()

SUSPECTS = open(args.suspects).read().split()
fam = {r["Species"]: r["Family"] for r in csv.DictReader(open(args.accessions, encoding="utf-8-sig"))}
for kv in args.family_override:
    sp, f = kv.split("=")
    fam[sp] = f

best = pd.read_csv(args.candidates, sep="\t")
P, C = {}, {}
for f in args.proteins:
    P.update({r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(f, "fasta")})
for f in args.cds:
    C.update({r.id: str(r.seq).upper() for r in SeqIO.parse(f, "fasta")})

pal = PairwiseAligner()
pal.substitution_matrix = substitution_matrices.load("BLOSUM62")
pal.open_gap_score = -11
pal.extend_gap_score = -1
pal.mode = "global"


def nt(p1, p2, c1, c2):
    """Nucleotide identity (%) and number of aligned codons, from a protein alignment mapped onto the CDSs."""
    a = pal.align(p1, p2)[0]
    m = n = 0
    for (s1, e1), (s2, e2) in zip(*a.aligned):
        for i, j in zip(range(s1, e1), range(s2, e2)):
            x, y = c1[3 * i:3 * i + 3], c2[3 * j:3 * j + 3]
            if len(x) == 3 and len(y) == 3:
                n += 3
                m += sum(p == q for p, q in zip(x, y))
    return (100 * m / n if n else np.nan), n // 3


out = []
for r in best.itertuples():
    if r.q in P and r.s in P and r.q in C and r.s in C:
        v, n = nt(P[r.q], P[r.s], C[r.q], C[r.s])
        out.append((r.q, r.s, v, n))
sc = pd.DataFrame(out, columns=["q", "s", "nt_id", "aligned_codons"])
res = best.merge(sc, on=["q", "s"], how="left")
print(f"scored {res.nt_id.notna().sum()} of {len(res)} candidate genes (missing sequences are left unscored)")

# suspect gene -> HOG, and copies per (HOG, species)
g2h, copies = {}, collections.Counter()
with open(args.n5_tsv) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {h: hdr.index(h) for h in SUSPECTS}
    for l in f:
        r = l.rstrip("\n").split("\t")
        for h, i in ix.items():
            if i < len(r) and r[i]:
                for g in r[i].split(", "):
                    g2h[h + "__" + g] = r[0]
                    copies[(r[0], h)] += 1
res["HOG"] = res.q.map(g2h)

KMAX = args.kmax
strong = (res.nt_id >= 99) & (res.k <= KMAX)
moderate = ((res.nt_id >= 97) & (res.nt_id < 99) & (res.k <= KMAX)) | ((res.nt_id >= 99) & (res.k > KMAX))
res["strength"] = np.where(strong, "strong", np.where(moderate, "moderate", "weak"))

genes = pd.DataFrame({
    "species": res.qsp, "gene_id": res.q.str.split("__").str[1], "HOG": res.HOG,
    "partner_species": res.ssp, "partner_gene": res.s.str.split("__").str[1], "partner_family": res.ssp.map(fam),
    "partner_is_suspect": res.ssp.isin(SUSPECTS), "protein_identity": res.pid.round(2), "aligned_aa": res.len,
    "nt_identity": res.nt_id.round(2), "aligned_codons": res.aligned_codons, "n_species_sharing": res.k,
    "strength": res.strength,
}).sort_values(["strength", "nt_identity"], ascending=[True, False])

os.makedirs(args.out_dir, exist_ok=True)
genes.to_csv(os.path.join(args.out_dir, "flag_table_genes.csv"), index=False)

cell = []
for (hog, sp), d in genes.dropna(subset=["HOG"]).groupby(["HOG", "species"]):
    n = copies[(hog, sp)]
    s = int((d.strength == "strong").sum())
    m = int(d.strength.isin(["strong", "moderate"]).sum())
    cell.append(dict(
        HOG=hog, species=sp, n_copies=n, n_flagged_strong=s, n_flagged_strong_or_moderate=m,
        remaining_if_strong_removed=max(n - s, 0), remaining_if_strong_moderate_removed=max(n - m, 0),
        cell_status_strong=("missing" if s > 0 and n - s <= 0 else ("reduced" if s > 0 else "ok")),
        cell_status_strong_moderate=("missing" if m > 0 and n - m <= 0 else ("reduced" if m > 0 else "ok"))))
cells = pd.DataFrame(cell)
cells.to_csv(os.path.join(args.out_dir, "flag_table_cells.csv"), index=False)

print(f"genes {len(genes)} | strong {int((genes.strength == 'strong').sum())} "
      f"moderate {int((genes.strength == 'moderate').sum())} weak {int((genes.strength == 'weak').sum())}")
st = cells[cells.n_flagged_strong > 0]
print(f"strong cells {len(st)} in {st.HOG.nunique()} HOGs; {int((cells.cell_status_strong == 'missing').sum())} missing")
