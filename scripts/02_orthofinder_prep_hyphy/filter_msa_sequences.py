#!/usr/bin/env python3
"""Remove truncated / frameshifted sequences from a MACSE codon alignment.

A sequence is removed if more than --max-frameshift-frac of its residues are frameshift
marks ("!"; isolated frameshifts are kept and masked later), it has more than
--max-internal-stops internal stops, or an ungapped length below --min-frac of the gene's median ungapped length.
Protected taxa (default Uloborus_diversus) are only removed if they have a frameshift/stop
or fall below the lower --protected-min-frac.

Safety nets (the gene is then written unfiltered and the reason is logged):
  - removal would exceed --max-removed-frac of the sequences
  - fewer than --min-seqs sequences would remain
  - all orbweavers or all non-orbweavers would be removed
Columns that are all gaps after removal are dropped.

Usage:
  filter_msa_sequences.py <in.fasta> <out.fasta> <removed.tsv> <orbweavers-list.txt> [options]
"""
import argparse
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_macse_sequences import internal_stops, read_fasta  # noqa: E402


def write_fasta(path, names, seqs):
    with open(path, "w") as fh:
        for n, s in zip(names, seqs):
            fh.write(f">{n}\n{s}\n")


def mask_frameshifts(seq):
    """Replace every codon containing a MACSE "!" with a gap codon."""
    if "!" not in seq:
        return seq
    return "".join("---" if "!" in seq[i:i + 3] else seq[i:i + 3] for i in range(0, len(seq), 3))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("in_fasta")
    ap.add_argument("out_fasta")
    ap.add_argument("removed_tsv")
    ap.add_argument("orb_list")
    ap.add_argument("--min-frac", type=float, default=0.5)
    ap.add_argument("--protected", default="Uloborus_diversus", help="species given the lower length floor")
    ap.add_argument("--protected-min-frac", type=float, default=0.3)
    ap.add_argument("--max-frameshift-frac", type=float, default=0.05, help="max fraction of '!' marks per ungapped residue")
    ap.add_argument("--max-internal-stops", type=int, default=0)
    ap.add_argument("--max-removed-frac", type=float, default=0.3)
    ap.add_argument("--min-seqs", type=int, default=4)
    a = ap.parse_args()

    with open(a.orb_list) as fh:
        orbs = {l.strip().replace("\r", "") for l in fh if l.strip()}

    names, seqs = read_fasta(a.in_fasta)
    lens = [sum(1 for ch in s if ch not in "-!") for s in seqs]
    med = statistics.median(lens)

    reasons = []
    for name, s, n in zip(names, seqs, lens):
        frac = n / med if med else 0.0
        floor = a.protected_min_frac if name.split("|")[0] == a.protected else a.min_frac
        why = []
        if frac < floor:
            why.append(f"short({frac:.2f}<{floor})")
        # "!" marks come in triples; the later "!" cleanup masks the affected codons
        fs_frac = s.count("!") / max(n, 1)
        if fs_frac > a.max_frameshift_frac:
            why.append(f"frameshift({fs_frac:.3f})")
        if internal_stops(s) > a.max_internal_stops:
            why.append("internal_stop")
        reasons.append(",".join(why))

    remove = [i for i, r in enumerate(reasons) if r]
    keep = [i for i in range(len(seqs)) if not reasons[i]]
    species = [n.split("|")[0] for n in names]

    skip = ""
    if len(remove) > a.max_removed_frac * len(seqs):
        skip = f"would remove {len(remove)}/{len(seqs)} (> {a.max_removed_frac})"
    elif len(keep) < a.min_seqs:
        skip = f"only {len(keep)} sequences would remain"
    elif not any(species[i] in orbs for i in keep) and any(s in orbs for s in species):
        skip = "no orbweavers would remain"
    elif not any(species[i] not in orbs for i in keep) and any(s not in orbs for s in species):
        skip = "no non-orbweavers would remain"

    with open(a.removed_tsv, "w") as fh:
        fh.write("sequence\tungapped_len\tfrac_of_median\treason\tstatus\n")
        status = f"NOT_REMOVED:{skip}" if skip else "removed"
        for i in remove:
            frac = lens[i] / med if med else 0.0
            fh.write(f"{names[i]}\t{lens[i]}\t{frac:.3f}\t{reasons[i]}\t{status}\n")

    # HyPhy can't read "!"; masked codons then count as gaps for trimming
    seqs = [mask_frameshifts(s) for s in seqs]

    if skip:
        write_fasta(a.out_fasta, names, seqs)
        print(f"Filter skipped ({skip}); wrote all {len(seqs)} sequences.")
        return

    kept_names = [names[i] for i in keep]
    kept_seqs = [seqs[i] for i in keep]
    # Drop codons that are all gaps in the remaining sequences
    L = len(kept_seqs[0])
    cols = []
    for c in range(0, L - L % 3, 3):
        if any(s[c:c + 3] != "---" for s in kept_seqs):
            cols.extend((c, c + 1, c + 2))
    kept_seqs = ["".join(s[c] for c in cols) for s in kept_seqs]
    write_fasta(a.out_fasta, kept_names, kept_seqs)
    print(f"Filter removed {len(remove)}/{len(seqs)} sequences; "
          f"alignment {len(seqs[0])} -> {len(cols)} columns.")


if __name__ == "__main__":
    main()
