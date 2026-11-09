#!/usr/bin/env python3
"""Remove named tips from a HOG's codon alignment (NEXUS, as written by remove-duplicates) and from its IQ-TREE treefile.

The alignment is NOT re-trimmed or re-aligned; the tree is the existing IQ-TREE tree with the tips pruned (single-child nodes
collapsed, branch lengths summed). HyPhy re-estimates branch lengths. Re-label the pruned tree afterwards with label-tree.bf.
Fails loudly if a requested tip is missing, or if fewer than --min-side tips remain.
usage: 06_prune_hog.py --nex IN.nex --tree IN.treefile --remove 'sp|id;sp|id' --out-nex OUT.nex --out-tree OUT.tree
"""
import argparse, re, sys, io
from Bio import Phylo

ap = argparse.ArgumentParser()
ap.add_argument('--nex', required=True); ap.add_argument('--tree', required=True)
ap.add_argument('--remove', required=True); ap.add_argument('--out-nex', required=True); ap.add_argument('--out-tree', required=True)
ap.add_argument('--min-tips', type=int, default=4)
a = ap.parse_args()
rm = [t for t in a.remove.split(';') if t]

txt = open(a.nex).read()
m = re.search(r"TAXLABELS\s*(.*?)\s*;\s*END;", txt, re.S)
labels = re.findall(r"'([^']*)'", m.group(1))
mm = re.search(r"MATRIX\s*\n(.*?)\n\s*;?\s*END;", txt, re.S)
seqs = [l.strip().rstrip(';') for l in mm.group(1).split('\n') if l.strip().rstrip(';')]
assert len(seqs) == len(labels), (len(seqs), len(labels))
missing = [t for t in rm if t not in labels]
if missing:
    sys.exit(f'ERROR: tips not in alignment: {missing}')
keep = [i for i, l in enumerate(labels) if l not in set(rm)]
if len(keep) < a.min_tips:
    sys.exit('ERROR: too few tips left')
nchar = len(seqs[0])
assert all(len(s) == nchar for s in seqs)
head = txt[:m.start()].replace(re.search(r"NTAX = \d+", txt).group(0), f"NTAX = {len(keep)}")
new = (head + "TAXLABELS\n\t\t" + ' '.join(f"'{labels[i]}'" for i in keep) + "\n;\nEND;\n")
rest = txt[m.end():]
# rebuild everything after the taxa block: characters block header up to MATRIX unchanged
pre = re.search(r"(BEGIN CHARACTERS;.*?MATRIX\s*\n)", txt, re.S).group(1)
new = head + "TAXLABELS\n\t\t" + ' '.join(f"'{labels[i]}'" for i in keep) + "\n;\nEND;\n\n" + pre + '\n'.join(seqs[i] for i in keep) + "\n;\nEND;\n"
open(a.out_nex, 'w').write(new)

tr = Phylo.read(a.tree, 'newick')
tips = {t.name for t in tr.get_terminals()}
assert tips == set(labels), 'tree tips differ from alignment taxa'
for t in rm:
    tr.prune(t)
assert {t.name for t in tr.get_terminals()} == {labels[i] for i in keep}
Phylo.write(tr, a.out_tree, 'newick', format_branch_length='%.10f')
print(f'kept {len(keep)} of {len(labels)} tips; removed {len(rm)}')
