"""Step 1: list retained tips (after PREQUAL/MACSE/filter/ClipKIT/dedup) and their HyPhy branch label per test.

Reads the cached HyPhy JSONs (results/hyphy_results_cache/{relax,busted_ph_orb,busted_ph_non_orb}/) and writes one long
table: HOG, test, tip_name (as in the JSON), species, gene_id, label.
Labels: RELAX 'Test'/'Reference'; BUSTED-PH 'test'/'background' (test = foreground).
Read-only on the JSONs; output is a TSV.
usage: 01_extract_tips.py OUT.tsv [N_PROCS]
"""
import json, sys, glob, os
from multiprocessing import Pool

ROOT = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'results', 'hyphy_results_cache')
TESTS = {'relax': 'relax', 'busted_ph_orb': 'busted_orb', 'busted_ph_non_orb': 'busted_nonorb'}

def one(path):
    test = TESTS[os.path.basename(os.path.dirname(path))]
    hog = os.path.basename(path).split('_')[0]
    j = json.load(open(path))
    lab = j['tested']['0']
    ba = j['branch attributes']['0']
    rows = []
    for tip, attrs in ba.items():
        if tip.startswith('Node'):
            continue
        sp, gid = attrs['original name'].split('|', 1)
        rows.append((hog, test, tip, sp, gid, lab.get(tip, 'NA')))
    return rows, j['input']['number of sequences']

if __name__ == '__main__':
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    files = sorted(glob.glob(os.path.join(ROOT, '*', 'N5.HOG*.json')))
    files = [f for f in files if os.path.basename(os.path.dirname(f)) in TESTS]
    print(len(files), 'files', flush=True)
    with Pool(n) as p, open(out, 'w') as o:
        o.write('HOG\ttest\ttip_name\tspecies\tgene_id\tlabel\tn_seq_json\n')
        for i, (rows, ns) in enumerate(p.imap_unordered(one, files, chunksize=4)):
            for r in rows:
                o.write('\t'.join(map(str, r)) + f'\t{ns}\n')
            if i % 1000 == 0:
                print(i, flush=True)
