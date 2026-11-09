#!/usr/bin/env python3
"""Run on the cluster: summarise the pruned HyPhy JSONs (one row per HOG x test) into a small TSV.
Classification mirrors hyphy_results_parser.py: RELAX p<=0.05 (k<1 relaxed, k>1 intensified); BUSTED-PH hit = test p<=0.05 & background p>0.05 & comparative p<=0.05.
usage: 09_extract_pruned.py OUT_DIR OUT.tsv   (OUT_DIR = hyphy_prune_261109/out)
"""
import sys, os, json, glob
out_dir, out = sys.argv[1], sys.argv[2]
def first(d, *keys):
    for k in keys:
        if k in d: return d[k]
rows = []
n_incomplete = 0
for f in sorted(glob.glob(os.path.join(out_dir, 'N5.HOG*', '*.pruned.json'))):
    if f.endswith('_intermediate.json'): continue
    hog = os.path.basename(f).split('_')[0]
    raw = open(f, 'rb').read() if os.path.getsize(f) else b''
    # HyPhy writes the JSON early (BUSTED-PH: after the plain BUSTED fit, before the comparative/background tests). Count only finished runs.
    need = b'relaxation or intensification parameter' if 'RELAX' in os.path.basename(f) else b'Comparative selection test results'
    if need not in raw or b'"p-value"' not in raw:
        n_incomplete += 1; continue
    try:
        j = json.load(open(f))
        n = j['input']['number of sequences']
        if 'RELAX' in os.path.basename(f):
            r = j['test results']; p = r['p-value']; k = r['relaxation or intensification parameter']
            res = 'ns' if p > 0.05 else ('relaxed' if k < 1 else 'intensified')
            rows.append((hog, 'relax', n, p, k, '', '', res))
        else:
            test = 'busted_orb' if 'orb_fg' in os.path.basename(f) and 'non_orb' not in f else 'busted_nonorb'
            tp = j['test results']['p-value']
            bg = first(j, 'Background selection test results', 'test results background')
            cm = first(j, 'Comparative selection test results', 'test results shared distributions')
            bp = bg['p-value'] if bg else float('nan'); sp = cm['p-value'] if cm else float('nan')  # absent when the unconstrained fit shows no selection (test p = 0.5)
            res = 'hit' if (tp <= .05 and bp > .05 and sp <= .05) else 'ns'
            rows.append((hog, test, n, tp, '', bp, sp, res))
    except Exception as e:
        rows.append((hog, os.path.basename(f), '', '', '', '', '', 'ERR:' + type(e).__name__))
with open(out, 'w') as o:
    o.write('HOG\ttest\tn_seq\tp\tk\tbackground_p\tshared_p\tresult\n')
    for r in rows: o.write('\t'.join(map(str, r)) + '\n')
print(len(rows), 'rows;', n_incomplete, 'incomplete files skipped; errors:', sum(r[-1].startswith('ERR') for r in rows))
