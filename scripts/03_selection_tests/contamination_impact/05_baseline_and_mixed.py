"""Step 5: (a) baseline for close mixed-label pairs: same measurement as step 4 for UNFLAGGED tips of the nine suspect species in a random
sample of HOGs with no flagged sequence; (b) HyPhy results for HOGs with a close mixed pair vs other affected HOGs vs unaffected.
usage: 05_baseline_and_mixed.py tips_long.tsv OUTDIR
"""
import sys, os, json, io, random
import numpy as np, pandas as pd
from scipy import stats
from Bio import Phylo
tips_f, OUT = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, '..', '..', '..', 'data'); R = os.path.join(HERE, '..', '..', '..', 'results')
buf = io.StringIO()
def P(*a):
    print(*a); print(*a, file=buf)
tips = pd.read_csv(tips_f, sep='\t'); orb = tips[tips.test == 'busted_orb']
lab = dict(zip(orb.tip_name + '|' + orb.HOG, orb.label))
suspects = [l.strip().upper() for l in open(os.path.join(D, 'hub_species_to_exclude.txt')) if l.strip()]
cnt = pd.read_csv(os.path.join(OUT, 'hog_flag_counts.csv'), index_col='HOG')
clean = [h for h in cnt.index[cnt.ret_any == 0]]
random.seed(1); samp = random.sample(clean, 359)
rows = []
for hog in samp:
    sub = orb[(orb.HOG == hog) & orb.species.str.upper().isin(suspects)]
    if sub.empty: continue
    j = json.load(open(os.path.join(R, 'hyphy_results_cache', 'relax', f'{hog}_RELAX.json')))
    tr = Phylo.read(io.StringIO(j['input']['trees']['0']), 'newick')
    terms = {c.name: c for c in tr.get_terminals()}
    for me in sub.tip_name:
        d = {n: tr.distance(terms[me], terms[n]) for n in terms if n != me}
        nn = min(d, key=d.get)
        rows.append((hog, me, d[nn], lab[me + '|' + hog] != lab[nn + '|' + hog]))
b = pd.DataFrame(rows, columns=['HOG', 'tip', 'nn_dist', 'mixed'])
f = pd.read_csv(os.path.join(OUT, 'flagged_tip_neighbours.csv'))
P('=== (a) baseline: unflagged suspect-species tips in 359 random HOGs with no flagged sequence (n=%d tips) vs flagged strong+moderate tips (n=%d) ===' % (len(b), len(f)))
for nm, d in [('baseline unflagged', b), ('flagged strong+moderate', f)]:
    dd = d.nn_dist; mx = d['mixed'] if 'mixed' in d else d.mixed_nn
    P(f'{nm:26s} nn<0.01: {(dd<0.01).mean():.1%} | mixed nn: {mx.mean():.1%} | close (<0.01) AND mixed: {((dd<0.01)&mx).mean():.1%} | median nn dist {dd.median():.4f}')

# (b) results by pair status
res = pd.read_csv(os.path.join(OUT, 'hog_flag_vs_results.csv'), index_col='HOG')
mixed_h = set(f[(f.nn_dist < 0.01) & f.mixed_nn].HOG)
aff_h = set(cnt.index[cnt.ret_sm > 0])
grp = pd.Series('unaffected (no strong/moderate)', index=res.index)
grp[res.index.isin(aff_h)] = 'affected, no close mixed pair'
grp[res.index.isin(mixed_h)] = 'affected, close mixed pair'
res['grp'] = grp
res['relax_sig_b'] = res.relax_sig.astype(bool)
res['borb_sig'] = res.borb_shared_pval <= 0.05; res['bnon_sig'] = res.bnonorb_shared_pval <= 0.05
P('\n=== (b) HyPhy outcomes by group (strong+moderate flags) ===')
t = res.groupby('grp').agg(n=('relax_sig_b', 'size'), relax_sig=('relax_sig_b', 'mean'), relax_med_p=('relax_p_value', 'median'),
    borb_cmp_sig=('borb_sig', 'mean'), borb_med_p=('borb_shared_pval', 'median'), bnon_cmp_sig=('bnon_sig', 'mean'), bnon_med_p=('bnonorb_shared_pval', 'median'), n_tips=('n_tips', 'median'))
P(t.round(3).to_string())
for c, pc in [('relax_sig_b', 'relax_p_value'), ('borb_sig', 'borb_shared_pval'), ('bnon_sig', 'bnonorb_shared_pval')]:
    a = res[res.grp == 'affected, close mixed pair']; u = res[res.grp == 'unaffected (no strong/moderate)']
    P(f'{c}: close-mixed vs unaffected Fisher p={stats.fisher_exact([[a[c].sum(), (~a[c]).sum()], [u[c].sum(), (~u[c]).sum()]])[1]:.3g}; p-value KS p={stats.ks_2samp(a[pc], u[pc]).pvalue:.3g}')
open(os.path.join(OUT, 'baseline_and_mixed_summary.txt'), 'w').write(buf.getvalue())
