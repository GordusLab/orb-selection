"""Step 3: compare cached HyPhy results for HOGs with vs without flagged (suspected contaminated) sequences.

Inputs : results/hyphy_contamination_impact/hog_flag_counts.csv (step 2), results/hyphy_results_cache/*_results.pkl
Outputs: results/hyphy_contamination_impact/comparison_summary.txt and hog_flag_vs_results.csv
Significance as in the paper code: RELAX p<=0.05 (relaxed/intensified); BUSTED-PH 'hit' = test p<=0.05 & background p>0.05 & shared p<=0.05.
Also reports the key BUSTED-PH comparative p (shared_pval), which is the foreground-vs-background test and is what the 'hit' rule hinges on.
Read-only on the caches.  Run with the repo's python (needs scripts/03_selection_tests on the path for the pickled classes).
"""
import sys, os, pickle, io
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
R = os.path.join(HERE, '..', '..', '..', 'results')
OUT = os.path.join(R, 'hyphy_contamination_impact')
cnt = pd.read_csv(os.path.join(OUT, 'hog_flag_counts.csv'), index_col='HOG')

res = {}
for n in ['relax', 'busted_ph_orb', 'busted_ph_non_orb']:
    res[n] = pickle.load(open(os.path.join(R, 'hyphy_results_cache', f'{n}_results.pkl'), 'rb')).results_df

buf = io.StringIO()
def P(*a):
    print(*a); print(*a, file=buf)

groups = {'strong': cnt.ret_strong > 0, 'strong+moderate': cnt.ret_sm > 0, 'any (incl. weak)': cnt.ret_any > 0}
P('HOGs analysed:', len(cnt))
for g, mask in groups.items():
    P(f'  {g}: {mask.sum()} HOGs with >=1 flagged sequence retained ({mask.mean():.1%})')
P('Retained-tip count (n_tips) median: unaffected(strong)', cnt[~groups['strong']].n_tips.median(), ' affected(strong)', cnt[groups['strong']].n_tips.median())

def logit_or(df, ycol, xmask):
    d = pd.DataFrame({'y': df[ycol].astype(int), 'aff': xmask.reindex(df.index).astype(int), 'lntips': np.log(cnt.n_tips.reindex(df.index))})
    X = sm.add_constant(d[['aff', 'lntips']])
    f = sm.Logit(d.y, X).fit(disp=0)
    return np.exp(f.params['aff']), np.exp(f.conf_int().loc['aff']).tolist(), f.pvalues['aff']

rows = []
def block(name, df, sigcol, pcol, effcol=None, efflab=''):
    P(f'\n=== {name} ===')
    df = df.copy(); df['sig'] = df[sigcol]
    P(f'total significant: {df.sig.sum()} of {len(df)} ({df.sig.mean():.1%})')
    for g, mask in groups.items():
        m = mask.reindex(df.index)
        a, u = df[m], df[~m]
        fa, fu = a.sig.mean(), u.sig.mean()
        tab = [[a.sig.sum(), (~a.sig).sum()], [u.sig.sum(), (~u.sig).sum()]]
        fp = stats.fisher_exact(tab)[1]
        orr, ci, lp = logit_or(df, 'sig', mask)
        ks = stats.ks_2samp(a[pcol], u[pcol])
        mw = stats.mannwhitneyu(a[pcol], u[pcol])
        s = (f'{g:18s} affected n={len(a):4d} sig={a.sig.sum():4d} ({fa:.1%}) | unaffected n={len(u):4d} sig={u.sig.sum():4d} ({fu:.1%}) '
             f'| Fisher p={fp:.3g} | adj OR (ln n_tips) {orr:.2f} [{ci[0]:.2f},{ci[1]:.2f}] p={lp:.3g} | median p aff/unaff {a[pcol].median():.3g}/{u[pcol].median():.3g} KS p={ks.pvalue:.3g}')
        if effcol:
            ea, eu = a[effcol], u[effcol]
            s += f' | median {efflab} aff/unaff {ea.median():.3f}/{eu.median():.3f} MW p={stats.mannwhitneyu(ea, eu).pvalue:.3g}'
        P(s)
        rows.append(dict(test=name, group=g, n_aff=len(a), sig_aff=int(a.sig.sum()), n_unaff=len(u), sig_unaff=int(u.sig.sum()), adj_OR=orr, adj_p=lp))
    return df

rx = res['relax'].copy(); rx['sig'] = (rx.p_value <= 0.05) & rx.result.isin(['relaxed', 'intensified']); rx['logk'] = np.log2(rx.k.astype(float))
block('RELAX (p<=0.05)', rx, 'sig', 'p_value', 'logk', 'log2 k')
rx['intens'] = rx.sig & (rx.result == 'intensified'); rx['relaxed'] = rx.sig & (rx.result == 'relaxed')
for c in ['relaxed', 'intensified']:
    rx['_s'] = rx.result.eq(c) & rx.sig
    P(f'\n  -- RELAX {c} only (Test = non-orb, Reference = orb; "intensified" = stronger in non-orb) --')
    for g, mask in groups.items():
        m = mask.reindex(rx.index)
        P(f'  {g:18s} affected {rx[m]._s.mean():.1%} (n={m.sum()})  unaffected {rx[~m]._s.mean():.1%}  adj OR {logit_or(rx, "_s", mask)[0]:.2f} p={logit_or(rx, "_s", mask)[2]:.3g}')

for n, lab in [('busted_ph_orb', 'BUSTED-PH orb fg'), ('busted_ph_non_orb', 'BUSTED-PH non-orb fg')]:
    d = res[n].copy(); d['sig'] = d.result == 'hit'
    block(f'{lab} (hit rule)', d, 'sig', 'shared_pval', 'MG94xREV_ω_test', 'omega_test')
    d['sig_shared'] = d.shared_pval <= 0.05
    block(f'{lab} (comparative p<=0.05 only)', d, 'sig_shared', 'shared_pval')

# flagged side: where do flagged tips sit?
P('\n=== position of retained flagged tips ===')
ft = pd.read_csv(os.path.join(OUT, 'flagged_tips_retained.csv'))
for tier, mr in [('strong', 3), ('strong+moderate', 2), ('any', 1)]:
    x = ft[ft['rank'] >= mr]
    P(f'{tier}: {len(x)} tips in {x.HOG.nunique()} HOGs; orb-weaver tips (orb-fg, RELAX Reference): {x.fg_orb.sum()}; non-orb tips (non-orb-fg, RELAX Test): {x.fg_nonorb.sum()}; partner copy also retained: {x.partner_retained.sum()}')
    P('   by species:', x.species.value_counts().to_dict())
x = ft[ft['rank'] == 3]
P('strong flagged tips: partner in another family retained too in', x.partner_retained.mean().round(3), 'of cases')

# affected hit lists
P('\n=== hit lists that feed GO enrichment: how many hits are affected HOGs ===')
hit_sets = {
    'RELAX relaxed (k<1 in non-orb Test)': rx.index[rx.sig & (rx.result == 'relaxed')],
    'RELAX intensified': rx.index[rx.sig & (rx.result == 'intensified')],
    'BUSTED-PH orb fg hits': res['busted_ph_orb'].index[res['busted_ph_orb'].result == 'hit'],
    'BUSTED-PH non-orb fg hits': res['busted_ph_non_orb'].index[res['busted_ph_non_orb'].result == 'hit'],
}
out = []
for k, ix in hit_sets.items():
    ix = pd.Index(ix)
    s = {g: int(mask.reindex(ix).sum()) for g, mask in groups.items()}
    exp = {g: mask.mean() * len(ix) for g, mask in groups.items()}
    P(f'{k}: {len(ix)} HOGs; affected strong {s["strong"]} (expected by chance {exp["strong"]:.1f}), +moderate {s["strong+moderate"]}, any {s["any (incl. weak)"]} (exp {exp["any (incl. weak)"]:.1f})')
    out.append((k, ix))
ba = pd.Index(res['busted_ph_orb'].index[res['busted_ph_orb'].result == 'hit'])
P('BUSTED-PH orb hit HOGs:', list(ba), cnt.loc[ba, ['n_tips', 'ret_strong', 'ret_sm', 'ret_any', 'n_flag_table_genes']].to_dict('index'))
for h in ['N5.HOG0052434', 'N5.HOG0068277']:
    if h in cnt.index:
        P(' raw-p candidate', h, cnt.loc[h, ['n_tips', 'ret_strong', 'ret_sm', 'ret_any']].to_dict(), 'shared_pval', res['busted_ph_orb'].loc[h, 'shared_pval'] if h in res['busted_ph_orb'].index else None)

j = cnt.join(rx[['p_value', 'k', 'result', 'sig']].add_prefix('relax_')).join(res['busted_ph_orb'][['shared_pval', 'result']].add_prefix('borb_')).join(res['busted_ph_non_orb'][['shared_pval', 'result']].add_prefix('bnonorb_'))
j.to_csv(os.path.join(OUT, 'hog_flag_vs_results.csv'))
pd.DataFrame(rows).to_csv(os.path.join(OUT, 'comparison_table.csv'), index=False)
open(os.path.join(OUT, 'comparison_summary.txt'), 'w').write(buf.getvalue())
