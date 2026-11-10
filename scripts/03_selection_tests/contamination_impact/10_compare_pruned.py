"""Step 10: compare original cached HyPhy results with the pruned re-run for the 359 HOGs.
Inputs : results/hyphy_contamination_impact/pruned_summary.tsv (from 09), hog_flag_vs_results.csv (step 3), hog_flag_counts.csv
Output : pruned_comparison_summary.txt, pruned_vs_original.csv
RELAX significance: p<=0.05 (relaxed/intensified). BUSTED-PH orb: 'hit' rule as in the paper code.
"""
import os, numpy as np, pandas as pd
from scipy import stats
O = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'results', 'hyphy_contamination_impact')
pr = pd.read_csv(os.path.join(O, 'pruned_summary.tsv'), sep='\t')
orig = pd.read_csv(os.path.join(O, 'hog_flag_vs_results.csv'), index_col='HOG')
buf = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); buf.append(s)
rx = pr[pr.test == 'relax'].set_index('HOG'); bo = pr[pr.test == 'busted_orb'].set_index('HOG')
P('pruned results available: RELAX', len(rx), ' BUSTED-PH orb', len(bo))
P('pruned BUSTED result values:', bo.result.value_counts().to_dict())

d = orig.join(rx[['n_seq', 'p', 'k', 'result']].add_prefix('pr_'), how='inner')
d['pr_sig'] = d.pr_result.isin(['relaxed', 'intensified']) & (d.pr_p <= 0.05)
d['o_sig'] = d.relax_sig.astype(bool)
d['dseq'] = d.n_tips - d.pr_n_seq
P('\n=== RELAX ===  n =', len(d), '; tips removed per HOG: median', d.dseq.median(), 'range', d.dseq.min(), '-', d.dseq.max())
P('significant original', d.o_sig.sum(), 'pruned', d.pr_sig.sum())
P(pd.crosstab(d.relax_result.where(d.o_sig, 'ns'), d.pr_result.where(d.pr_sig, 'ns')).to_string())
P('lost significance:', (d.o_sig & ~d.pr_sig).sum(), ' gained:', (~d.o_sig & d.pr_sig).sum(), ' changed direction:', (d.o_sig & d.pr_sig & (d.relax_result != d.pr_result)).sum())
lk, pk = np.log2(d.relax_k.astype(float)), np.log2(d.pr_k.astype(float))
ok = np.isfinite(lk) & np.isfinite(pk)
P('log2 k original vs pruned: Spearman', round(stats.spearmanr(lk[ok], pk[ok])[0], 3), ' median |diff|', round((lk - pk)[ok].abs().median(), 3))
lp, pp = -np.log10(d.relax_p_value), -np.log10(d.pr_p)
P('-log10 p: Spearman', round(stats.spearmanr(lp, pp)[0], 3), ' median |diff|', round((lp - pp).abs().median(), 3))
# baseline expectation: how often does an unrelated change flip significance near threshold?
near = d[(d.relax_p_value > 0.02) & (d.relax_p_value < 0.1)]
P('HOGs with original p in (0.02,0.1):', len(near), ' of which sig status changed:', (near.o_sig != near.pr_sig).sum())
big = d[(d.o_sig != d.pr_sig)]
P('status-changed HOGs: original p/k -> pruned p/k')
P(big[['relax_p_value', 'relax_k', 'relax_result', 'pr_p', 'pr_k', 'pr_result', 'ret_strong', 'ret_sm', 'dseq']].round(4).to_string())

b = orig.join(bo[['n_seq', 'p', 'shared_p', 'background_p', 'result']].add_prefix('pr_'), how='inner')
b['o_hit'] = b.borb_result == 'hit'; b['pr_hit'] = b.pr_result == 'hit'
b['o_sh'] = b.borb_shared_pval <= 0.05; b['pr_sh'] = b.pr_shared_p <= 0.05
P('\n=== BUSTED-PH orb fg ===  n =', len(b))
P('hit original', b.o_hit.sum(), 'pruned', b.pr_hit.sum(), '; comparative p<=0.05 original', b.o_sh.sum(), 'pruned', b.pr_sh.sum(),
  '; lost', (b.o_sh & ~b.pr_sh).sum(), 'gained', (~b.o_sh & b.pr_sh).sum())
s = b.dropna(subset=['pr_shared_p'])
P('comparative p: Spearman', round(stats.spearmanr(s.borb_shared_pval, s.pr_shared_p)[0], 3))
P(b[b.o_hit | b.pr_hit | (b.o_sh != b.pr_sh)][['borb_shared_pval', 'pr_shared_p', 'borb_result', 'pr_result', 'ret_strong', 'dseq' if 'dseq' in b else 'n_tips']].round(4).to_string())

# does effect size of change depend on mixed-pair / strong count?
d['absdlogk'] = (lk - pk).abs()
P('\nmedian |dlog2k| by strong-tips removed (ret_strong): ', d.groupby(d.ret_strong.clip(upper=3)).absdlogk.median().round(3).to_dict())
d.to_csv(os.path.join(O, 'pruned_vs_original.csv'))
open(os.path.join(O, 'pruned_comparison_summary.txt'), 'w').write('\n'.join(buf) + '\n')
