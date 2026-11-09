"""Step 2: per HOG, count flagged (suspected cross-sample contamination) sequences among the tips HyPhy actually used.

Inputs : tips_long.tsv from 01_extract_tips.py; data/flag_table_genes.csv; data/N5.udiv.o75_list.txt
Outputs: <outdir>/hog_flag_counts.csv  (one row per HOG; counts for each test, fg/bg)
         <outdir>/flagged_tips_retained.csv (one row per flagged gene retained in a HyPhy alignment)
Tiers are cumulative: 'strong' = strong; 'sm' = strong+moderate; 'any' = strong+moderate+weak.
Partner side: flag-table partner_gene (the near-identical copy in another family), counted when it is retained and is not itself
flagged ('partner_only'); 'pair_retained' = a flagged gene AND its partner are both retained in the same HOG alignment.
usage: 02_count_flagged.py tips_long.tsv OUTDIR
"""
import sys, os
import pandas as pd

tips_f, outdir = sys.argv[1], sys.argv[2]
D = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'data')
os.makedirs(outdir, exist_ok=True)

hogs = [l.strip() for l in open(os.path.join(D, 'N5.udiv.o75_list.txt')) if l.strip()]
tips = pd.read_csv(tips_f, sep='\t')
tips['gid'] = tips.gene_id.str.upper()
tips['sp'] = tips.species.str.upper()

fl = pd.read_csv(os.path.join(D, 'flag_table_genes.csv'))
fl['sp'] = fl.species.str.upper(); fl['gid'] = fl.gene_id.str.upper()
fl['psp'] = fl.partner_species.str.upper(); fl['pgid'] = fl.partner_gene.str.upper()
rank = {'weak': 1, 'moderate': 2, 'strong': 3}
fl['rank'] = fl.strength.map(rank)

# which HOG ids are in the 4,756 list
assert set(tips.HOG) == set(hogs), (len(set(tips.HOG)), len(hogs))

# retained flagged genes
key = tips[tips.test == 'relax'][['HOG', 'sp', 'gid', 'tip_name']]
m = fl.merge(key, left_on=['HOG', 'sp', 'gid'], right_on=['HOG', 'sp', 'gid'], how='left', indicator=True)
m['retained'] = m._merge == 'both'
print('flag-table genes in list HOGs:', m.HOG.isin(hogs).sum(), 'of', len(m), '; retained in HyPhy alignments:', m.retained.sum())

# labels per tip per test
lab = tips.pivot_table(index=['HOG', 'gid'], columns='test', values='label', aggfunc='first').reset_index()
ret = m[m.retained].merge(lab, on=['HOG', 'gid'], how='left')
ret['fg_orb'] = ret.busted_orb == 'test'          # orb-weaver foreground in BUSTED-PH orb
ret['fg_nonorb'] = ret.busted_nonorb == 'test'    # non-orb foreground in BUSTED-PH non-orb
ret['relax_ref_orb'] = ret.relax == 'Reference'   # RELAX Reference = orb-weavers, Test = non-orb
# partner retained?
pk = set(zip(tips[tips.test == 'relax'].HOG, tips[tips.test == 'relax'].gid))
ret['partner_retained'] = [(h, g) in pk for h, g in zip(ret.HOG, ret.pgid)]
flagged_keys = set(zip(fl.HOG, fl.gid))
ret.to_csv(os.path.join(outdir, 'flagged_tips_retained.csv'), index=False)

# per-HOG table
base = tips[tips.test == 'relax'].groupby('HOG').agg(n_tips=('gid', 'size'))
base['n_orb_tips'] = tips[tips.test == 'busted_orb'].query("label=='test'").groupby('HOG').size()
base = base.reindex(hogs).fillna(0).astype(int)
tab = base.copy()
tab['n_flag_table_genes'] = m[m.HOG.isin(hogs)].groupby('HOG').size().reindex(hogs).fillna(0).astype(int)
for tier, minrank in [('strong', 3), ('sm', 2), ('any', 1)]:
    r = ret[ret['rank'] >= minrank]
    g = r.groupby('HOG')
    tab[f'ret_{tier}'] = g.size().reindex(hogs).fillna(0).astype(int)
    tab[f'ret_{tier}_orbfg'] = g.fg_orb.sum().reindex(hogs).fillna(0).astype(int)      # on orb-fg branches (orb-weaver tips)
    tab[f'ret_{tier}_nonorbfg'] = g.fg_nonorb.sum().reindex(hogs).fillna(0).astype(int)  # on non-orb-fg tips
    tab[f'ret_{tier}_pair'] = g.partner_retained.sum().reindex(hogs).fillna(0).astype(int)
tab['ret_strong_species'] = ret[ret['rank'] == 3].groupby('HOG').species.agg(lambda s: ';'.join(sorted(set(s)))).reindex(hogs).fillna('')
tab['affected_strong'] = tab.ret_strong > 0
tab['affected_sm'] = tab.ret_sm > 0
tab['affected_any'] = tab.ret_any > 0
tab.index.name = 'HOG'
tab.to_csv(os.path.join(outdir, 'hog_flag_counts.csv'))
print(tab[['affected_strong', 'affected_sm', 'affected_any']].sum())
print(tab[['ret_strong', 'ret_sm', 'ret_any']].describe().T)
