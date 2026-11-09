"""Step 4: mechanisms. (a) retention of flagged vs unflagged genes of the nine suspect species through the stage-02 filtering/dedup;
(b) for flagged tips that are retained: nearest neighbour in the HyPhy input tree, whether the near-identical partner copy is
retained, and whether the pair straddles the foreground/background label (orb vs non-orb).
Inputs: tips_long.tsv (step 1), results/hyphy_contamination_impact/flagged_tips_retained.csv (step 2), data/N5.tsv, flag table, cached RELAX JSONs
        (input tree = the IQ-TREE tree HyPhy started from; patristic distance computed on that tree).
usage: 04_pairs_and_retention.py tips_long.tsv OUTDIR
"""
import sys, os, json, io
import numpy as np, pandas as pd
from Bio import Phylo

tips_f, OUT = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, '..', '..', '..', 'data'); R = os.path.join(HERE, '..', '..', '..', 'results')
buf = io.StringIO()
def P(*a):
    print(*a); print(*a, file=buf)

tips = pd.read_csv(tips_f, sep='\t'); tips['gid'] = tips.gene_id.str.upper(); tips['sp'] = tips.species.str.upper()
orbt = tips[tips.test == 'busted_orb'][['HOG', 'gid', 'tip_name', 'label']]
fl = pd.read_csv(os.path.join(D, 'flag_table_genes.csv')); fl['gid'] = fl.gene_id.str.upper(); fl['pgid'] = fl.partner_gene.str.upper()
hogs = [l.strip() for l in open(os.path.join(D, 'N5.udiv.o75_list.txt')) if l.strip()]
suspects = [l.strip() for l in open(os.path.join(D, 'hub_species_to_exclude.txt')) if l.strip()]
retained = set(zip(orbt.HOG, orbt.gid))

# (a) retention through filtering for all genes of the nine species in the 4,756 HOGs
n5 = pd.read_csv(os.path.join(D, 'N5.tsv'), sep='\t', usecols=['HOG'] + suspects, dtype=str, low_memory=False).set_index('HOG').reindex(hogs)
rows = []
for sp in suspects:
    for hog, cell in n5[sp].dropna().items():
        for g in str(cell).split(', '):
            g = g.strip()
            if g:
                rows.append((hog, sp, g.upper()))
allg = pd.DataFrame(rows, columns=['HOG', 'species', 'gid'])
allg['retained'] = [(h, g) in retained for h, g in zip(allg.HOG, allg.gid)]
fk = fl.set_index(['HOG', 'gid']).strength
allg['strength'] = [fk.get((h, g), 'unflagged') for h, g in zip(allg.HOG, allg.gid)]
P('=== (a) retention of nine suspect species\' genes in the 4,756 HOGs (in OrthoFinder HOG -> in HyPhy alignment) ===')
t = allg.groupby('strength').retained.agg(['size', 'sum', 'mean']).rename(columns={'size': 'n_genes', 'sum': 'n_retained', 'mean': 'frac_retained'})
P(t.to_string())
P('(a2) flagged genes by nt identity: retained fraction')
x = fl.merge(allg[['HOG', 'gid', 'retained']], on=['HOG', 'gid'])
x['ident_bin'] = pd.cut(x.nt_identity, [0, 97, 99, 99.99, 100.0], labels=['<97', '97-99', '99-<100', '100'])
x['partner_retained'] = [(h, g) in retained for h, g in zip(x.HOG, x.pgid)]
P(x.groupby('ident_bin', observed=True).agg(n=('retained', 'size'), frac_retained=('retained', 'mean'), frac_partner_retained=('partner_retained', 'mean')).to_string())
both = x[x.retained & x.partner_retained]
P('(a3) among flagged genes retained AND partner retained, nt identity bins:', both.ident_bin.value_counts().sort_index().to_dict())
P('     flagged retained, partner (same HOG listing) retained, by ident_bin: both-retained fraction of 100%-identical pairs =',
  round(((x.ident_bin == '100') & x.retained & x.partner_retained).sum() / max((x.ident_bin == '100').sum(), 1), 3))
# of identical (100%) pairs whose two members are both present in the HOG listing: how often are both kept?
pk = set(zip(allg.HOG, allg.gid))
# partner may be a non-suspect species, so check partner presence in the HOG via tips of ALL HyPhy tips or the N5 listing is unavailable here:
P('     (partner presence in the pre-filter HOG is not checked for non-suspect partners; see (b) for retained partners)')

# (b) pairs in the input tree
ft = pd.read_csv(os.path.join(OUT, 'flagged_tips_retained.csv'))
ft = ft[ft['rank'] >= 2]  # strong + moderate
lab_orb = {(h, g): l for h, g, l in zip(orbt.HOG, orbt.gid, orbt.label)}
name = {(h, g): n for h, g, n in zip(orbt.HOG, orbt.gid, orbt.tip_name)}
out = []
for hog, sub in ft.groupby('HOG'):
    j = json.load(open(os.path.join(R, 'hyphy_results_cache', 'relax', f'{hog}_RELAX.json')))
    tr = Phylo.read(io.StringIO(j['input']['trees']['0']), 'newick')
    terms = {c.name: c for c in tr.get_terminals()}
    names = list(terms)
    for _, r in sub.iterrows():
        me = name[(hog, r.gid)]
        # nearest neighbour by patristic distance
        d = {n: tr.distance(terms[me], terms[n]) for n in names if n != me}
        nn = min(d, key=d.get)
        row = dict(HOG=hog, species=r.species, gid=r.gid, strength=r.strength, nt_identity=r.nt_identity, tip_len=terms[me].branch_length,
                   nn=nn, nn_dist=d[nn], my_label=lab_orb[(hog, r.gid)], partner_species=r.partner_species,
                   partner_retained=r.partner_retained)
        if r.partner_retained:
            pn = name[(hog, r.pgid)]
            row['partner_dist'] = d[pn]; row['partner_label'] = lab_orb[(hog, r.pgid)]
            row['partner_is_nn'] = (pn == nn)
        # label of nn
        nnid = [g for (h, g), n in name.items() if h == hog and n == nn] if False else None
        out.append(row)
res = pd.DataFrame(out)
# label of nearest neighbour
nlab = dict(zip(orbt.tip_name + '|' + orbt.HOG, orbt.label))
res['nn_label'] = [nlab.get(n + '|' + h) for n, h in zip(res.nn, res.HOG)]
res['mixed_nn'] = res.my_label != res.nn_label
res.to_csv(os.path.join(OUT, 'flagged_tip_neighbours.csv'), index=False)
P(f'\n=== (b) strong+moderate flagged tips retained: {len(res)} in {res.HOG.nunique()} HOGs ===')
P('tip branch length (subst/codon, IQ-TREE tree): median %.4f, 90th pct %.4f; all-tips median in these trees is larger or similar (not computed)' % (res.tip_len.median(), res.tip_len.quantile(.9)))
P('nearest-neighbour distance: median %.4f; <0.01: %d (%.0f%%); <0.001: %d' % (res.nn_dist.median(), (res.nn_dist < 0.01).sum(), 100 * (res.nn_dist < 0.01).mean(), (res.nn_dist < 0.001).sum()))
P('nearest neighbour has the opposite orb/non-orb label (mixed pair): %d of %d (%.0f%%)' % (res.mixed_nn.sum(), len(res), 100 * res.mixed_nn.mean()))
close = res[res.nn_dist < 0.01]
P('  of the close pairs (nn dist <0.01): %d, mixed-label %d' % (len(close), close.mixed_nn.sum()))
pr = res[res.partner_retained == True]
P('partner copy retained in same alignment: %d tips; partner is nearest neighbour: %d; patristic dist to partner median %.4f; partner label differs from flagged tip: %d' % (len(pr), pr.partner_is_nn.sum(), pr.partner_dist.median(), (pr.my_label != pr.partner_label).sum()))
P('  partner dist <0.01: %d (mixed label: %d)' % ((pr.partner_dist < 0.01).sum(), ((pr.partner_dist < 0.01) & (pr.my_label != pr.partner_label)).sum()))
mx = res[res.mixed_nn & (res.nn_dist < 0.01)]
P('HOGs with >=1 close mixed-label pair (strong+moderate):', mx.HOG.nunique())
open(os.path.join(OUT, 'pairs_and_retention_summary.txt'), 'w').write(buf.getvalue())
