# HyPhy (stage 03) impact of suspected cross-sample contamination, 2026-11-09

Scripts: `scripts/03_selection_tests/contamination_impact/01..05_*.py` (run in order; step 1 needs the cached JSONs in
`results/hyphy_results_cache/`, step 1 output `tips_long.tsv` is not committed, 1.5M rows). Nothing was re-run on HyPhy or Rockfish.
Flag table: `data/flag_table_genes.csv` (strong / moderate / weak; flags are suspicions, not proof).
Raw text outputs: `comparison_summary.txt`, `pairs_and_retention_summary.txt`, `baseline_and_mixed_summary.txt`.

## Where the data are (evidence)
- Retained sequences and fg/bg labels for all 4,756 HOGs were read from the cached HyPhy JSONs (`input`, `branch attributes`, `tested`).
  The three tests use the identical tip set per HOG; labels match `orbweavers-list.txt` exactly (orb-fg: orb = test; non-orb-fg: non-orb = test;
  RELAX: Reference = orb, Test = non-orb).
- Missing locally: the alignments (`*_NT.trim.dedup.nex`) and the pre-filter alignments/removal logs. They are on Rockfish
  (`/scratch4/agordus1/crunnel2/hyphy_wd_260929/<HOG>/`); `inspect/` holds one HOG only. A pruned re-run needs them.

## Exposure
- 2,363 of the 16,063 flagged genes sit in the 4,756 HOGs; 2,048 survived into the HyPhy alignments.
- HOGs with >=1 retained flagged sequence: strong 331 (7.0%), strong+moderate 359 (7.5%), any incl. weak 660 (13.9%).
- Strong-flagged retained tips: 548 (Araneus inustus 312 = orb-weaver, so foreground in orb tests / RELAX Reference; the eight non-orb carriers 236 =
  foreground in non-orb BUSTED-PH / RELAX Test).
- Weak flags are mostly 78-85% nt identity ordinary orthologs (example: N5.HOG0068277); I treat strong+moderate as the contamination-relevant set.

## Mechanism (evidence)
- Survival: unflagged genes of the nine species are retained 87.6%, strong 80.6%, moderate 69.7%; flagged genes with 100% nt identity only 61.7%, and both
  copies of a 100%-identical pair survive in 6 of 214 cases. Exact duplicates are mostly collapsed by remove-duplicates (supports the instinct), but which of
  the pair survives is arbitrary, so the retained label can be either side. Not checked: the removal logs.
- Near-identical pairs that remain: of 631 retained strong+moderate tips, 64% have a nearest neighbour <0.01 subst/codon (baseline unflagged
  tips of the same species 18%) and 40% have that close neighbour on the OPPOSITE orb/non-orb label (baseline 0.0%). 132 HOGs have >=1 such
  close mixed-label pair. So the concern is real in structure: short branches whose sister is on the other label.

## Results (evidence)
- Group comparison (cached results, raw p<=0.05): affected HOGs are NOT enriched for hits. RELAX 40.5% (strong) vs 46.6% unaffected (adj. OR 0.77 for log tip count);
  BUSTED-PH comparative p<=0.05 26% vs 36% (orb), 28% vs 34% (non-orb). BUSTED-PH hit rule: 0/331 vs 2 hits; non-orb 1/331 vs 20.
- The 132 close-mixed-pair HOGs: RELAX 48.5% sig vs 46.6% unaffected (p=0.72); BUSTED-PH orb comparative 30% vs 36% (p=0.17). No detectable shift,
  but with n=132 an effect under about 10 percentage points is not detectable.
- Hit lists feeding GO: RELAX relaxed 609 HOGs contain 42 strong-affected (42.4 expected by chance), intensified 1,588 contain 92 (110.5 expected);
  BUSTED-PH orb 2 hits contain 0 strong/moderate-flagged sequences (HOG0068277 has 8 weak flags, ordinary orthologs); non-orb 21 hits contain 1.
  The enrichment hit lists are neither depleted nor enriched in a way that points to contamination.

## Guess / unproven
- The lower hit rate in "any" HOGs is probably confounded: weak-flag HOGs are conserved genes (median omega_test 0.04 vs 0.12), which have less power. Not tested.
- Group-level null results cannot show that individual HOGs are unchanged: a mixed close pair could push p either way. Only a pruned re-run answers that per HOG.
- Direction of transfer (which sample is the source) is unknown, so "flagged side" is not necessarily the contaminated side.
