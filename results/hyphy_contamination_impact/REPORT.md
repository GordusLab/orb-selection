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

## Pruned re-run (359 HOGs, flagged strong/moderate tips removed; RELAX + BUSTED-PH orb fg), 2026-11-10
Scripts `06..10_*`; table `pruned_summary.tsv`, comparison `pruned_comparison_summary.txt`, `pruned_vs_original.csv`.
Status at the time of this comparison: 357 of 359 RELAX and 358 of 359 BUSTED-PH orb runs finished (3 tier-B/C tasks still running or resubmitted; they are not included). Only the first-fit tips were removed (median 1 per HOG, max 8); trees were pruned, not re-inferred.

**Evidence**
- RELAX (n=357): 144 significant before, 141 after. 25 HOGs changed status (14 lost, 11 gained, so no net direction); one HOG switched relaxed -> intensified. Of the 25, 23 are in HOGs with a strong flagged tip, 2 in moderate-only HOGs (strong HOGs are most of the set, 329 of 357, so this is not informative).
- Effect sizes barely move: Spearman of log2 k 0.915, median |change| 0.015; -log10 p Spearman 0.88.
- Most status changes are near the threshold or at k within about 0.1 of 1 (e.g. p 0.03 -> 0.15). Of 41 HOGs with an original p between 0.02 and 0.1, 8 changed status.
- A few large changes: HOG0058602, 0059807, 0068520 (not significant -> relaxed with k 0.19-0.36), HOG0066754 and 0060452 (strongly relaxed -> not significant), HOG0043408, 0048515, 0047291 (k shifts). HOG0066754 and 0060452 had 1 and 4 strong-flagged tips removed.
- BUSTED-PH orb: 0 hits before and 0 after in these HOGs (matches the earlier 0/331). Comparative p <= 0.05: 95 -> 91 (9 lost, 5 gained). Spearman of comparative p 0.865.

**Guess / unproven**
- I did not re-run unpruned HOGs under identical settings, so the run-to-run noise of HyPhy optimisation is unknown. Part or all of the near-threshold flips could be optimisation noise rather than a contamination effect.
- The large single-HOG changes (p of 1.0 or ~0 with extreme k, three with k 0.19-0.36 after pruning, original k 1.3583 with p 1.0 for HOG0068520) look like unstable fits (the runs use TOLERATE_NUMERICAL_ERRORS=1); they should be checked individually before being interpreted.
- Overall reading: no sign of a systematic effect on RELAX or BUSTED-PH conclusions; at most a handful of individual HOGs (the large-change list above) deserve a closer look.
- Not done: IQ-TREE re-inference on pruned alignments (prune before dedup); the two remaining RELAX tasks; the non-orb BUSTED-PH test.
