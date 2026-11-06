# ORT TODO

Open items for the odds-ratio permulation test (`odds_ratio_test.py`). Nothing here is done unless
marked. Origin: discussion of where the ORT sits relative to other methods (phyloGLM, OrthoGLMM,
Forward Genomics, RERconverge permulations) and how it should work on transcriptome- or genome-based data.

## 1. Transcriptome copy number reflects expression (not only gene copy number)

Currently handled at the species level only: stage 01 CD-HIT/TransDecoder/BLASTP, BUSCO weights
(loss: fraction complete; duplication: fraction single-copy), duplication = count > 1, cell-level
contamination flags. Not handled: gene-specific detectability (expression, tissue/stage sampling),
and near-identical isoforms/alleles within a species and HOG.

- [ ] Collapse near-identical same-species copies within a HOG (e.g. >= 98% identity) before
      counting; rerun the duplication test as a sensitivity check.
- [ ] Compare transcriptome-derived HOG counts with genome-annotation counts for species that have
      both (*P. tepidariorum*, *U. diversus*) to measure how much counts reflect expression.
- [ ] Sensitivity run restricted to high-BUSCO species; check the hits hold.
- [ ] If reads are available, optional minimum-expression requirement for a copy to count.
- [ ] Longer term: gene-level detectability model (e.g. probability a gene is detected given
      expression breadth and tissues sampled).

## 2. Transcriptome- or genome-based input

- [ ] Document "genome mode": BUSCO weights optional (`buscos_filename=None` gives weights of 1),
      flags optional, input is just a gene-count matrix + phenotype lists + permulations.
- [ ] Make species-count and occupancy bounds configurable (currently tuned to the 98-species N5 set).
- [ ] Test on a genome-based example dataset.

## 3. Calibration and power (what the ORT must show to be taken seriously)

- [ ] Simulation study: ORT vs phyloGLM calibration and power by number of origins, effect size,
      and loss rate.
- [ ] Positive-control benchmark with a known loss/convergence signal (e.g. marine mammal
      chemoreceptor or other published gene losses).
- [ ] Permulation-model sensitivity: transition-count restriction, ARD vs ER, single-origin and
      easy-loss scenarios.
- [ ] Pooled calibration in prevalence bins (as in OrthoGLMM) tried on the existing 2,000-labeling
      phyloGLM nulls and on the ORT.

## 4. Positioning relative to other methods

- [ ] Write a short comparison (ORT vs phyloGLM, OrthoGLMM, Forward Genomics, treeWAS/Scoary,
      Pagel/fitPagel, CAFE5, RERconverge permulations): what is new (permulation null applied to
      weighted loss/duplication odds ratios, BUSCO weighting, contamination-flag handling, cheap
      per-test cost) and what is not.
- [ ] Record ORT vs phyloGLM cost (~5 orders of magnitude cheaper per test) with the measured numbers.

## 5. Tails and significance thresholds

Added 2026-11-06 from the ORT statistics session (the full-run equivalence check and results are
in the lab notebook entry and the session archive `~/claude-sessions/2026-11-04_ort_permulations_mtc`).

- [ ] Bring back a tail option (e.g. `--tails prespecified|both|less|greater`; `prespecified` as the
      default). The April 2026 script had `alternative="less"/"greater"` (one test per run, runs named
      `Loss_LT` / `Dup_RT`); the current code always computes all four tails.
      `prespecified` = loss in non-orb-weavers and duplication in orb-weavers only, each at the full
      alpha; `both` = the current behavior (alpha / 2 per direction).
- [ ] Decide which alpha the hit lists use. Today `fltrd_hits.csv`, the per-direction counts, the
      UpSet inputs and the GO lists test all four directions at alpha / 2 = 0.025, including the two
      directions that were not pre-specified (loss more likely in fg, duplication more likely in bg).
      The directional q-values (BH within each pre-specified tail, q < 0.05) are unaffected by alpha.
      A tail option would change which HOGs count as hits at the alpha level, not the directional q.
- [ ] Keep the stacked-tails q available under `both` (N5.HOG0071799 has stacked q = 0.066 but
      directional q = 0.033 at 500k and 1M, so the choice matters for that gene).

## 6. Provenance and housekeeping

- [ ] Record the exact batch range of the 500k permulation set. The 1M set is array batches 46-495
      (2,000 each) on top of the 100k set; the 500k batch range was not in the logged submit
      commands (likely the first ~245 extension batches). Check by comparing the first 500,000 rows of
      `data/perms_tip_values_500k.csv` with the 1M file, or the combine command on Rockfish
      (`combine_permulation_batches.py` has `--max-batch`).
- [ ] Widen the x-range of the 3a-3c and 5a plots together in `Odds Ratio Test Plots.ipynb` so the one
      HOG with log odds ratio below -6.6 is no longer clipped (both sets must be widened to stay
      comparable).
- [ ] First real use of the edited SLURM scripts (`run_odds_ratio_test_slurm.sh`,
      `permulations_batch_slurm.sh`): the `reports/<YYMMDD>_<job-name>/` folder must exist before
      `sbatch`, and the `PATH_TO_*` placeholders must be overridden on the command line (see
      `CLAUDE.md`).
