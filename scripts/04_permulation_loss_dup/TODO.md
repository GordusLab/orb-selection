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
