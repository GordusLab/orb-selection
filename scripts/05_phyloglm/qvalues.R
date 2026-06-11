library(qvalue)
results <- read.csv(here("results", "odds_ratio_test/Results_Apr16/Run1_occ_30-88_10000x_all_orb/perm_loss_pvalues.csv"))
pvals_fg <- results$pval_fg
qobj_fg <- qvalue(p = pvals_fg)
results$qval_fg <- qobj_fg$qvalues
pvals_bg <- results$pval_bg
qobj_bg <- qvalue(p = pvals_bg)
results$qval_bg <- qobj_bg$qvalues
write.csv(results, here("results", "odds_ratio_test/Results_Apr16/Run1_occ_30-88_10000x_all_orb/perm_loss_qvals.csv"), row.names = FALSE)
