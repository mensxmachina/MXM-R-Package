max_min_assoc = function(target, dataset, test, wei, threshold, max_k, selectedVars, pvalues, stats, remainingVars, 
                         univariateModels, selectedVarsOrder, hash, stat_hash, pvalue_hash) {
  #Initialize
  selected_var = -1;
  selected_pvalue = 2;
  selected_stat = 0;
  
  varsToIterate = which(remainingVars==1);
  
  for (cvar in varsToIterate) {
    mma_res = min_assoc(target, dataset, test, max_k, cvar, wei, selectedVars, pvalues, stats, univariateModels, selectedVarsOrder, hash, stat_hash, pvalue_hash);
    pvalues = mma_res$pvalues;
    stats = mma_res$stats;
    stat_hash = mma_res$stat_hash;
    pvalue_hash = mma_res$pvalue_hash;
    if (mma_res$pvalue > threshold)  remainingVars[[cvar]] = 0;
    if ( compare_p_values(mma_res$pvalue, selected_pvalue, mma_res$stat, selected_stat) ) {
      selected_var = cvar;
      selected_pvalue = mma_res$pvalue;
      selected_stat = mma_res$stat;
    }
  }
  
  results <- list(selected_var = selected_var , selected_pvalue = selected_pvalue, remainingVars = remainingVars , pvalues = pvalues , stats = stats, stat_hash=stat_hash, pvalue_hash = pvalue_hash);
  return(results); 
}
