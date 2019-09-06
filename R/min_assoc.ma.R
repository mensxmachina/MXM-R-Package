min_assoc.ma = function(target, dataset, test, max_k, cvar, statistic, selectedVars, pvalues, stats, univariateModels , selectedVarsOrder, hash, stat_hash, pvalue_hash)
{

  ma_pvalue = pvalues[[cvar]]; #CHANGE
  ma_stat = stats[[cvar]]; #CHANGE
  selectedVars = which(selectedVars==1);
  #max size of the condiotioning test
  k = min(c(max_k , length(selectedVars)));
  
  ck = 1;
  while ( ck <= k ) {
    #lastvar = unique(which(selectedVarsOrder == max(selectedVarsOrder)));
    lastvar = which(selectedVarsOrder == max(selectedVarsOrder))[1]; #CHANGE
    tempCS = setdiff(selectedVars, lastvar) #CHANGE
    if ( ck == 1 ) {   #CHANGE
      subsetcsk = as.matrix(lastvar); #CHANGE
    } else {
      subsetcsk = as.matrix( nchoosek(tempCS, ck - 1) )
      numSubsets = dim(subsetcsk)[2]; #CHANGE
      subsetcsk = rbind(subsetcsk, lastvar*rep(1,numSubsets)); #CHANGE
    }
    
    
    for (i in 1:ncol(subsetcsk) ) {
      s = subsetcsk[,i];
      cur_results = test(target = target, dataset = dataset, xIndex = cvar, csIndex = s, statistic = statistic, univariateModels = univariateModels, hash = hash, stat_hash = stat_hash, pvalue_hash = pvalue_hash);
      stat_hash = cur_results$stat_hash;
      pvalue_hash = cur_results$pvalue_hash;
      #check if the pvalues and stats should be updated
      if ( !compare_p_values(cur_results$pvalue, ma_pvalue, cur_results$stat, ma_stat) ) {
        ma_pvalue = cur_results$pvalue;
        pvalues[[cvar]] = cur_results$pvalue;
        ma_stat = cur_results$stat;
        stats[[cvar]] = cur_results$stat;
      }
    }
    ck = ck + 1;
  }
  results <- list(pvalue = ma_pvalue , stat = ma_stat , pvalues = pvalues , stats = stats, stat_hash=stat_hash, pvalue_hash  = pvalue_hash);
  return(results);
}