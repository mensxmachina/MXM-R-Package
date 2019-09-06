IdentifyEquivalence.ma = function(queues, target, dataset, test, threshold, statistic, max_k, selectedVars, pvalues, stats, remainingVars, univariateModels, selectedVarsOrder, hash, stat_hash, pvalue_hash)
{ 
  varsToBeConsidered = which(selectedVars==1 | remainingVars==1); #CHANGE
  lastvar = which(selectedVarsOrder == max(selectedVarsOrder))[1]; #CHANGE
  
  #for every variable to be considered
  for (cvar in varsToBeConsidered) {
    #if var is the last one added, no sense to perform the check
    if (cvar == lastvar)  next;
    #the maximum conditioning set
    selectedVarsIDs = which(selectedVars == 1);
    cs = setdiff(selectedVarsIDs , cvar);
    k = min( c(max_k , length(cs)) );
    #for all possible subsets at most k size
    #problem in 1:k if K=0 - solve with while temporary
    klimit = 1;
    while (klimit <= k) {
      #set to investigate
      tempCS = setdiff(cs, lastvar)#CHANGE
      if (klimit == 1) {  #CHANGE
        subsetcsk = as.matrix(lastvar); #CHANGE
      } else {
        subsetcsk = as.matrix( nchoosek(tempCS, klimit - 1) )   
        numSubsets = dim(subsetcsk)[2]; 
        subsetcsk = rbind(subsetcsk, lastvar*rep(1,numSubsets));#CHANGE
      }
      #flag to get out from the outermost loop
      breakFlag = FALSE;
      
      for (i in 1:ncol(subsetcsk) ) {
        z = subsetcsk[,i];
        z = t(t(z));
        
        cur_results = test(target = target, dataset = dataset, xIndex = cvar, csIndex = z, statistic = statistic, univariateModels = univariateModels, hash = hash, stat_hash = stat_hash, pvalue_hash = pvalue_hash);
        stat_hash = cur_results$stat_hash;
        pvalue_hash = cur_results$pvalue_hash;
        #check if the pvalues and stats should be updated
        if ( compare_p_values(pvalues[[cvar]], cur_results$pvalue, stats[[cvar]], cur_results$stat) )  {
          pvalues[[cvar]] = cur_results$pvalue;
          stats[[cvar]] = cur_results$stat;
        }
        #if there is any subset that make var independent from y,
        #then let's throw away var; moreover, we also look for
        #equivalent variables. Note that we stop after the first 
        #z such that pvalue_{var, y | z} > threshold
        if ( cur_results$pvalue > threshold ) {
          remainingVars[[cvar]] = 0;
          selectedVars[[cvar]] = 0;
          queues = identifyTheEquivalent.ma(queues, target, dataset, cvar, z, test, threshold, statistic, univariateModels, pvalues, hash, stat_hash, pvalue_hash);
          breakFlag = TRUE;
          break;
        }
      }
      
      if( breakFlag ) { 
        break;
      } else {
        klimit = klimit + 1;
      }
    }
  }
  results <- list(pvalues = pvalues, stats = stats, queues = queues, selectedVars = selectedVars, remainingVars = remainingVars, stat_hash = stat_hash, pvalue_hash = pvalue_hash);
  return(results);
}