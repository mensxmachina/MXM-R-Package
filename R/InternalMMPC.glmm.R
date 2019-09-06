InternalMMPC.glmm = function(target, reps, group, dataset, max_k, threshold, test = NULL, ini, wei, user_test = NULL,  
                                 hash=FALSE, varsize, stat_hash, pvalue_hash, targetID, slopes, ncores) {
  #get the current time
  #######################################################################################
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]
  #univariate feature selection test
  if ( is.null(ini) ) {
    univariateModels <- glmm.univregs(target = target, reps = reps, id = group, dataset = dataset, targetID = targetID, test = test, wei = wei, 
                              slopes = slopes, ncores = ncores)
  } else  univariateModels <- ini
  
  pvalues <- univariateModels$pvalue;
  stats <- univariateModels$stat;
  #if we dont have any associations , return
  if ( min(pvalues, na.rm = TRUE) > threshold ) {   #or min(pvalues, na.rm=TRUE)
    #cat('No associations!');
    results = NULL;
    results$selectedVars = c();
    class(results$selectedVars) = "numeric";
    results$selectedVarsOrder = c();
    class(results$selectedVarsOrder) = "numeric";
    results$hashObject = NULL;
    class(results$hashObject) = 'list';
    class(results$univ) = 'list';
    results$pvalues = pvalues;
    results$stats = stats;
    results$univ = univariateModels
    results$max_k = max_k;
    results$threshold = threshold;
    results$slope = slopes
    results$n.tests = length(stats)
    
    return(results);
  }
  #Initialize the data structs
  selectedVars = numeric(varsize);
  selectedVarsOrder = numeric(varsize);
  #select the variable with the highest association
  selectedVar = which.max(stats)
  selectedVars[selectedVar] = 1;
  selectedVarsOrder[selectedVar] = 1; #CHANGE
  #lets check the first selected var
  #cat('First selected var: %d, p-value: %.6f\n', selectedVar, pvalues[selectedVar]);
  #remaining variables to be considered
  remainingVars = numeric(varsize) + 1;
  remainingVars[selectedVar] = 0;
  remainingVars[pvalues > threshold] = 0;
  if (targetID > 0)    remainingVars[targetID] = 0;
  # main MMPC.glmm loop
  # loop until there are not remaining vars
  loop = any(as.logical(remainingVars));
  
  while (loop) {
    #lets find the variable with the max min association
    max_min_results = max_min_assoc.glmm(target, reps, group, dataset, test, wei, threshold, max_k, selectedVars, pvalues, stats, remainingVars, 
                                         univariateModels, selectedVarsOrder, hash=hash, stat_hash=stat_hash, pvalue_hash=pvalue_hash, slopes = slopes);
    selectedVar = max_min_results$selected_var;
    selectedPvalue = max_min_results$selected_pvalue;
    remainingVars = max_min_results$remainingVars;
    pvalues = max_min_results$pvalues;
    stats = max_min_results$stats;
    stat_hash=max_min_results$stat_hash;
    pvalue_hash=max_min_results$pvalue_hash;
    #if the selected variable is associated with target , add it to the selected variables
    if (selectedPvalue <= threshold) {
      selectedVars[selectedVar] = 1;
      selectedVarsOrder[selectedVar] = max(selectedVarsOrder) + 1;
      remainingVars[selectedVar] = 0;
    }
    loop = any(as.logical(remainingVars));
  }
  
  selectedVarsOrder[which(!selectedVars)] = varsize;
  numberofSelectedVars = sum(selectedVars); 
  selectedVarsOrder = sort(selectedVarsOrder);
  #   selectedVars = selectedVarsOrder[1:numberofSelectedVars];
  #adjusting the results
  if (targetID > 0)  {
    toAdjust <- which(selectedVars > targetID);
    selectedVars[toAdjust] = selectedVars[toAdjust] + 1;
  }
  
  results = NULL;
  results$selectedVars = which(selectedVars == 1);
  svorder = sort(pvalues[results$selectedVars] , index.return = TRUE);
  svorder = results$selectedVars[svorder$ix];
  results$selectedVarsOrder = svorder;
  hashObject = NULL;
  hashObject$stat_hash = stat_hash;
  hashObject$pvalue_hash = pvalue_hash;
  results$hashObject = hashObject;
  class(results$hashObject) = 'list';
  class(results$univ) = 'list';
  results$pvalues = pvalues;
  results$stats = stats;
  results$univ = univariateModels
  results$max_k = max_k;
  results$threshold = threshold;
  results$slope = slopes
  results$n.tests = length(stats) + length( hashObject$stat_hash )
  
  return(results);
}
