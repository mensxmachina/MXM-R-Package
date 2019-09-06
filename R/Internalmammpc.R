Internalmammpc = function(targ, data, statistic, max_k, threshold, test = NULL, ini = NULL, user_test = NULL, hash=FALSE, 
                          varsize, stat_hash, pvalue_hash, targetID, ncores)  {
  #get the current time
  runtime = proc.time();
  
  #######################################################################################
  
  D = length(targ) 
  
  if ( is.null(ini) ) { 
    
    cols = ncol( data[[ 1 ]] )  
    sta = pval = matrix(0, D, cols)   
    univariateModels = list()
    
    if ( identical(test, testIndFisher) )  { ## Pearson's correlation 
      
      for ( l in 1:D ) {  
        target = targ[[ l ]]
        dataset = data[[ l ]]
        rows = length(target)
        a = as.vector( cor(target, dataset) )
        dof = rows - 3; #degrees of freedom
        wa = abs( 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof) )
        if ( targetID != - 1 )   wa[ targetID ] = 0  
        sta[l, ] = wa
        pval[l, ] = log(2) + pt(-wa, dof, log.p = TRUE)
      }  
      univariateModels$stat = - 2 * Rfast::colsums(pval);
      univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
      univariateModels$stat_hash = stat_hash;
      univariateModels$pvalue_hash = pvalue_hash;
      
    } else if ( identical(test, testIndSpearman) ) {  ## Spearman's correlation
      
      for ( l in 1:D ) {  
        
        target = targ[[ l ]]
        dataset = data[[ l ]]
        rows = length(target) 
        a = as.vector( cor(target, dataset) )
        dof = rows - 3; #degrees of freedom
        wa = abs( 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof) / 1.029563 )
        if ( targetID != - 1 )  wa[ targetID ] = 0
        sta[l, ] = wa
        pval[l, ] = log(2) + pt(-wa, dof, lower.tail = FALSE, log.p = TRUE)
      }  
      univariateModels$stat = - 2 * Rfast::colsums(pval);
      univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
      univariateModels$stat_hash = stat_hash;
      univariateModels$pvalue_hash = pvalue_hash;    
      
    } else if ( identical(test, testIndMMReg) ) {  ## MM (Robust) linear regression
      
      univariateModels = list();
      fit1 = MASS::rlm(target ~ 1, maxit = 2000, method = "MM")
      lik2 = numeric(cols)
      dof = numeric(cols)
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
        for (l in 1:D) {
          
          target = targ[[ l ]]
          dataset = data[[ l ]]
          
          for ( i in 1:cols ) {
            
            if ( i != targetID ) {
              fit2 = MASS::rlm(target ~ dataset[, i], method = "MM", maxit = 2000 )
              lik2[i] = as.numeric( logLik(fit2) )
              dof[i] = length( coef(fit2) ) - 1
            } else {
              lik2[i] = lik1
            }
            
          } 
          lik1 = as.numeric( logLik(fit1) )
          sta[l, ] = as.vector( 2 * abs(lik1 - lik2) )
          pval[l, ] = pchisq(sta[l, ], dof, lower.tail = FALSE, log.p = TRUE)
        } 
        
        univariateModels$stat = - 2 * Rfast::colsums(pval)
        univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        
        for ( l in 1:D )	{
          
          target = targ[[ l ]]
          dataset = data[[ l ]]
          
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mod <- foreach(i = 1:cols, .combine = rbind, .packages = "MASS") %dopar% {
            
            if ( i != targetID ) {
              fit2 = rlm(target ~ dataset[, i] )
              lik2 = as.numeric( logLik(fit2) )
              
              return( c(lik2, length( coef(fit2) ) ) )
            } else{
              return( c(0, 0) )
            }
            
          }
          stopCluster(cl)
          
          lik1 = as.numeric( logLik(fit1) )
          dof = as.vector( mod[, 2] ) - 1 
          sta[l, ] = as.vector( 2 * abs(lik1 - mod[, 1]) )
          pval[l, ] = pchisq(sta[l, ], dof, lower.tail = FALSE, log.p = TRUE)
          
        }
        
        univariateModels$stat = - 2 * Rfast::colsums(pval)
        univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }   
      
      
    } else if ( identical(test, testIndLogistic)  &  is.ordered(target) ) {  ## 
      
      lik2 = numeric(cols)
      dof = numeric(cols)
      univariateModels = list();
      
      fit1 = ordinal::clm(target ~ 1)
      df1 = length( coef(fit1) )
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
        for ( l in 1:D )	{
          
          target = targ[[ l ]]
          dataset = data[[ l ]]
          
          for ( i in 1:cols ) {
            
            if (i != targetID){
              
              x <- model.matrix(target ~ dataset[, i] )
              fit2 <- ordinal::clm.fit(target, x)
              lik2[i] <- as.numeric( fit2$logLik )
              dof[i] <- length( coef(fit2) ) - df1
            } else {
              lik2[i] <- lik1
            }   
          }
          
          lik1 = as.numeric( logLik(fit1) )
          sta[l, ] = as.vector( 2 * abs(lik1 - lik2) )
          pval[l, ] = pchisq(sta[l, ], dof, lower.tail = FALSE, log.p = TRUE)
          
        }
        
        univariateModels$stat = - 2 * Rfast::colsums(pval)
        univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        
        for ( l in 1:D )	{
          
          target = targ[[ l ]]
          dataset = data[[ l ]]
          
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mod <- foreach(i = 1:cols, .combine = rbind, .packages = "ordinal") %dopar% {
            ## arguments order for any CI test are fixed
            if ( i != targetID ) {
              x <- model.matrix(target ~ dataset[, i] )
              fit2 <- ordinal::clm.fit(target, x)
              lik2 <- as.numeric( fit2$logLik )
              
              return( c(lik2, length( coef(fit2) ) ) )
            } else{
              return( c(0, 0) )
            }
          }
          stopCluster(cl)
          
          lik1 = as.numeric( logLik(fit1) )
          dof = as.vector( mod[, 2] ) - df1 
          sta[l, ] = as.vector( 2 * abs(lik1 - mod[, 1]) )
          pval[l, ] = pchisq(sta[l, ], dof, lower.tail = FALSE, log.p = TRUE)
          
        }
        
        univariateModels$stat = - 2 * Rfast::colsums(pval)
        univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }
      
    } else {  
      univariateModels = univariateScore.ma(targ, data, test, statistic = FALSE, hash=hash, stat_hash=stat_hash, pvalue_hash=pvalue_hash, targetID=targetID, ncores=ncores);
    }
    
  } else {
    univariateModels = ini
  } 
  
  pvalues = univariateModels$pvalue;      
  stats = univariateModels$stat;
  #if we dont have any associations , return
  if ( min(pvalues , na.rm = TRUE) > threshold ) {
    results = NULL;
    results$selectedVars = c();
    class(results$selectedVars) = "numeric";
    results$selectedVarsOrder = c();
    class(results$selectedVarsOrder) = "numeric";
    results$hashObject = NULL;
    class(results$hashObject) = 'list';
    class(results$univ) = 'list';
    results$pvalues = exp(pvalues);
    results$stats = stats;
    results$univ = univariateModels
    results$max_k = max_k;
    results$threshold = threshold;
    runtime = proc.time() - runtime;
    results$runtime = runtime;
    return(results);
  }
  #Initialize the data structs
  selectedVars = numeric(varsize);
  selectedVarsOrder = numeric(varsize)
  #select the variable with the highest association
  selectedVar = which( pvalues == pvalues[[which.min(pvalues)]] );
  selectedVars[selectedVar] = 1;
  selectedVarsOrder[selectedVar] = 1; #CHANGE
  #remaining variables to be considered
  remainingVars = numeric(varsize) + 1;
  remainingVars[selectedVar] = 0;
  remainingVars[pvalues > threshold] = 0;
  if (targetID > 0)  remainingVars[targetID] = 0;
  ################ main MMPC loop ################
  #loop until there are not remaining vars
  loop = any(as.logical(remainingVars));
  #rep = 1;
  while(loop) {
    #lets find the variable with the max min association
    max_min_results = max_min_assoc.ma(target, dataset, test, threshold, statistic, max_k, selectedVars, pvalues, stats, remainingVars, univariateModels, selectedVarsOrder, hash=hash, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    selectedVar = max_min_results$selected_var;
    selectedPvalue = max_min_results$selected_pvalue;
    remainingVars = max_min_results$remainingVars;
    pvalues = max_min_results$pvalues;
    stats = max_min_results$stats;
    stat_hash=max_min_results$stat_hash;
    pvalue_hash=max_min_results$pvalue_hash;
    #if the selected variable is associated with target , add it to the selected variables
    if(selectedPvalue <= threshold) {
      #print(paste("rep: ",rep,", selected var: ",selectedVar,", pvalue = ",exp(selectedPvalue)))
      #rep = rep + 1;
      selectedVars[selectedVar] = 1;
      selectedVarsOrder[selectedVar] = max(selectedVarsOrder) + 1;
      remainingVars[selectedVar] = 0;
    }
    loop = any(as.logical(remainingVars));
  }
  
  
  selectedVarsOrder[which(!selectedVars)] = varsize;#
  numberofSelectedVars = sum(selectedVars);#
  selectedVarsOrder = sort(selectedVarsOrder);#
  #adjusting the results
  if(targetID > 0)  {
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
  results$pvalues = pvalues;
  results$stats = stats;
  results$univ = univariateModels
  results$max_k = max_k;
  results$threshold = threshold;
  runtime = proc.time() - runtime;
  results$runtime = runtime;
  
  return(results);
}


