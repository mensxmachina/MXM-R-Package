identifyTheEquivalent.gee = function(queues, target, reps, group, dataset, cvar, z, test, wei, threshold, univariateModels, 
                                     pvalues, hash, stat_hash, pvalue_hash, correl, se) {
  z = t(z);
  #if we have more than one equivalent vars in z , we select the first one
  #for loop if we dont use the below lapply function
  for(i in 1:ncol(z)) {
    w = z[, i];
    w = t(t(w));
    zPrime = c(setdiff(z, w), cvar);
    cur_results = test(target, reps, group, dataset, w, zPrime, wei = wei, univariateModels, hash = hash, stat_hash, pvalue_hash, correl= correl, se = se);
    
    if ( cur_results$pvalue > threshold ) {  
      queues[[w]] = as.matrix( c(queues[[w]], queues[[cvar]]) );
      break;
      #equalsY = equalsY+1;
    }
  }
  #cat("\nEquals Ys in the current z subset = %d",equalsY);
  return(queues);
}