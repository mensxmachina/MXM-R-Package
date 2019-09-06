apply_ideq.ma = function(i, queues, target, dataset, cvar, z, test, threshold, statistic, univariateModels, hash, stat_hash, pvalue_hash)
{
  w = z[,i];
  w = t(t(w));
  zPrime = c(setdiff(z , w) , cvar);
  cur_results = test(target = target, dataset = dataset, xIndex = w, csIndex = zPrime, statistic = statistic, univariateModels = univariateModels, hash = hash, stat_hash = stat_hash, pvalue_hash = pvalue_hash);
  
  if ( cur_results$pvalue > threshold ) {
    queues[[w]] = as.matrix(c(queues[[w]] , queues[[cvar]]));
    return(queues[[w]]);
  } else {
    return(NA);
  }
}