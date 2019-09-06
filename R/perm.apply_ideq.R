perm.apply_ideq = function(i, queues, target, dataset, cvar, z, test, wei, threshold, univariateModels, hash, stat_hash, pvalue_hash, R = R)
{
  w = z[,i];
  w = t(t(w));
  zPrime = c(setdiff(z , w) , cvar);
  cur_results = test(target = target, dataset = dataset, xIndex = w, csIndex = zPrime, wei = wei, univariateModels = univariateModels, hash = hash, stat_hash = stat_hash, pvalue_hash = pvalue_hash, threshold = threshold, R = R);
  if ( cur_results$pvalue > threshold ) {
    queues[[w]] = as.matrix(c(queues[[w]] , queues[[cvar]]));
    return(queues[[w]]);
  } else  return(NA);
}