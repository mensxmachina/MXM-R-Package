SES.timeclass <- function(target, reps, id, dataset, max_k = 3, threshold = 0.05, ini = NULL, wei = NULL,  
                          hash = FALSE, hashObject = NULL, ncores = 1) {
  ##############################
  # initialization part of MMPC 
  #############################
  runtime <- proc.time()
  stat_hash <- NULL;
  pvalue_hash <- NULL;
  
  if ( hash )  {
    if (is.null(hashObject) )  {
      stat_hash <- Rfast::Hash();
      pvalue_hash <- Rfast::Hash();
    } else if ( class(hashObject) == "list" ) {
      stat_hash <- hashObject$stat_hash;
      pvalue_hash <- hashObject$pvalue_hash;
    } else   stop('hashObject must be a list of two hash objects (stat_hash, pvalue_hash)')
  }
  ################################
  # test checking and initialize #
  ################################
  len <- length( unique(target) ) 
  if ( len == 2 ) {
    target <- as.numeric( as.factor(target) ) - 1
    ci_test <- "testIndTimeLogistic"
    test <- testIndTimeLogistic
  } else  {
    ci_test <- "testIndTimeMultinom"
    test <- testIndTimeMultinom  
  }
  
  dataset <- group.mvbetas(dataset, id, reps)
  la <- length( unique(id) )
  tar <- numeric(la)
  for(i in 1:la)   tar[i] <- unique( target[id == i] )
  target <- tar
  tar <- NULL
  
  if ( is.null(ini) )   ini <- univariateScore.timeclass(target = target, dataset = dataset, test = test, wei = wei, ncores = ncores)
  ###################################
  # options checking and initialize #
  ###################################
  max_k <- floor(max_k);
  varsize <- ncol(dataset);
  #option checking
  if ( (typeof(max_k)!="double") || max_k < 1 )   stop('invalid max_k option');
  if ( max_k > varsize )   max_k = varsize;
  if ( (typeof(threshold) != "double") || threshold < 0 || threshold >= 1 )   stop('invalid threshold option');
  #######################################################################################
  oop <- options(warn = -1) 
  on.exit( options(oop) )
  results <- InternalSES.timeclass(target, dataset, max_k, log(threshold), test, ini, wei, hash, varsize, stat_hash, pvalue_hash);
  #backward phase
  #varsToIterate <- results$selectedVarsOrder
  #if ( backward  & length( varsToIterate ) > 0  ) {
  #  varsOrder <- results$selectedVarsOrder
  #  bc <- mmpcbackphase(target, dataset[, varsToIterate, drop = FALSE], test = test, wei = wei, max_k = max_k, threshold = threshold)
  #  met <- bc$met
  #  results$selectedVars <- varsOrder[met]
  #  results$selectedVarsOrder = varsOrder[met]
  #  results$pvalues[varsToIterate] = bc$pvalues
  #  results$n.tests <- results$n.tests + bc$counter
  #}
  
  runtime <- proc.time() - runtime
  SES.timeclass.output <- new("SESoutput", selectedVars = results$selectedVars, selectedVarsOrder = results$selectedVarsOrder, queues = results$queues, 
                   signatures = results$signatures, hashObject = results$hashObject, pvalues = results$pvalues, stats = results$stats, univ = results$univ,
                   max_k = results$max_k, threshold = results$threshold, n.tests = results$n.tests, runtime = runtime, test = ci_test);
  
  return(SES.timeclass.output);
}


