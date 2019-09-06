univariateScore.gee <- function(target, reps = NULL, group, dataset, test, wei, targetID, correl, se, ncores) {
  univariateModels <- list();
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]
  ind <- 1:cols
  
  if (targetID != -1 ) {
    target <- dataset[, targetID]
    dataset[, targetID] <- rnorm(rows)
  }
    
  poia <- Rfast::check_data(dataset)
  if ( sum(poia) > 0 )   ind[poia] <- 0

  univariateModels$pvalue <- numeric(cols) 
  univariateModels$stat <- numeric(cols)
  test_results = NULL;
  if ( ncores == 1 | is.null(ncores) | ncores <= 0 ) {
    
    for(i in ind) {
      test_results <- test(target, reps, group, dataset, i, 0, wei = wei, correl = correl, se = se)
      univariateModels$pvalue[[i]] <- test_results$pvalue;
      univariateModels$stat[[i]] <- test_results$stat;
    } 
  } else {
    #require(doParallel, quiet = TRUE, warn.conflicts = FALSE)  
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    test = test
    mod <- foreach(i = ind, .combine = rbind, .export = c("geeglm", "ordgee"), .packages = "geepack") %dopar% {
      test_results = test(target, reps, group, dataset, i, 0, wei = wei, correl = correl, se = se)
      return( c(test_results$pvalue, test_results$stat) )
    }
    stopCluster(cl)
    univariateModels$pvalue[ind] <- as.vector( mod[, 1] )
    univariateModels$stat[ind] <- as.vector( mod[, 2] )
  }
  
  if (targetID != - 1) {
    univariateModels$stat[targetID] <- 0
    univariateModels$pvalue[targetID] <- log(1)
  }

  univariateModels
}
