univariateScore.glmm = function(target, reps = NULL, group, dataset, test, wei, targetID, slopes, ncores) {
  univariateModels <- list();
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]
  ind <- 1:cols
  
  la <- length( unique(target) )
  if ( la > 2  &  sum( round(target) - target ) != 0  &  !slopes  &  is.null(wei) ) {
    group <- as.numeric(group)
    if ( !is.null(reps) )   reps <- as.numeric(reps)
    univariateModels <- MXM::rint.regs(target = target, dataset = dataset, targetID = targetID, id = group, reps = reps, tol = 1e-08)

  } else {
  
    if (targetID != -1 ) {
      target <- dataset[, targetID]
      dataset[, targetID] <- rnorm(rows)
    }
  
    poia <- Rfast::check_data(dataset)
    if ( sum(poia) > 0 )   ind[poia] <- 0
    
    univariateModels$pvalue <- numeric(cols) 
    univariateModels$stat <- numeric(cols)
    
    if ( ncores == 1 | is.null(ncores) | ncores <= 0 ) {
    
      for(i in ind) {
        test_results <- test(target, reps, group, dataset, i, 0, wei = wei, slopes = slopes)
        univariateModels$pvalue[[i]] <- test_results$pvalue;
        univariateModels$stat[[i]] <- test_results$stat;
      } 
    } else {
      #require(doParallel, quiet = TRUE, warn.conflicts = FALSE)  
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      test <- test
      mod <- foreach(i = ind, .combine = rbind, .export = c("lmer", "glmer"), .packages = "lme4") %dopar% {
        test_results <- test(target, reps, group, dataset, i, 0, wei = wei, slopes = slopes)
        return( c(test_results$pvalue, test_results$stat) )
      }
      stopCluster(cl)
      univariateModels$pvalue[ind] <- as.vector( mod[, 1] )
      univariateModels$stat[ind] <- as.vector( mod[, 2] )
    }

  }
  if ( !is.null(univariateModels) )  {
    if (targetID != - 1) {
      univariateModels$stat[targetID] <- 0
      univariateModels$pvalue[targetID] <- log(1)
    }
  }
  
  univariateModels
}
