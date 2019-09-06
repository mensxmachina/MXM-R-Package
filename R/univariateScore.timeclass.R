univariateScore.timeclass = function(target, dataset, test, wei = NULL, ncores = 1) {

    univariateModels <- list();
    dm <- dim(dataset)
    rows <- dm[1]
    cols <- dm[2]

    poia <- Rfast::check_data(dataset)
    if ( sum(poia) > 0 )   dataset[, poia] <- rnorm(rows * length(poia) )    
    nTests = cols
    univariateModels = NULL;
    univariateModels$pvalue = numeric(nTests) 
    univariateModels$stat = numeric(nTests)
    if ( ncores == 1 | is.null(ncores) | ncores <= 0 ) {
      
      for(i in 1:nTests) {
        test_results = test(target, dataset, i, 0, wei = wei)
        univariateModels$pvalue[[i]] = test_results$pvalue;
        univariateModels$stat[[i]] = test_results$stat;
      } 
    } else {
      #require(doParallel, quiet = TRUE, warn.conflicts = FALSE)  
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      test = test
      mod <- foreach(i = 1:nTests, .combine = rbind, .export = c("multinom"), .packages = "nnet") %dopar% {
        test_results = test(target, dataset, i, 0, wei = wei)
        return( c(test_results$pvalue, test_results$stat) )
      }
      stopCluster(cl)
      univariateModels$pvalue = as.vector( mod[, 1] )
      univariateModels$stat = as.vector( mod[, 2] )
    }
    
    if ( sum(poia>0) > 0 ) {
      univariateModels$stat[poia] = 0
      univariateModels$pvalue[poia] = log(1)
    }

    univariateModels
  }
  