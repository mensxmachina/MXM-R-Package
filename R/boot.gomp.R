boot.gomp <- function(target, dataset, tol = qchisq(0.95, 1), test = "testIndLogistic", method = "ar2", B = 500, ncores = 1) { 
  
  runtime <- proc.time()
  
  sel <- NULL
  n <- dim(dataset)[1]
  if ( !is.matrix(target) )  dim(target) <- c(n, 1)
  
  if ( ncores <= 1 ) {
    
    for (i in 1:B) {
      ina <- sample(n, n, replace = TRUE)
      sel <- c(sel, gomp(target[ina, ], dataset[ina, ], tol = tol, test = test, method = method)$res[-1, 1])
    }
    
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:B, .packages = "MXM", .export = "gomp") %dopar% {
      ina <- sample(n, n, replace = TRUE)
      sel <- MXM::gomp(target[ina, ], dataset[ina, ], tol = tol, test = test, method = method)$res[-1, 1]
      return( sel )
    }  
    stopCluster(cl)
    sel <- unlist(mod)
  }  
    
  res <- cbind(unique(sel), Rfast::Table(sel)/B )
  colnames(res) <- c("Variable", "Selection proportion")
  runtime <- proc.time() - runtime
  list(runtime = runtime, res = res)
}  
  
