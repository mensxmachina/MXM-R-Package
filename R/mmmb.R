mmmb = function(target, dataset, max_k = 3, threshold = 0.05, test = "testIndFisher", user_test = NULL, ncores = 1) {
  
  durat <- proc.time()  
  
  mmpcobject <- MMPC(target, dataset, max_k = max_k, threshold = threshold, test = test, user_test = user_test, ncores = ncores, backward = TRUE)
  varsToIterate <- mmpcobject@selectedVars;
  pct <- varsToIterate
  met <- 1:length(pct)
  d <- dim(dataset)[2] 
  ci_test <- test <- mmpcobject@test
  lista <- list()
  aa <- NULL
  
  if ( length(pct) > 0 ) {
    
    for ( i in met) {
      tar <- dataset[, varsToIterate[i] ];
      datas <- cbind( dataset[, -varsToIterate[i] ], target)
      res <- MMPC(tar, datas, max_k = max_k, threshold = threshold, test = test, user_test = user_test, ncores = ncores, backward = TRUE) 
      poies <- sort( res@selectedVars )
      poies <- poies[poies != d ]
      poies[ poies >= varsToIterate[i] ] = poies[ poies >= varsToIterate[i] ] + 1
      lista[[ i ]] <- poies     
    }

  }
  runtime <- proc.time() - durat   
  
  list( mb = sort( c(mmpcobject@selectedVars[met], aa) ), ci_test = ci_test, runtime = runtime )
}