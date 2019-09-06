fbed.gee.reg <- function(target, dataset, id, reps = NULL, ini = NULL, threshold = 0.05, wei = NULL, K = 0, test = "testIndGEEReg", correl = "exchangeable", se = "jack") { 
  

  if ( length(K) > 1 ) {
    
    res <- kfbed.gee.reg(y = target, x = dataset, id = id, reps = reps, univ = ini, alpha = threshold, wei = wei, K = K, test = test, correl = correl, se = se) 

  } else {
  
  dataset <- dataset[order(id), ]
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any(is.na(dataset) ) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
  }

  if ( is.null(reps) ) {
    if (test == "testIndGEEReg") {
      res <- fbed.geelm(y = target, x = dataset, id = id, univ = ini, alpha = threshold, wei = wei, K = K, correl = correl, se = se) 
    } else  res <- fbed.geeglm(y = target, x = dataset, id = id, univ = ini, alpha = threshold, wei = wei, K = K, test = test, correl = correl, se = se) 

  } else if ( !is.null(reps) ) {
    if (test == "testIndGEEReg") {
      res <- fbed.geelm.reps(y = target, x = dataset, id = id, reps = reps, univ = ini, alpha = threshold, wei = wei, K = K, correl = correl, se = se) 
    } else  res <- fbed.geeglm.reps(y = target, x = dataset, id = id, reps = reps, univ = ini, alpha = threshold, wei = wei, K = K, test = test, correl = correl, se = se) 
  }
  
  }  ## end if ( length(K) > 1 )
  
  res
}
  
  