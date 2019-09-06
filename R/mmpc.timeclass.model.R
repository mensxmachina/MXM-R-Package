mmpc.timeclass.model <- function(target, dataset, id, reps, wei = NULL, mmpctimeclass.Object) {
  
  signature <- mmpctimeclass.Object@selectedVars
  if ( sum( is.na(signature) ) > 0 | length(signature) == 0 )  {
    mod = paste("No associations were found, hence no model is produced.")
    signature = NULL

  } else {
    
    if ( any(is.na(dataset) ) ) {
      warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    }
    
    dat <- dataset[, signature, drop = FALSE]
    x <- group.mvbetas(dat, id, reps)
    n <- 0.5 * dim(x)[1]
    xx <- cbind(x[1:n, ], x[(n + 1):(2 * n), ]) 
    la <- length( unique(id) )
    tar <- numeric(la)
    for(i in 1:la)   tar[i] <- unique( target[id == i] )
    target <- tar
    tar <- NULL
  
    if ( is.null( colnames(dataset) ) ) {  
      names(signature) = paste("Var", signature, sep = " ")
    } else  names(signature) = colnames(dataset)[signature]
  
    if ( mmpctimeclass.Object@test == "testIndTimeLogistic" ) {
      mod <- glm(target ~ xx, binomial, weights = wei)
    } else  mod <- nnet::multinom(target ~ xx, weights = wei, trace = FALSE)
    
    signature <- c( signature, BIC(mod) )
    names(signature)[length(signature)] = "bic" 
  }  ## end if ( sum( is.na(signature) ) > 0 )
  list(mod = mod, signature = signature)  
  
}