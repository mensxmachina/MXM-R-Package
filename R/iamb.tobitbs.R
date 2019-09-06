iamb.tobitbs <- function(target, dataset, threshold = 0.05, wei = NULL) {
  
  threshold <- log(threshold)
  dm <- dim(dataset)
  n <- dm[1]  ## sample size 
  p <- dm[2]  ## number of variables
  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. No backward procedure was attempted")
  } else {

    a1 <- internaliamb.tobitbs( target = target, dataset = dataset, threshold = threshold, wei = wei, p = p ) 
    ind <- 1:p
    a2 <- list()
    poies <- a1$mat[, 1]
    if ( length(poies) > 0 ) {
      ind[-poies] <- 0
      ind <- ind[ind > 0]
      dat <- dataset[, poies, drop = FALSE ]
      a2 <- internaliamb.tobitbs( target = target, dataset = dat, threshold = threshold, wei = wei, p = length(ind) ) 
      poies <- a2$mat[, 1]
      ind[-poies] <- 0
      ind <- ind[ind > 0]
      dat <- dat[, poies, drop = FALSE]
      i <- 2
    } else {
      ind <- NULL
      a2$mat <- NULL  
    }
    while ( length(a1$mat[, 1]) - length(a2$mat[, 1]) != 0 ) {
      i <- i + 1
      a1 <- a2
      a2 <- internaliamb.tobitbs( target = target, dataset = dat, threshold = threshold, wei = wei, p = length(ind) ) 
      poies <- a2$mat[, 1]
      if ( length(poies) > 0 ) {
        ind[-poies] <- 0
        ind <- ind[ind > 0]
        dat <- dat[, poies, drop = FALSE]
      } else {
        ind <- NULL
        dat <- NULL  
      }  
    }
    
    res <- list(info = ind, mat = a2$mat, ci_test = "testIndTobit", final = a2$final ) 
  }
  
  res
}  









