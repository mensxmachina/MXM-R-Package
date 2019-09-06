partialcor <- function(R, indx, indy, indz, n) {
  ## R is a correlation matrix
  ## i and j denote the two variables whose conditional correlation is to be estimated 
  ## k denotes the set of conditioning variables
  if ( is.null(indz) || indz == 0  )  {## classical correlation
     r <- R[indx, indy]
  } else if ( length(indz) == 1  &  all( indz > 0 ) ) {
    a1 <- R[indx, indy]     ;    a2 <- R[indx, indz]
    a3 <- R[indy, indz]
    r <- (a1 - a2 * a3) / sqrt( (1 - a3^2) * (1 - a2^2) )
    
  } else if ( length(indz) > 1 ) {
     rho <- solve( R[c(indx, indy, indz), c(indx, indy, indz)] )
     r <-  - rho[1, 2] / sqrt(rho[1, 1] * rho[2, 2])
  }
    if ( abs(r) > 1 ) r <- 0.99999
    z <- 0.5 * log( (1 + r) / (1 - r) ) * sqrt( n - sum(indz > 0) - 3 )
    pvalue <- pt( abs(z), n - sum(indz > 0) - 3, lower.tail = FALSE, log.p = TRUE )
    res <- c(r, pvalue)
    names(res) <- c("partial cor", "p-value")
    res
}

