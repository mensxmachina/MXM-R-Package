ord.resid <- function(y, est) {
  y <- as.numeric(y)
  res <- y
  dm <- dim(est)
  n <- dm[1]
  d <- dm[2]
  ind1 <- which(y == 1)
  indd <- which(y == d)  
  ind <- setdiff( 1:n, c(ind1, indd) )
  res[ind1] <- est[ind1, 1] - 1
  res[indd] <- 1 - est[indd, d]
  w <- 1 - Rfast::design_matrix(y, ones = FALSE)
  w <- (w * est)[ind, ]
  a <- numeric( dim(w)[1] ) 
  for ( i in 1:dim(w)[1] ) {
    a[i] <- sum(w[i, 1:(y[ind[i]] - 1)]) - sum(w[i, (y[ind[i]] + 1):d])
  }
  res[ind] <- a
  res
}
