group.mvbetas <- function(x, id, reps) {
 
  dm <- dim(x)
  p <- dm[2]
  len <- length( unique(id) )
  B <- matrix(0, len, ncol = 2 * p )
  
  if ( p > 1 ) {
    for ( i in 1:len ) {
      mod <- Rfast::mvbetas( x[id == i, ], reps[id==i] )
      B[i, 1:p] <- mod[, 1]
      B[i, c(p + 1):c(2 * p)] <- mod[, 2]
    }
  } else {
    for ( i in 1:len ) {
      y <- x[id ==i, , drop = FALSE]
      z <- reps[id == i]
      r <- as.vector( cov(z, y) )
      n <- length(z)
      my <- Rfast::colmeans(y)
      mz <- sum(z)/n
      sz <- (sum(z^2) - sum(z)^2/n)/(n - 1)
      be <- r/sz
      a <- my - be * mz
      B[i, 1] <- a
      B[i, 2] <- be
    }
  }
  B <- rbind(B[, 1:p, drop = FALSE], B[, c(p + 1):c(2 * p), drop = FALSE])
}  