rdag2 <- function(n, A = NULL, p, nei, low = 0.1, up = 1) {
  
  if ( is.null(A) ) {
    prob <- nei/(p - 1)
    A <- matrix(0, p, p)
    qa <- rbinom( 0.5 * p * (p - 1), 1,  prob )
    A[upper.tri(A)] <- qa
  }
  V <- colnames(A)
  if ( is.null(V) )   V <- paste("X", 1:p, sep = "")
  
  x <- matrix(0, n, p)
  x[, 1] <- rnorm(n)
   for (i in 2:p) {
    if ( sum( A[, i] != 0 ) == 0 ) {
      x[, i] <- rnorm(n)
    } else {
      id <- which(A[, i] == 1) 
      wa <- x[, id, drop = FALSE]
      ub <- runif( dim(wa)[2] )
      b <- runif( dim(wa)[2], -up, -low) * (ub<0.5) + runif( dim(wa)[2], low, up) * (ub>0.5)
      x[, i] <- rnorm(n, wa %*% b, 1)
      x[, i] <- ( x[, i] - mean(x[, i]) ) / Rfast::Var(x[, i], std = TRUE)
    }
  }
  colnames(A) <- rownames(A) <- V
  colnames(x) <- V
  list(G = A, x = x)
}
