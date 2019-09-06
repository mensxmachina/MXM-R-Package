rdag <- function(n, p, s, a = 0, m = NULL, A = NULL, seed = FALSE) {
  ## n in the sample size
  ## p is the number of (nodes or) variables
  ## s is the success probability of the binomial distribution
  ## a is the percentage of outliers, is set to zero by default
  ## a number between 0 and 1
  ## m is the mean vector which is used only if you want outliers, i.e. if a > 0
  ## A is an adjancey matrix given by the user 
  if ( is.null(A) ) { ## no adjacency matrix is given
    if ( s > 1 || s < 0 )  s <- 0.5
    if ( a > 1 || a < 0 )  a <- 0
    if ( seed )  set.seed(1234567)
    A <- matrix( 0, p, p )
    nu <- 0.5 * p * (p - 1)
    A[ lower.tri(A) ] <- rbinom(nu, 1, s)
    A[ A == 1 ] <- runif( sum(A), 0.1, 1 )
  } else {
    A <- A
    p <- ncol(A)
  }
     
  Ip <- diag(p)
  sigma <- solve( Ip - A )
  sigma <- tcrossprod( sigma ) 
  nout <- 0
  if ( seed )  set.seed(1234567)
  
  if (a > 0) {
    y <- Rfast::rmvnorm( n - nout, numeric(p), sigma)
    nout <- round( a * n )
    yout <- Rfast::rmvnorm(nout, m, sigma)
    x <- rbind(y, yout)  
  } else  x <- Rfast::rmvnorm(n, numeric(p), sigma)
  
  G <- t( A )
  G[ G > 0 ] <- 2
  ind <- which( t(G) == 2 )
  G[ind] <- 3
  
  V <- colnames(A)
  if ( is.null(V) )   V <- paste("X", 1:p, sep = "")
  colnames(x) <- V
  colnames(G) <- rownames(G) <- V
  colnames(A) <- rownames(A) <- V
  list(nout = nout, G = G, A = A, x = x)
}