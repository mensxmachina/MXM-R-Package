rmdag <- function(n, A = NULL, p, nei, low = 0.1, up = 1) {
  
  if ( is.null(A) ) {
    prob <- nei/(p - 1)
    G <- matrix(0, p, p)
    qa <- rbinom( 0.5 * p * (p - 1), 1, prob )
    G[upper.tri(G)] <- qa
    V <- paste("X", 1:p, sep = "")
  } else {
    G <- A
    V <- colnames(G)
  }  
  
  x <- data.frame(matrix(vector(), n, p) )
  u <- runif(1)
  if ( u < 0.5 ) {
    x[, 1] <- rnorm(n)
  } else {
    k <- sample(1:4, 1)
    x[, 1] <- factor( rbinom(n, k, 0.5), ordered = TRUE)
  }
  for (i in 2:p) {
    if ( sum( G[, i] != 0 ) == 0 ) {
      u <- runif(1)
      if ( u < 0.5 ) {
        x[, i] <- rnorm(n)
      } else {
        k <- sample(1:4, 1)
        x[, i] <- factor( rbinom(n, k, 0.5), ordered = TRUE )
      }
    } else {
      id <- which(G[, i] == 1) 
      wa <- x[, id, drop = FALSE]
      wa <- model.matrix(~., data.frame(wa))[1:n, ]  
      ub <- runif( dim(wa)[2] )
      b <- runif( dim(wa)[2], -1, -0.1) * (ub<0.5) + runif( dim(wa)[2], 0.1, 1) * (ub>0.5)
      x[, i] <- rnorm(n, wa %*% b, 1)
      x[, i] <- ( x[, i] - mean(x[, i]) ) / Rfast::Var(x[, i], std = TRUE)
      u <- runif(1)
      if ( u > 0.5 ) { 
        k <- sample(1:4, 1)
        pou <- rep(0.15, k + 1)
        u <- runif(k + 1)  
        u <- u / sum(u) * ( 1 - sum(pou) )
        pou <- round( (pou + u) * n )
        pou <- cumsum(pou) / n
        pou <- pou/sum(pou)
        ep <- quantile(x[, i], pou)
        poia <- findInterval(x[, i], ep)
        x[, i] <- factor(poia, ordered = TRUE)
      }
    }
  }
  if ( is.null(V) )   V <- paste("X", 1:p, sep = "")
  colnames(x) <- V
  colnames(G) <- rownames(G) <- V
  list(G = G, x = x)
}