big.score.univregs <- function(target = NULL, dataset, test) {
  
  univariateModels <- list();
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]
  id <- NULL
  if ( is.null(target) )  {
    y <- dataset[, 1]
    x <- bigmemory::sub.big.matrix(dataset, firstCol = 2)
  } else {
    y <- target
    x <- dataset
  }

  ## Beta regression 
  if ( identical(test, testIndBeta) ) {
    
    param <- Rfast::beta.mle(y)$param
    m1 <- digamma(param[1]) - digamma(param[2])
    z <- log(y) - log(1 - y)
    u <- colSums(x[] * (z - m1))
    m2 <- trigamma(param[1]) + trigamma(param[2])
    seu <- colSums(x[]^2) * m2
    stat <- u^2/seu
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
    
    ## Negative Binomial 
  } else if ( identical(test, testIndNB) ) {
    mod <- Rfast::negbin.mle(y)
    r <- mod$param[2]
    p <- mod$param[1]
    my <- mod$param[3]
    sxy <- colSums(y * x[])
    u <- sxy - (1 - p) * ( sxy + r * colSums( x[]) )
    vu <- colSums(x[]^2) * p^2 * (my + my^2/r)
    stat <- u^2/vu
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
    
    ## Poisson
  } else if ( identical(test, testIndPois) ) {
    n <- length(y)
    sx2 <- colSums(x[]^2)
    my <- sum(y)/n
    sx <- colSums(x[])
    up <- as.vector( cov(y, x[]) ) * (n - 1)
    down <- (sx2 - sx^2/n) * my
    stat <- up/sqrt(down)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pvalue <- log(2) + pt(abs(stat), n - 2, lower.tail = FALSE, log.p = TRUE)
    
    ## logistic regression
  } else if ( identical(test, testIndLogistic) ) { 
    n <- length(y)
    sx2 <- colSums(x[]^2)
    my <- sum(y)/n
    sx <- colSums(x[])
    up <- as.vector( cov(y, x[]) ) * (n - 1)
    down <- (sx2 - sx^2/n) * (my - my^2)
    stat <- up/sqrt(down)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pvalue <- log(2) + pt(abs(stat), n - 2, lower.tail = FALSE, log.p = TRUE)
    
    ## multinomial regression
  } else if ( identical(test, testIndMultinom) ) { 
    n <- length(y)
    p <- dim(x)[2]
    dof <- Rfast::sort_unique.length(y) - 1
    if (dof == 1) {
      sx2 <- colSums(x[]^2)
      my <- sum(y)/n
      sx <- colSums(x[])
      up <- as.vector( cov(y, x[]) ) * (n - 1)
      down <- (sx2 - sx^2/n) * (my - my^2)
      stat <- up/sqrt(down)
      pvalue <- log(2) + pt(abs(stat), n - 2, lower.tail = FALSE, log.p = TRUE)
    } else {
      m0 <- numeric(dof)
      y1 <- Rfast::design_matrix(y)[, -1]
      m <- Rfast::colmeans(y1)
      sx <- colSums(x[])
      sx2 <- colSums(x[]^2)
      vp <- diag(m) - tcrossprod(m)
      mx <- matrix(rep(m, rep(p, dof)), ncol = dof)
      ni <- tabulate(y)
      u <- t(rowsum(x[], y))[, -1] - sx * mx
      stat <- Rfast::mahala(u, m0, vp)/sx2
      pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    }
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
    
    ## Gamma regression
  } else if ( identical(test, testIndGamma) ) { 
    pa <- Rfast::gammamle(y)$param
    m <- pa[1]/pa[2]
    u <- colSums(x[]) - colSums(y * x[])/m
    vb <- colSums(x[]^2)/pa[1]
    stat <- u^2/vb
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
    
    # Weibull
  } else if ( identical(test, censIndWR) ) {
    mod <- Rfast::weibull.mle(y)
    k <- mod$param[1]
    lam <- mod$param[2]
    yk <- y^k
    u <- k/lam^k * colSums(yk * x[]) - k * colSums(x[])
    vu <- k^2 * colSums(x[]^2)
    stat <- u^2/vu
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
    
  } else univariateModels <- NULL
  
  if ( !is.null(univariateModels) )  {
    id <- which( is.na(univariateModels$stat) )
    univariateModels$stat[id] <- 0
    univariateModels$pvalue[id] <- 1
  }
  
  univariateModels
  
}