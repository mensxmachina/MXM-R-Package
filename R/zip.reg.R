zip.reg <- function(target, dataset, wei = NULL, lgy = NULL) {
  
  n <- length(target)
  oop <- options(warn = -1) 
  on.exit( options(oop) ) 
  if ( NCOL(dataset) == 0 ) {
    if ( is.null(wei) ) {
      mod <- Rfast::zip.mle(target)
      lam <- mod$param[1]
      prop <- mod$param[2]
      res <- list(be = lam, prop = prop, loglik = mod$loglik, est = (1 - prop) * lam)
    } else {
      mod <- zipmle.wei(target, wei)
      res <- list(be = mod$lam, prop = mod$prop, loglik = mod$loglik, est = (1 - mod$prop) * mod$lam)
    }
    
  } else {
    x <- model.matrix(target ~ ., as.data.frame(dataset) )
    poia <- which(target == 0)
    n0 <- length(poia)   ;    n1 <-  n - n0
    target1 <- target[ -poia ]    
    
    if ( is.null(wei) ) {
      mod <- glm.fit(x[-poia, ], target1, family = poisson(log) )
      p1 <- ( n0 - sum( exp( - mod$fitted.values ) ) ) / n
      g1 <- log( p1 / ( 1 - p1 ) )
      pa <- c( g1, mod$coefficients)
      pa[is.na(pa)] <- rnorm( sum(is.na(pa)) )
      lik <- nlm( regzip, pa, y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
      lik2 <- optim( lik$estimate, regzip, y1 = target1, x = x, n1 = n1, poia = poia, control = list(maxit = 10000) )
      if ( is.null(lgy) )  lgy <- sum( lgamma(target1 + 1) )  
     
    } else {  
      wei <- wei / sum(wei)
      w0 <- wei[poia]     ;    w1 <- wei[-poia]
      mod <- glm.fit(x[-poia, ], target1, family = poisson(log), weights = w1 )
      p1 <- ( n0 - sum( exp( - mod$fitted.values ) ) ) / n
      g1 <- log( p1 / ( 1 - p1 ) )
      pa <- c( g1, mod$coefficients)
      pa[is.na(pa)] <- rnorm( sum(is.na(pa)) )
      lik <- nlm( regzipwei, pa, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, iterlim = 10000 )
      lik2 <- optim( lik$estimate, regzipwei, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, control = list(maxit = 10000) )  
      if ( is.null(lgy) )  lgy <- sum( w1 * lgamma(target1 + 1) )  
    }
    prop <- exp(lik2$par[1]) / ( 1 + exp(lik2$par[1]) )
    res <- list(be = lik2$par[-1], prop = prop, loglik = -lik2$value - lgy)
  }  
  res
}