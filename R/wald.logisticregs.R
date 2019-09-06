wald.logisticregs <- function(y, x, tol = 1e-09, wei = NULL, check = FALSE, logged = FALSE, ncores = 1) { 
  
  if ( is.factor(y) )   y <- as.numeric(y) - 1
  if (check) {
    id <- Rfast::check_data(x, ina = y)
    if ( sum(id>0) > 0 )  x[, id] <- rnorm(length(y) * length(id) )
  }  
  
  if ( is.null(wei) ) {
    
  if ( ncores <= 1 ) { 

    dm <- dim(x)    
    sy <- sum(y)
    n <- dm[1]
    p0 <- sy / n
    y0 <- 1 - y
    D <- dm[2]
    stat <- numeric(D) 
    a0 <- log( p0 / (1 - p0) )
    p0 <- rep(p0, n)
    w0 <- p0 * (1 - p0)
    de <- y - p0
    dera20 <- sum(w0)
    dera0 <- sy - n * p0[1]
    
    for (j in 1:D) {
      X <- x[, j] 
      derb <- sum(de * X) 
      derab <- sum(w0 * X)
      derb2 <- sum(w0 * X^2)
      bold <- c(a0, 0)
      bnew <- bold + c( derb2 * dera0 - derab * derb, - derab * dera0 + dera20 * derb ) / ( dera20 * derb2 - derab^2 )
        
      while ( sum( abs(bnew - bold) )  > tol ) {
        bold <- bnew
        a <-  - bold[1]   ;  b <-  - bold[2]
        m <- exp( + a + b * X)
        p <- 1 / (1 + m)
        de <- y - p
        dera <- sy - sum(p)
        derb <- sum(de * X) 
        w <- p * (1 - p)
        dera2 <- sum(w)
        derab <- sum(w * X)
        derb2 <- sum(w * X^2)
        bnew <- bold + c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
      }        
      stat[j] <- bnew[2]^2 / dera2 * ( dera2 * derb2 - derab^2 )   
    }
    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    res <- cbind(stat, pvalue)
    
  } else {
    
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)

    dm <- dim(x)    
    sy <- sum(y)
    n <- dm[1]
    p0 <- sy / n
    y0 <- 1 - y
    D <- dm[2]
    stat <- numeric(D) 
    a0 <- log( p0 / (1 - p0) )
    p0 <- rep(p0, n)
    w0 <- p0 * (1 - p0)
    de <- y - p0
    dera20 <- sum(w0)
    dera0 <- sy - n * p0[1]
    
    stat <- foreach( j = 1:D, .combine = rbind ) %dopar% {
      X <- x[, j] 
      derb <- sum(de * X) 
      derab <- sum(w0 * X)
      derb2 <- sum(w0 * X^2)
      bold <- c(a0, 0)
      bnew <- bold + c( derb2 * dera0 - derab * derb, - derab * dera0 + dera20 * derb ) / ( dera20 * derb2 - derab^2 )
      while ( sum( abs(bnew - bold) )  > tol ) {
        bold <- bnew
        a <-  - bold[1]   ;  b <-  - bold[2]
        m <- exp( + a + b * X)
        p <- 1 / (1 + m)
        de <- y - p
        dera <- sy - sum(p)
        derb <- sum(de * X) 
        w <- p * (1 - p)
        dera2 <- sum(w)
        derab <- sum(w * X)
        derb2 <- sum(w * X^2)
        bnew <- bold + c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
      }        
      ta <- bnew[2]^2 / dera2 * ( dera2 * derb2 - derab^2 )   
      return( ta )
    }
    stopCluster(cl)
    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    res <- cbind(stat, pvalue)
  }
    
  } else {
    if (ncores <= 1) {
      D <- ncol(x)
      stat <- numeric(D)
      for (i in 1:D) {
        mod <- glm(y ~ x[, i], binomial, weights = wei)
        stat[i] <- summary(mod)[[ 12 ]][2, 3]^2
      }
    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      stat <- foreach( j = 1:D, .combine = rbind ) %dopar% {
        mod <- glm(y ~ x[, i], binomial, weights = wei)
        return( summary(mod)[[ 12 ]][2, 3]^2)
      }    
    }
    pvalue <- pchisq(stat, 1, lower.tail = FALSE)
    res <- cbind(stat, pvalue)
     
  }  
  
  res
}



#ela2 <- function(y, x) {
#  D <- ncol(x)
#  stat <- dof <- bic <- numeric(D)
#  for (i in 1:D) {
#    mod = glm(y ~ x[, i], binomial)
#    stat[i] <- mod$deviance
#    dof[i] <- length( coef(mod) ) - 1
#    bic[i] <- BIC(mod)
#  }
#    stat <- glm(y ~ 1, binomial)$null.dev - stat
#    pvalue <- pchisq(stat, dof, lower.tail = FALSE)
#    cbind(stat, pvalue, bic)
#}

