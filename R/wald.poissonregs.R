wald.poissonregs <- function(y, x, tol = 1e-09, wei = NULL, check = FALSE, logged = FALSE, ncores = 1) { 
   
  if (check) {
    id <- Rfast::check_data(x)
    if ( sum(id>0) > 0 )  x[, id] <- rnorm(length(y) * length(id) )
  }
  
  if ( is.null(wei) ) {
    
  if ( ncores <= 1 ) {
    
    dm <- dim(x)
    sy <- sum(y)
    n <- dm[1]
    m0 <- sy / n
    a0 <- log( m0 )
    de <- y - m0
    D <- dm[2]
    stat <- numeric(D) 
    dera20 <- sy
    dera0 <- sy - dera20
 
    for (j in 1:D) { 
      X <-  x[, j]
      derb <- sum(de * X) 
      derab <- sum(m0 * X)
      derb2 <- sum(m0 * X^2)
      bold <- c(a0, 0)
      bnew <- bold + c( derb2 * dera0 - derab * derb, - derab * dera0 + dera20 * derb ) / ( dera20 * derb2 - derab^2 )
      while ( sum( abs(bnew - bold) )  > tol ) {
        bold <- bnew
        a <- bold[1]   ;   b <- bold[2]
        m <- exp(a + b * X)
        de <- y - m
        dera <- sy - sum(m)
        derb <- sum(de * X) 
        dera2 <- sum(m)
        derab <- sum(m * X)
        derb2 <- sum(m * X^2)
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
    m0 <- sy / n
    a0 <- log( m0 )
    de <- y - m0
    D <- dm[2]
    stat <- numeric(D) 
    dera20 <- sy
    dera0 <- sy - dera20
  
    stat <- foreach( j = 1:D, .combine = rbind ) %dopar% {
      X <-  x[, j]
      derb <- sum(de * X) 
      derab <- sum(m0 * X)
      derb2 <- sum(m0 * X^2)
      bold <- c(a0, 0)
      bnew <- bold + c( derb2 * dera0 - derab * derb, - derab * dera0 + dera20 * derb ) / ( dera20 * derb2 - derab^2 )
      while ( sum( abs(bnew - bold) )  > tol ) {
        bold <- bnew
        a <- bold[1]   ;   b <- bold[2]
        m <- exp(a + b * X)
        de <- y - m
        dera <- sy - sum(m)
        derb <- sum(de * X) 
        dera2 <- sum(m)
        derab <- sum(m * X)
        derb2 <- sum(m * X^2)
        bnew <- bold + c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
      }
      ta <- bnew[2]^2 / dera2 * ( dera2 * derb2 - derab^2 )   
      return(ta)
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
      mod <- glm(y ~ x[, i], poisson, weights = wei)
      stat[i] <- summary(mod)[[ 12 ]][2, 3]^2
    }
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    stat <- foreach( j = 1:D, .combine = rbind ) %dopar% {
      mod <- glm(y ~ x[, i], poisson, weights = wei)
      return( summary(mod)[[ 12 ]][2, 3]^2)
    }    
  }
  pvalue <- pchisq(stat, 1, lower.tail = FALSE)
  res <- cbind(stat, pvalue)
  
}  
  
  res
}


# ela3 <- function(y, x) {
#   D <- ncol(x)
#   stat <- dof <- bic <- numeric(D)
#   for (i in 1:D) {
#     mod = glm(y ~ x[, i], poisson)
#     stat[i] <- mod$deviance
#     dof[i] <- length( coef(mod) ) - 1
#     bic[i] <- BIC(mod)
#   }
#     stat <- glm(y ~ 1, poisson)$null.dev - stat
#     pvalue <- pchisq(stat, dof, lower.tail = FALSE)
#     cbind(stat, pvalue, bic)
# }

