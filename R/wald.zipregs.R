wald.zipregs <- function(target, dataset, wei = NULL, check = FALSE, logged = FALSE, ncores = 1) {
  
  if (check) {
    id <- which( Rfast::colrange(dataset) == 0 )
    if ( length(id) > 0 )  dataset[, id] <- rnorm( length(target) * length(id) )
  }
  if ( ncores <= 1 ) {
    D <- ncol(dataset)
    n <- length(target)
    stat <- numeric(D)   
    poia <- which( target == 0 ) 
    n0 <- length(poia)    ;     n1 <-  n - n0
    target1 <- target[ -poia ]   

    if ( is.null(wei) ) {
      lgy <- sum( lgamma(target1 + 1) )  

      for (i in 1:D) { 
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )[1:n, ]
        mod <- glm.fit(x[ - poia, ], target1, family = poisson(log) ) 
        p1 <- ( n0 - sum( exp( - fitted(mod) ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) ) 
        lik <- nlm( regzip, c(g1, coef(mod) ), y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
        lik2 <- optim(lik$estimate, regzip, y1 = target1, x = x, n1 = n1, poia = poia, control = list(maxit = 10000), hessian = TRUE )
        s <- diag( solve( lik2$hessian ) )
        stat[i] <- lik2$par^2/s[3]
      }
    } else { 
      wei <- wei / sum(wei)
      w0 <- wei[poia]    ;   w1 <- wei[-poia] 
      lgy <- sum( lgamma(target1 + 1) )  
      ini <- 2 * zipmle.wei(target, wei)$loglik

      for (i in 1:D) {
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )[1:n, ]
        mod <- glm.fit( x[ -poia, ], target1, family = poisson(log), weights = w1 ) 
        p1 <- ( n0 - sum( exp( - fitted(mod) ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) )
        lik <- nlm( regzipwei, c(g1, coef(mod) ), y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, poia = poia, iterlim = 10000 )
        lik2 <- optim( lik$estimate,regzipwei, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, poia = poia, control = list(maxit = 10000), hessian = TRUE)  
        s <- diag( solve( lik2$hessian ) )
        stat[i] <- lik2$par^2/s[3]
      }
    }

    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged) 
    
  } else {
    
    if ( is.null(wei) ) {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      poia <- which( target == 0 ) 
      D <- ncol(dataset)
      n <- length(target)
      n0 <- length(poia)    
      n1 <-  n - n0
      target1 <- target[ -poia ]   
      lgy <- sum( lgamma(target1 + 1) )  
      ini <- 2 * Rfast::zip.mle(target)[3]

      mod <- foreach( i = 1:D, .combine = rbind, .export = "regzip") %dopar% {
        x <- cbind(1, dataset[, i]) 
        mod2 <- glm.fit(x[ - poia, ], target1, family = poisson(log) ) 
        p1 <- ( n0 - sum( exp( - fitted(mod2) ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) )  
        lik <- nlm( regzip, c(g1, coef(mod2) ), y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
        lik2 <- optim(lik$estimate, regzip, y1 = target1, x = x, n1 = n1, poia = poia, control = list(maxit = 10000), hessian = TRUE )
        s <- diag( solve( lik2$hessian ) )
        return( lik2$par^2/s[3] )
      }

      stopCluster(cl)

    } else {  
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      poia <- which( target == 0 ) 
      D <- ncol(dataset)
      n <- length(target)
      n0 <- length(poia)    
      n1 <-  n - n0
      target1 <- target[ -poia ]   
      wei <- wei / sum(wei)
      w0 <- wei[poia]    ;   w1 <- wei[-poia] 
      lgy <- sum( w1 * lgamma(target1 + 1) )  
         
      mod <- foreach( i = 1:D, .combine = rbind, .export = "regzipwei" ) %dopar% {
        x <- cbind(1, dataset[, i]) 
        mod2 <- glm.fit( x[ -poia, ], target1, family = poisson(log), weights = w1 ) 
        p1 <- ( n0 - sum( exp( - fitted(mod2) ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) )
        lik <- nlm( regzipwei, c(g1, coef(mod2) ), y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, iterlim = 10000 )
        lik2 <- optim( lik$estimate, regzipwei, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, control = list(maxit = 10000), hessian = TRUE )  
        s <- diag( solve( lik2$hessian ) )
        return( lik2$par^2/s[3] )
      }

    }

    pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = logged) 
  }  

  if ( length(id) > 0 ) {
    stat[id] = 0
    pvalue[id] = log(1)
    bic[id] = NA
  }  
  cbind(stat, pvalue)
}





   
   