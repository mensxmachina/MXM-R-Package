zip.regs <- function(target, dataset, wei = NULL, check = FALSE, logged = FALSE, ncores = 1) {
  
  dm <- dim(dataset)
  if (check) {
    id <- Rfast::check_data(dataset)
    if ( sum(id>0) > 0 )  dataset[, id] <- rnorm(n * length(id) )
  }
  
  if ( ncores <= 1 ) {
    oop <- options(warn = -1) 
	on.exit( options(oop) ) 
    D <- ncol(dataset)
    n <- length(target)
    dof <- loglik <- numeric(D)   
    poia <- which( target == 0 ) 
    n0 <- length(poia)    ;     n1 <-  n - n0
    target1 <- target[ -poia ]   

    if ( is.null(wei) ) {
      lgy <- sum( lgamma(target1 + 1) )  
      ini <- 2 * Rfast::zip.mle(target)$loglik

      for (i in 1:D) { 
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )
        mod <- glm.fit(x[-poia, ], target1, family = poisson(log) ) 
        p1 <- ( n0 - sum( exp( - mod$fitted.values ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) ) 
        pa <- c( g1, mod$coefficients)
        pa[is.na(pa)] <- rnorm( sum(is.na(pa)) )
        lik <- nlm( regzip, pa, y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
        lik2 <- optim( lik$estimate, regzip, y1 = target1, x = x, n1 = n1, poia = poia, control = (maxit = 10000) )
        loglik[i] <-  - lik2$value 
        dof[i] <- length(lik2$par) - 2
      }
    } else { 
      wei <- wei / sum(wei)
      w0 <- wei[poia]    ;   w1 <- wei[-poia] 
      lgy <- sum( w1 * lgamma(target1 + 1) )  
      ini <- 2 * zipmle.wei(target, wei)$loglik

      for (i in 1:D) {
        mod <- glm.fit(x[-poia, ], target1, family = poisson(log), weights = w1 )
        p1 <- ( n0 - sum( exp( - mod$fitted.values ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) )
        pa <- c( g1, mod$coefficients)
        pa[is.na(pa)] <- rnorm( sum(is.na(pa)) )
        lik <- nlm( regzipwei, pa, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, poia = poia, iterlim = 10000 )
        lik2 <- optim( lik$estimate, regzipwei, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, poia = poia, control = list(maxit = 10000) )  
        dof[i] <- length(lik2$par) - 2
        loglik[i] <-  - lik2$value
      }
    }
    
    bic <-  - 2 *loglik + 2 * lgy + (dof + 2) * log(n)
    stat <- 2 * loglik - 2 * lgy - ini
    pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = logged) 
    
  } else {
    
    if ( is.null(wei) ) {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      poia <- which( target == 0 ) 
      D <- ncol(dataset)
      n <- length(target)
      n0 <- length(poia)    
      n1 <-  n - n0
      target1 <- target[ -poia ]   
      lgy <- sum( lgamma(target1 + 1) )  
      ini <- 2 * Rfast::zip.mle(target)$loglik
      oop <- options(warn = -1) 
	  on.exit( options(oop) ) 
	  
      mod <- foreach::foreach( i = 1:D, .combine = rbind, .export = "regzip") %dopar% {
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )
        mod <- glm.fit(x[-poia, ], target1, family = poisson(log) )
        p1 <- ( n0 - sum( exp( - mod$fitted.values ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) )  
        pa <- c( g1, mod$coefficients)
        pa[is.na(pa)] <- rnorm( sum(is.na(pa)) )
        lik <- nlm( regzip, pa, y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
        lik2 <- optim( lik$estimate, regzip, y1 = target1, x = x, n1 = n1, poia = poia, control = list(maxit = 10000) )
        return( c( - lik2$value, length(lik2$par) - 2 ) )
      }
      parallel::stopCluster(cl)

    } else {  
      cl <-  parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      poia <- which( target == 0 ) 
      D <- ncol(dataset)
      n <- length(target)
      n0 <- length(poia)    
      n1 <-  n - n0
      target1 <- target[ -poia ]   
      wei <- wei / sum(wei)
      w0 <- wei[poia]    ;   w1 <- wei[-poia] 
      ini <- 2 * zipmle.wei(target, wei)$loglik
      lgy <- sum( w1 * lgamma(target1 + 1) )  
      oop <- options(warn = -1) 
	  on.exit( options(oop) ) 
	  
      mod <- foreach::foreach( i = 1:D, .combine = rbind, .export = "regzipwei" ) %dopar% {
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )
        mod <- glm.fit(x[-poia, ], target1, family = poisson(log), weights = w1 )
        p1 <- ( n0 - sum( exp( - mod$fitted.values ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) )
        pa <- c( g1, mod$coefficients)
        pa[is.na(pa)] <- rnorm( sum(is.na(pa)) )
        lik <- nlm( regzipwei, pa, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, iterlim = 10000 )
        lik2 <- optim( lik$estimate, regzipwei, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, control = list(maxit = 10000) )  
        return( c( - lik2$value, length(lik2$par) - 2 ) )
      }
      parallel::stopCluster(cl)
    }
    
    bic <-  - 2 * mod[, 1] +2 * lgy + (mod[, 2] + 2) * log(n) 
    stat <- 2 * mod[, 1] - 2 * lgy - ini
    pvalue <- pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = logged) 
  }  

  if (check) {
    if ( sum(id>0) > 0 ) {
      stat[id] <- 0
      pvalue[id] <- log(1)
      bic[id] <- NA
    }
  }  
  cbind(stat, pvalue, bic)
}




#zip.regs(y, dataset, logged = T, wei = NULL, ncores = 2)

#ela <- function(y, x) {
#  ini <- 2 * as.numeric( logLik( zeroinfl(y~1|1) ) )
#  p <- ncol(x)
#  stat <- numeric(p)
#  for (i in 1:p) { 
#    mod <- zeroinfl(y~x[, i]|1)
#    stat[i] <- as.numeric( logLik(mod) )
#  }
  
#  stat <- 2 * stat - ini
#  pval <- pchisq(stat, 1, lower.tail=F, log.p = T)
#  cbind(stat, pval)
#}
  
  

