perm.betaregs <- function(target, dataset, wei = NULL, check = FALSE, logged = FALSE, threshold = 0.05, R = 999, ncores = 1) {
  
  dm <- dim(dataset)
  n <- dm[1]
  D <- dm[2]
  if (check) {
    id <- Rfast::check_data(dataset)
    if ( sum(id>0) > 0 )  dataset[, id] <- rnorm(n * length(id) )
  }
  thres <- threshold * R + 1

  if ( ncores <= 1 ) {
    oop <- options(warn = -1) 
    on.exit( options(oop) )    
    if ( is.null(wei) ) {
      ly1 <- log(1 - target)     ;    ly <- log(target) - ly1           
      sly1 <- sum(ly1)           ;    sly2 <- sum( log( target ) ) + sly1  
      a <- Rfast::beta.mle(target)
      ini <- 2 * a$loglik
    } else {
      w <- wei / sum(wei)
      ly1 <- w * log(1 - target)     ;    ly <- w * log(target) - ly1
      sly1 <- sum(ly1)               ;    sly2 <- sum( w * log( target * (1 - target) ) ) 
      a <- betamle.wei(target, w)    
      ini <- 2 * a$loglik
    }
    iniphi <- log( sum(a$param) )
    pvalue <- dof <- loglik <- numeric(D)   
    
    if ( is.null(wei) ) {
      for (i in 1:D) {
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )
        mod1 <- nlm(regbeta, c(iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
        mod2 <- nlm(regbeta, mod1$estimate, ly = ly, sly1 = sly1, x = x, n = n,  iterlim = 10000 )
        loglik[i] <- mod2$minimum 
        dof[i] <- length(mod2$estimate)
        step <- 0
        j <- 1		
        while (j <= R & step < thres ) {
          xb <- x[sample(n, n), ]  
          bit2 <- nlm(regbeta, mod2$estimate, ly = ly, sly1 = sly1, x = xb, n = n,  iterlim = 10000 )
          step <- step + ( bit2$minimum < loglik[i] )
          j <- j + 1
        }
        pvalue[i] <- (step + 1) / (R + 1) 
      }
    } else {

      for (i in 1:D) {
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )
        mod1 <- nlm(regbetawei, c(iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
        mod2 <- nlm(regbetawei, mod1$estimate, ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
        loglik[i] <-  mod2$minimum 
        dof[i] <- length(mod2$estimate)
        step <- 0
        j <- 1		
        while (j <= R & step < thres ) {
          xb <- x[sample(n, n), ]  
          bit2 <- nlm(regbeta, mod2$estimate, ly = ly, sly1 = sly1, x = xb, n = n,  iterlim = 10000 )
          step <- step + ( bit2$minimum < loglik[i] )
          j <- j + 1
        }
        pvalue[i] <- (step + 1) / (R + 1) 
      }
    }
    lik <-  - loglik - sly2
    bic <-  - 2 * lik + dof * log(n) 
    stat <- 2 * lik - ini

  } else {
    
    if ( is.null(wei) ) {
      
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      ly1 <- log(1 - target)     ;    ly <- log(target) - ly1           
      sly1 <- sum(ly1)           ;    sly2 <- sum( log( target ) ) + sly1  
      a <- Rfast::beta.mle(target)
      ini <- 2 * a$loglik
      n <- length(target)
      iniphi <- log( sum(a$param) )
      oop <- options(warn = -1) 
      on.exit( options(oop) )
      mod <- foreach::foreach( i = 1:D, .combine = rbind, .export = "regbeta" ) %dopar% {
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )
        mod1 <- nlm(regbeta, c(iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
        mod2 <- nlm(regbeta, mod1$estimate, ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
        loglik <- mod2$minimum
        dof <- length(mod2$estimate)
        step <- 0
        j <- 1		
        while (j <= R & step < thres ) {
          xb <- x[sample(n, n), ]  
          bit2 <- nlm(regbeta, mod2$estimate, ly = ly, sly1 = sly1, x = xb, n = n,  iterlim = 10000 )
          step <- step + ( bit2$minimum < loglik )
          j <- j + 1
        }
        pvalue <- (step + 1) / (R + 1) 
        return( c( loglik, pvalue, dof ) )
      }
      parallel::stopCluster(cl)
    } else {
      
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      n <- length(target)
      w <- wei / sum(wei)
      ly1 <- w * log(1 - target)     ;    ly <- w * log(target) - ly1
      sly1 <- sum(ly1)               ;    sly2 <- sum( w * log( target * (1 - target) ) )    
      a <- betamle.wei(target ,w)
      ini <- 2 * a$loglik
      iniphi <- log( sum( a$param ) )
      oop <- options(warn = -1) 
      on.exit( options(oop) )      
      mod <- foreach::foreach( i = 1:D, .combine = rbind, .export = "regbetawei" ) %dopar% {
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )
        mod1 <- nlm(regbetawei, c(iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
        mod2 <- nlm(regbetawei, mod1$estimate, ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
        loglik <- mod2$minimum
        dof <- length(mod2$estimate)
        step <- 0
        j <- 1		
        while (j <= R & step < thres ) {
          xb <- x[sample(n, n), ]  
          bit2 <- nlm(regbetawei, mod1$estimate, ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
          step <- step + ( bit2$minimum < loglik )
          j <- j + 1
        }
        pvalue <- (step + 1) / (R + 1) 
        return( c( loglik, pvalue, dof ) )
      }
      parallel::stopCluster(cl)
    }
    
    lik <-  - mod[, 1] - sly2
    bic <-  - 2 * lik + mod[, 3] * log(n)
    stat <-  2 * lik - ini
    pvalue <- mod[, 2]
  }  
  
  if (check) {
    if ( sum(id>0) > 0 ) {
      stat[id] <- 0
      pvalue[id] <- log(1)
      bic[id] <- NA
    }
  } 
  cbind(stat, log(pvalue), bic)
}




#ela = function(y, x) {
#  ini = as.numeric( logLik( betareg(y ~ 1) ) )
#  D <- dim(x)[2]
#  lik = dof = numeric(D)
#  for (i in 1:D) {
#    mod <- betareg(y ~ x[, i])
#    lik[i] = as.numeric( logLik(mod) )
#    dof[i] = length( coef(mod) )
#  } 
#  stat <- 2 * (lik - ini)
#  pvalue <- pchisq(stat, dof - 2, lower.tail =F)
#  cbind(stat, pvalue)
#}