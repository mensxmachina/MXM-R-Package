wald.betaregs <- function(target, dataset, wei = NULL, check = FALSE, logged = FALSE, ncores = 1) {
  
  dm <- dim(dataset) 
  D <- dm[2]
  if (check) {
    id <- which( Rfast::colrange(dataset) == 0 )
    if ( length(id) > 0 )  dataset[, id] <- rnorm( length(target) * length(id) )
  }
  
  if ( ncores <= 1 ) {
    oop <- options(warn = -1) 
    on.exit( options(oop) )
    if ( is.null(wei) ) {
      ly1 <- log(1 - target)     ;    ly <- log(target) - ly1           
      sly1 <- sum(ly1)           ;    sly2 <- sum( log( target ) ) + sly1  
    } else {
      w <- wei / sum(wei)
      ly1 <- w * log(1 - target)     ;    ly <- w * log(target) - ly1
      sly1 <- sum(ly1)               ;    sly2 <- sum( w * log( target * (1 - target) ) ) 
    }
    
    n <- dm[1]
    a <- Rfast::beta.mle(target)
    iniphi <- log( sum(a$param) )
    stat <- numeric(D)
    id <- which( Rfast::colVars(dataset) == 0 )
    if ( length(id) > 0 )  dataset[, id] <- target
    
    if ( is.null(wei) ) {
      for (i in 1:D) {
        x <- model.matrix(y~., data.frame( dataset[, i]) )
        mod1 <- nlm(regbeta, c(iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
        mod2 <- optim(mod1$estimate, regbeta, ly = ly, sly1 = sly1, x = x, n = n, control = list(maxit = 10000), hessian = TRUE )
        s <- diag( solve(mod2$hessian) ) 
        stat[i] <- mod2$par[3]^2/s[3]
      }

    } else {
      for (i in 1:D) {
        x <- model.matrix(y~., data.frame( dataset[, i]) )
        mod1 <- nlm(regbetawei, c(iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
        mod2 <- optim(mod1$estimate, regbetawei, ly = ly, sly1 = sly1, x = x, w = w, control = list(maxit= 10000), hessian = TRUE )
        s <- diag( solve(mod2$hessian) ) 
        stat[i] <- mod2$par[3]^2/s[3]
      }
    }
    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged) 
    
  } else {

   if ( is.null(wei) ) {
     cl <- parallel::makePSOCKcluster(ncores)
     doParallel::registerDoParallel(cl)
     ly1 <- log(1 - target)     ;    ly <- log(target) - ly1           
     sly1 <- sum(ly1)           ;    sly2 <- sum( log( target ) ) + sly1  
     n <- length(target)
     a <- Rfast::beta.mle(target)
     iniphi <- log( sum(a$param) )
     oop <- options(warn = -1) 
     on.exit( options(oop) )
     mod <- foreach::foreach( i = 1:D, .combine = rbind, .export = "regbeta" ) %dopar% {
       x <- model.matrix(y~., data.frame( dataset[, i]) )
       mod1 <- nlm(regbeta, c(iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
       mod2 <- optim(mod1$estimate, regbeta, ly = ly, sly1 = sly1, x = x, n = n, control = list(maxit = 10000), hessian = TRUE)
       s <- diag( solve(mod2$hessian) ) 
       return( mod2$par[3]^2/s[3] )
     }
     parallel::stopCluster(cl)

   } else {
     cl <- parallel::makePSOCKcluster(ncores)
     doParallel::registerDoParallel(cl)
     n <- length(target)
     w <- wei / sum(wei)
     ly1 <- w * log(1 - target)     ;    ly <- w * log(target) - ly1
     sly1 <- sum(ly1)               ;    sly2 <- sum( w * log( target * (1 - target) ) ) 
     a <- Rfast::beta.mle(target)
     iniphi <- log( sum( a$param ) )
     oop <- options(warn = -1) 
     on.exit( options(oop) )
     mod <- foreach::foreach( i = 1:D, .combine = rbind, .export = "regbetawei" ) %dopar% {
        x <- model.matrix(y~., data.frame( dataset[, i]) ) 
        mod1 <- nlm(regbetawei, c(iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
        mod2 <- optim(mod1$estimate, regbetawei, ly = ly, sly1 = sly1, x = x, w = w, control = list(maxit= 10000), hessian = TRUE )
        s <- diag( solve(mod2$hessian) ) 
        return( mod2$par[3]^2/s[3] )
     }
     parallel::stopCluster(cl)
   }
    
    pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = logged) 
  } 
  
  if ( length(id) > 0 ) {
    stat[id] <- 0
    pvalue[id] <- log(1)
    bic[id] <- NA
  }  
  cbind(stat, pvalue)
}





