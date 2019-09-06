#######################
### K-fold cross validation for selecting the optimal lambda in ridge regression
#######################
 
### usage:  ridgereg.cv( target, dataset, K = 10, lambda = seq(0, 1, by = 0.01), 
### auto = FALSE, seed = FALSE, ncores = 1, mat = NULL )

ridgereg.cv <- function( target, dataset, K = 10, lambda = seq(0, 2, by = 0.1), 
                        auto = FALSE, seed = FALSE, ncores = 1, mat = NULL ) {
  ## target is a dependent continuous variable or a matrix with continuous variables
  ## dataset is a matrix with continuous variables only
  ## K is the number of folds for the K-fold cross validation, set to 10 by default
  ## lambda is a vector containing a grid of values
  ## auto is a boolean variable. If TRUE, the GCV criterion is used to return the best lambda automatically.
  ## Otherwise, a K-fold cross validation is performed
  ## seed is boolean. Should the same K-folds be used always or not?
  ## ncores specifies how many cores to be used
  target <- as.vector(target)
  dataset <- as.matrix(dataset)
  n <- length(target)
  p <- ncol(dataset)
  mspe <- NULL
  performance <- NULL

  if ( auto ) {
    runtime <- proc.time()
    mod <- MASS::lm.ridge( target ~ dataset, lambda = lambda ) 
    gcv <- min(mod$GCV)
    plot(lambda, mod$GCV, type = "b", xlab = expression(paste(lambda, " values")), ylab = "GCV")
    lam <- lambda[ which.min(mod$GCV) ]
    runtime <- proc.time() - runtime

  } else {
    if ( is.null(mat) ) { ## no folds were given by the user
      if ( seed )  set.seed(1234567) ## the folds will always be the same
      nu <- sample(1:n, min( n, round(n / K) * K ) )
      ## It may be the case this new nu is not exactly the same
      ## as the one specified by the user
      oop <- options(warn = -1)   # if the length of nu does not fit 
      ## to a matrix a warning message should appear 
      on.exit( options(oop) )
      mat <- matrix( nu, ncol = K ) 
    } else mat <- mat
      rmat <- dim(mat)[1]
      
    if ( ncores == 1 ) {
      
      runtime <- proc.time()
      mi <- length(lambda)
      per <- matrix( nrow = K, ncol = mi )
      
      for (vim in 1:K) {
        ytest <- as.vector( target[ mat[, vim] ] )  ## test set dependent vars
        ytrain <- as.vector( target[ -mat[, vim] ] )  ## train set dependent vars
        my <- mean(ytrain) 
		    yy <- ytrain - my
        xtrain <- dataset[ -mat[, vim], , drop = FALSE]  ## train set independent vars
        xtest <- dataset[ mat[, vim], , drop = FALSE]  ## test set independent vars
        sa <- svd(xtrain)
        d <- sa$d    ;    v <- t(sa$v)    ;     tu <- t(sa$u)  
        d2 <- d^2    ;    A <- d * tu %*% yy
        
        for (i in 1:mi) {
          ## betas <- ( v %*% (tu * ( d / ( d^2 + lambda[i] ) ) ) ) %*% yy 
          betas <- crossprod( v / ( d2 + lambda[i] ), A )
          est <- xtest %*% betas + my 
          per[vim, i] <- sum( (ytest - est)^2 ) / rmat
        }
        
      }
      
      runtime <- proc.time() - runtime

    } else {
      
      runtime <- proc.time()
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mi <- length(lambda)
      pe <- numeric(mi)
	  
      per <- foreach::foreach(vim = 1:K, .combine = rbind) %dopar% {
        ytest <- as.vector( target[mat[, vim] ] )  ## test set dependent vars
        ytrain <- as.vector( target[-mat[, vim] ] )  ## train set dependent vars
        my <- sum(ytrain) / (n - rmat)
        xtrain <- dataset[-mat[, vim], , drop = FALSE]  ## train set independent vars
        xtest <- dataset[mat[, vim], , drop = FALSE]  ## test set independent vars
        yy <- ytrain - my  ## center the dependent variables
        sa <- svd(xtrain)
        d <- sa$d    ;    v <- t(sa$v)    ;     tu <- t(sa$u)  
        d2 <- d^2    ;    A <- d * tu %*% yy
        
        for ( i in 1:mi ) {
          ## betas <- ( v %*% (tu * ( d / ( d^2 + lambda[i] ) ) ) ) %*% yy 
          betas <- crossprod( v / ( d2 + lambda[i] ), A )
          est <- xtest %*% betas + my 
          pe[i] <- sum( (ytest - est)^2 ) / rmat
        }
        return(pe)
      }
      parallel::stopCluster(cl)
      runtime <- proc.time() - runtime
    }

    per <- as.matrix(per)
    mspe <- as.vector( Rfast::colmeans(per) )
    bias <- per[ , which.min(mspe)] - apply(per, 1, min)  ## TT estimate of bias
    estb <- mean( bias )  ## TT estimate of bias
    names(mspe) <- lambda
    lam <- lambda[ which.min(mspe) ]
    plot(lambda, mspe, xlab = expression(paste(lambda, " values")), ylab = "MSPE", type = "b")
    names(mspe) <- lambda
    performance <- c( min(mspe) + estb, estb)
    names(performance) <- c("Estimated MSPE", "Estimated bias")
  }
  
  list(mspe = mspe, lambda = lam, performance = performance, runtime = runtime)
}
