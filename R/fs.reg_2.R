fs.reg_2 <- function(target, dataset, iniset = NULL, threshold = 0.05, wei = NULL, test = NULL, user_test = NULL, stopping = "BIC", tol = 2, ncores = 1) {
  
  threshold <- log(threshold)
  
  dm <- dim(dataset) 
  if ( is.null(dm) ) {
    n <- length(dataset)
    p <- 1
  } else {
    n <- dm[1]  ## sample size 
    p <- dm[2]  ## number of variables
  }  
  devi <- dof <- numeric( p )  
  moda <- list()
  k <- 1   ## counter
  tool <- numeric( min(n, p) )
  con <- log(n)

  pa <- NCOL(iniset)
  da <- 1:pa
  dataset <- cbind(iniset, dataset)
  dataset <- as.data.frame(dataset)
  ci_test <- test

  if ( ( test == "testIndLogistic" )  ||  test == "testIndBinom"  || test == "testIndPois" ) {
    
    result <- glm.fsreg_2(target, dataset, iniset = iniset, wei = wei, threshold = exp(threshold), tol = tol, 
                          ncores = 1) 
      
  } else if ( test == "testIndReg" || test == "testIndFisher" ) {
      result <- lm.fsreg_2( target, dataset, iniset = iniset, wei = wei, threshold = exp(threshold), stopping = stopping, tol = tol, 
                            ncores = ncores ) 
    
  } else if ( test == "testIndBeta" ) {
    result <- beta.fsreg_2(target, dataset, iniset = iniset, threshold = exp(threshold), wei = wei, tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndZIP" ) {
    result <- zip.fsreg_2(target, dataset, iniset = iniset, threshold = exp(threshold), wei = wei, tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndGamma" ) {
    result <- gammafsreg_2(target, dataset, iniset = iniset, threshold = exp(threshold), wei = wei, tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndNormLog" ) {
    result <- normlog.fsreg_2(target, dataset, iniset = iniset, threshold = exp(threshold), wei = wei, tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndTobit" ) {
    result <- tobit.fsreg_2(target, dataset, iniset = iniset, threshold = exp(threshold), wei = wei, tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndClogit" ) {
    result <- clogit.fsreg_2(target, dataset, iniset = iniset, threshold = exp(threshold), tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndQBinom" ) {
    result <- quasibinom.fsreg_2(target, dataset, iniset = iniset, threshold = exp(threshold), tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndQPois" ) {
    result <- quasipois.fsreg_2(target, dataset, iniset = iniset, threshold = exp(threshold), tol = tol, ncores = ncores) 
    
  } else {

    if ( test == "censIndCR" ) {
      test <- survival::coxph 
      stopping <- "BIC"
      
    } else if ( test == "censIndWR" ) {
      test <- survival::survreg
      stopping <- "BIC"
      
    } else if ( test == "testIndOrdinal" ) {
      test <- ordinal::clm
      stopping <- "BIC"
    } else if (test == "testIndMultinom") {
      test <- nnet::multinom
      stopping <- "BIC"
      
    } else if ( test == "testIndNB" ) {
      test <- MASS::glm.nb
      stopping <- "BIC"
      
    } else if ( test == "testIndRQ" ) {
      test <- quantreg::rq
      stopping <- "BIC"
    } 
    
    if ( !is.null(user_test) )  {
      test <- user_test 
      ci_test <- "user_test"
    }  
  
  dataset <- as.data.frame(dataset)
  runtime <- proc.time()
  
  devi = dof = numeric(p)
  ini <- test(target ~., data = as.data.frame( iniset ), weights = wei)  ## residual deviance
  do <- length( coef(ini) )
  if (ci_test == "censIndCR") {
    ini <- 2 * ini$loglik 
  }  else  ini <-  2 * as.numeric( logLik(ini) )  ## initial 

  
  if (ncores <= 1) {
    for (i in 1:p) {
      mi <- test( target ~ . , data = dataset[, c(da, pa + i)], weights = wei )
      devi[i] <- 2 * logLik(mi)
      dof[i] = length( coef( mi ) ) 
    }
    stat <- devi - ini
    pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
    
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach( i = 1:p, .combine = rbind, .packages = c("MASS", "quantreg", "nnet", "survival", "ordinal") ) %dopar% {
      ww <- test( target ~., data = dataset[, c(da, pa + i)], weights = wei )
      return( c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) ) ) )
    }
    
    stopCluster(cl)
    stat <- abs(mod[, 1] - ini)
    pval <- pchisq( stat, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE )
  }
  
  mat <- cbind(1:p, pval, stat) 
  colnames(mat)[1] <- "variables"
  rownames(mat) <- 1:p
  sel <- which.min(pval)
  info <- matrix( c( 1e300, 0, 0 ), ncol = 3 )
  sela <- pa + sel
  
  if ( mat[sel, 2] < threshold ) {
    info[k, ] <- mat[sel, ]
    mat <- mat[-sel, , drop = FALSE] 
    mi <- test( target ~., data = dataset[, c(da, sela) ], weights = wei )
    la <- logLik(mi)
    tool[1] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
    moda[[ 1 ]] <- mi
  }
  
  ############
  ###       k equals 2
  ############ 
  
  if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
    
    k <- k + 1
    pn <- p - k + 1   
    ini <- 2 * logLik(moda[[ 1 ]])
    do <- length( coef( moda[[ 1 ]]  ) ) 
    devi <- dof <- numeric( pn )  
    
    if ( ncores <= 1 ) {
      devi <- dof <- numeric(pn)
      for ( i in 1:pn ) {
        ww <- test( target ~., data = dataset[, c(da, sela, pa + mat[i, 1]) ], weights = wei )
        devi[i] <- 2 * as.numeric( logLik(ww) )
        dof[i] <- length( coef( ww ) )          
      }
      stat <- abs( devi - ini )
      pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
      
    } else {
      
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mod <- foreach( i = 1:pn, .combine = rbind, .packages = c("MASS", "quantreg", "nnet", "survival", "ordinal") ) %dopar% {
        ww <- test( target ~., data = dataset[, c(da, sela, pa + mat[i, 1]) ], weights = wei )
        return( c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) ) ) )
      }
      
      stopCluster(cl)
      
      stat <- abs(mod[, 1] - ini)
      pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
      
    }
    
    mat[, 2:3] <- cbind(pval, stat)
    ina <- which.min(mat[, 2])
    sel <- pa + mat[ina, 1]    
    
    if ( mat[ina, 2] < threshold ) {
      ma <- test( target ~ dataset[, sela] + dataset[, sel], weights = wei )
      la <- logLik(ma)
      tool[2] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
      
      if ( tool[ k - 1 ] - tool[ k ] <= tol ) {
        info <- info
        
      } else { 
        info <- rbind(info, c( mat[ina, ] ) )
        sela <- c(sela, sel)
        mat <- mat[-ina , , drop = FALSE] 
        moda[[ k ]] <- ma
      }
      
    } else   info <- info
    
  }
  
  ############
  ###       k greater than 2
  ############ 
  
  if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
    while ( info[k, 2] < threshold &  k < n - 15  &  tool[ k - 1 ] - tool[ k ] > tol  &  nrow(mat) > 0  )  {
      
      ini =  2 * as.numeric( logLik( moda[[ k ]] ) ) 
      do = length( coef( moda[[ k ]]  ) )      
      k <- k + 1   
      pn <- p - k  + 1
      devi <- dof <- numeric( pn )  
      
      if (ncores <= 1) {  
        for ( i in 1:pn ) {
          ma <- test( target ~., data = dataset[, c(da, sela, pa + mat[i, 1] ) ], weights = wei )
          devi[i] <-  2 * as.numeric( logLik(ma) ) 
          dof[i] <- length( coef( ma ) ) 
        }
        
        stat <- devi - ini
        pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        devi <- dof <- numeric(pn)
        mod <- foreach( i = 1:pn, .combine = rbind, .packages = c("MASS", "quantreg", "nnet", "survival", "ordinal") ) %dopar% {
          ww <- test( target ~., data = dataset[, c(da, sela, pa + mat[i, 1]) ], weights = wei )
          return( c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) ) ) )
        }
        
        stopCluster(cl)
  
        stat <- mod[, 1] - ini
        pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
      }
      
      mat[, 2:3] <- cbind(pval, stat)
      ina <- which.min(mat[, 2])
      sel <- pa + mat[ina, 1]    
      
      if ( mat[ina, 2] < threshold ) {
        ma <- test( target ~., data = dataset[, c(da, sela, sel) ], weights = wei )
        la <- logLik(ma)
        tool[k] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
        if ( tool[ k - 1 ] - tool[ k  ] < tol ) {
          info <- rbind(info, c( 1e300, 0, 0 ) )
          
        } else { 
          info <- rbind( info, mat[ina, ] )
          sela <- c(sela, sel)
          mat <- mat[-ina, , drop = FALSE]
          moda[[ k ]] <- ma
        } 
        
      } else   info <- rbind(info, c( 1e300, 0, 0 ) )
      
    }
    
  } 
  
  runtime <- proc.time() - runtime
  
  d <- length(moda)
  final <- NULL
  
  if ( d == 0 ) {
    final <- test( target ~., data = as.data.frame( iniset ), weights = wei )
    info <- NULL
    
  } else {
    final <- test( target ~., data = dataset[, c(da, sela) ], weights = wei )
    info <- info[1:d, , drop = FALSE]
    info <- cbind( info, tool[ 1:d ] ) 
    colnames(info) <- c( "variables", "log.p-value", "stat", "BIC" )
    rownames(info) <- info[, 1]
  }    
  
  result <- list( runtime = runtime, mat = t(mat), info = info, ci_test = ci_test, final = final ) 
  
  }  
  
}    










