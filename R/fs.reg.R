fs.reg <- function(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, test = NULL, user_test = NULL, stopping = "BIC", tol = 2, ncores = 1) {
  ## target can be Real valued (normal), binary (binomial) or counts (poisson)
  ## dataset is a matrix or a data.frame with the predictor variables
  ## test is the test used, but not really required because depending on the data this will be decided.
  ## there is no hrm is psecifying this though
  ## threshold is the level of significance
  ## method can be either BIC or adjrsq (for linear models only). The BIC is a consistent method for selecting
  ## models, thus we also use it to avoid overfitting
  ## stopping is based on "BIC"
  ## tol is the tolerance value for the method. If BIC is used as the stopping rule, the default is 2, but usually can be 2 or 4.
  ## If BIC is used as a way to proceed, the tol is 0.
  ## ncores is for parallel processing 
  
  if ( !is.null(ini) ) {
    result <- fs.reg_2(target, dataset, iniset = ini, test = test, wei = wei, threshold = threshold, tol = tol, ncores = ncores) 
    
  } else { 
  
  threshold <- log(threshold)  ## log of the significance level
  dm <- dim(dataset)
  p <- dm[2]  ## number of variables
  devi <- dof <- numeric( p )  
  moda <- list()
  k <- 1   ## counter
  n <- dm[1]  ## sample size
  tool <- numeric( min(n, p) )
  result <- NULL
  con <- log(n)
  sela <- NULL
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any( is.na(dataset) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
      for( i in poia )  {
        xi <- dataset[, i]
        if ( is.numeric(xi) ) {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) )     xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        dataset[, i] <- xi
      }
    }
  }
  ##################################
  # target checking and initialize #
  ##################################
  if ( is.null(test)  &  is.null(user_test) ) {
    
    ## surival data
    if ( sum( class(target) == "Surv" ) == 1 ) {
      ci_test <- test <- "censIndCR"
      ## ordinal, multinomial or perhaps binary data
    } else if ( length( unique(target) ) == 2 ) {
      ci_test <- test <- "testIndLogistic"
      ## count data
    } else if ( length( unique(target) ) > 2  &  !is.factor(target) ) {
      if ( sum( round(target) - target ) == 0 ) {
         ci_test <- test <- "testIndPois"
      } else  ci_test <- test <- "testIndReg"  
    }
  }

  ci_test <- test
  #cat(test)
  if ( test == "testIndLogistic" | test == "testIndBinom"  | test == "testIndPois" ) {
    result <- glm.fsreg( target, dataset, wei = wei, threshold = exp(threshold), tol = tol, ncores = ncores) 

  } else if ( test == "testIndFisher"  &  is.matrix(dataset)  &  is.null(wei) ) {
    if (stopping == "adjrsq")   stopping = "ar2"
    result <- Rfast::cor.fsreg(target, dataset, threshold = exp(threshold), tolb = tol, tolr = tol, stopping = stopping)
    
  } else if ( test == "testIndReg" | ( test == "testIndFisher"  &  !is.matrix(dataset) ) ) {
    result <- lm.fsreg( target, dataset, wei = wei, threshold = exp(threshold), stopping = stopping, tol = tol, ncores = ncores ) 
  
  } else if ( test == "testIndMMReg" ) {
    result <- mm.fsreg( target, dataset, wei = wei, threshold = exp(threshold), tol = tol, ncores = ncores) 
	
  } else if ( test == "testIndBeta" ) {
    result <- beta.fsreg(target, dataset, threshold = exp(threshold), wei = wei, tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndZIP" ) {
    result <- zip.fsreg(target, dataset, threshold = exp(threshold), wei = wei, tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndGamma" ) {
    result <- gammafsreg(target, dataset, threshold = exp(threshold), wei = wei, tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndNormLog" ) {
    result <- normlog.fsreg(target, dataset, threshold = exp(threshold), wei = wei, tol = tol, ncores = ncores) 

  } else if ( test == "testIndTobit" ) {
    result <- tobit.fsreg(target, dataset, threshold = exp(threshold), wei = wei, tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndClogit" ) {
    result <- clogit.fsreg(target, dataset, threshold = exp(threshold), tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndQBinom" ) {
    result <- quasibinom.fsreg(target, dataset, threshold = exp(threshold), tol = tol, ncores = ncores) 
    
  } else if ( test == "testIndQPois" ) {
    result <- quasipois.fsreg(target, dataset, threshold = exp(threshold), tol = tol, ncores = ncores) 
      
  } else {
    
    dataset <- as.data.frame(dataset)
    
    if ( test == "censIndCR" ) {
      test <- survival::coxph 
      stopping <- "BIC"
      
    } else if ( test == "censIndWR" ) {
      test <- survival::survreg
      stopping <- "BIC"
      
    } else if ( test == "testIndOrdinal" ) {
	  test <- ordinal::clm
      stopping <- "BIC"
      
	} else if ( test == "testIndMultinom" )  {
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
    
    runtime <- proc.time()
  
    devi <- dof <- numeric(p)
    ini <- test( target ~ 1, weights = wei ) 
    dof1 <- length( coef(ini) )
    if (ci_test == "censIndCR") {
      ini <- 2 * ini$loglik 
    }  else  ini <-  2 * as.numeric( logLik(ini) )  ## initial 
    
    if (ncores <= 1) {
      for (i in 1:p) {
        mi <- test( target ~ dataset[, i], weights = wei )
        devi[i] <-  2 * as.numeric( logLik(mi) )
        dof[i] <- length( coef( mi ) ) 
      }
      
      stat <- abs( devi - ini )
      if ( ci_test == "censIndCR" )  {
        dof <- dof + 1 
        dof1 <- 0
      }  
      pval <- pchisq( stat, dof - dof1, lower.tail = FALSE, log.p = TRUE )
      
    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mod <- foreach( i = 1:p, .combine = rbind, .packages = c("MASS", "quantreg", "nnet", "survival", "ordinal") ) %dopar% {
        ww <- test( target ~ dataset[, i], weights = wei )
        return( c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) ) ) )
      }
      stopCluster(cl)
      
      stat <-  abs( mod[, 1] - ini )
      if ( ci_test == "censIndCR" )  {
        mod[, 2] <- mod[, 2] + 1
        dof1 <- 0
      }  
      pval <- pchisq( stat, mod[, 2] - dof1, lower.tail = FALSE, log.p = TRUE )
    }
    
    mat <- cbind(1:p, pval, stat) 
    colnames(mat) <- c( "variables", "log.p-value", "stat" )
    rownames(mat) <- 1:p
    sel <- which.min(pval)
    info <- matrix( numeric(3), ncol = 3 )

    if ( mat[sel, 2] < threshold ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE] 
      sela <- sel
      mi <- test( target ~ dataset[, sel], weights = wei )
      la <- logLik(mi)
      tool[1] <-  - 2 * as.numeric( la ) + attr(la, "df") * con
      moda[[ 1 ]] <- mi
    }  else  {
       info <- info  
       sela <- NULL
    }
    ############
    ###       k equals 2
    ############ 
    if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
      
      k <- 2
      pn <- p - k + 1   
      ini <-  2 * as.numeric( la ) 
      do <- length( coef( moda[[ 1 ]]  ) ) 
      devi <- dof <- numeric( pn )  
      
      if ( ncores <= 1 ) {
        devi <- dof <- numeric(pn)
        for ( i in 1:pn ) {
          ww <- test( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei )
          devi[i] <-  2 * as.numeric( logLik(ww) )
          dof[i] <- length( coef( ww ) )          
        }
        stat <- abs( devi - ini )
        pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach( i = 1:pn, .combine = rbind, .packages = c("MASS", "quantreg", "nnet", "survival", "ordinal") ) %dopar% {
          ww <- test( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei )
          return( c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) ) ) )
        }
        stopCluster(cl)
        stat <- abs( mod[, 1] - ini )
        pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
      }
    
    mat[, 2:3] <- cbind(pval, stat)
    ina <- which.min(mat[, 2])
    sel <- mat[ina, 1]    
    
    if ( mat[ina, 2] < threshold ) {
      ma <- test( target ~ dataset[, sela] + dataset[, sel], weights = wei )
      la <- logLik(ma)
      tool[2] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
      
      if ( tool[ 1 ] - tool[ 2 ] <= tol ) {
        info <- info
        
      } else { 
        info <- rbind(info, c( mat[ina, ] ) )
        sela <- info[, 1]
        mat <- mat[-ina, , drop = FALSE ] 
        moda[[ 2 ]] <- ma
      }
    } else  info <- info  
  }
  ############
  ###       k greater than 2
  ############ 
  if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
    while ( info[k, 2] < threshold   & k < n - 15  & tool[ k - 1 ] - tool[ k ] > tol  &  nrow(mat) > 0 )  {
      
      ini =  2 * as.numeric( logLik( moda[[ k ]] ) ) 
      do = length( coef( moda[[ k ]]  ) ) 
      k <- k + 1   
      pn <- p - k  + 1
      devi <- dof <- numeric( pn )  
      
      if (ncores <= 1) {  
        devi = dof = numeric(pn) 
        for ( i in 1:pn ) {
          ma <- test( target ~., data = dataset[, c(sela, mat[i, 1] ) ], weights = wei )
          devi[i] <-  2 * as.numeric( logLik(ma) ) 
          dof[i] = length( coef( ma ) ) 
        }
        
        stat = abs( devi - ini )
        pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mod <- foreach( i = 1:pn, .combine = rbind, .packages = c("MASS", "quantreg", "nnet", "survival", "ordinal") ) %dopar% {
            ww <- test( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei )
            return( c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) ) ) )
          }
          stopCluster(cl)
          stat = abs( mod[, 1] - ini )
          pval = pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
        }
        
        mat[, 2:3] <- cbind(pval, stat)
        ina <- which.min(mat[, 2])
        sel <- mat[ina, 1]    
        
        if ( mat[ina, 2] < threshold ) {
            ma <- test( target ~., data = dataset[, c(sela, sel) ], weights = wei )
            la <- logLik(ma)
            tool[k] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
 
          if ( tool[ k - 1 ] - tool[ k  ] <= tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )
            
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- info[, 1]
            mat <- mat[-ina , , drop = FALSE]
            moda[[ k ]] <- ma
          } 
          
        } else   info <- rbind(info, c( 1e300, 0, 0 ) )
    } 
  }  
    
    runtime <- proc.time() - runtime
    
    d <- length(sela)
    final <- NULL
    if ( d >= 1 ) {
      final <- test( target ~., data = dataset[, sela, drop = FALSE], weights = wei )
      info <- info[1:d,  , drop = FALSE]
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "log.p-values", "stat", "BIC" )
      rownames(info) <- info[, 1]
    }
    
    result = list(runtime = runtime, mat = mat, info = info, ci_test = ci_test, final = final ) 
  }  ## end else if not other forward regressions       
  
  }
  
  result
  
}    










