mm.fsreg_2 <- function(target, dataset, iniset = NULL, threshold = 0.05, wei = NULL, stopping = "BIC", tol = 2, ncores = 1 ) {
  
  threshold = log(threshold)
  p <- dim(dataset)[2]  ## number of variables
  pval <- stat <- dof <- numeric( p )  
  moda <- list()
  ## percentages
  if ( min( target ) > 0  &  max( target ) < 1 )   target <- log( target / (1 - target) ) 
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any(is.na(dataset)) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
      for ( i in poia )  {
        xi <- dataset[, i]
        if ( is.numeric(xi) )
        {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) )   xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        dataset[, i] <- xi
      }
    }
  }
  
  ## is there already an initial set of variables to start with?
  pa <- NCOL(iniset)
  da <- 1:pa
  dataset <- cbind(iniset, dataset)
  dataset <- as.data.frame(dataset)
  
  n <- length(target)  ## sample size
  tool <- numeric( length( min(n, p) ) )
  info <- matrix( c( 1e300, 0, 0 ), ncol = 3 )
  k <- 1   ## counter k is 1, step 1
  mi <- lm(target ~. , data = as.data.frame(iniset), weights = wei, y = FALSE, model = FALSE)
  runtime <- proc.time()
  
  if (ncores <= 1) {
    for (i in 1:p) {
      ww <- MASS::rlm( target ~., data = dataset[, c(da, pa + i)], maxit = 2000, method = "MM") 
      stat[i] <- 2 * as.numeric( logLik(ww) )
      dof[i] <- length( coef(ww) )
    }
    fit0 <- MASS::rlm( target ~ ., data = data.frame(iniset), maxit = 2000, method = "MM") 
    stat0 <- 2 * as.numeric( logLik(fit0) )
    difa <- stat - stat0
    pval <- pchisq(difa, dof - length( coef(fit0) ), lower.tail = FALSE, log.p = TRUE)
    mat <- cbind(1:p, pval, difa)
    
  } else {
    fit0 = MASS::rlm( target ~., data = data.frame(iniset), maxit = 2000, method = "MM") 
    stat0 = 2 * logLik(fit0)
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach( i = 1:p, .combine = rbind, .export = "rlm", .packages = "MASS" ) %dopar% {
      ww = MASS::rlm( target ~., data = dataset[, c(da, pa + i)], maxit = 2000, method = "MM") 
      return( c(2 * as.numeric( logLik(ww) ), length( coef(ww) ) ) ) 
    }
    stopCluster(cl)
    stat = mod[, 1] - stat0
    pval = pchisq(stat, mod[, 2] - length( coef(fit0) ), lower.tail = FALSE, log.p = TRUE)
    mod = cbind( pval, stat)
    mat <- cbind(1:p, mod)      
  }
  
  colnames(mat) <- c( "variables", "log.p-value", "stat" )
  rownames(mat) <- 1:p
  sel <- which.min(mat[, 2])
  sela <- pa + sel
  
  if ( mat[sel, 2] < threshold ) {
    
    info[k, ] <- mat[sel, ]
    mat <- mat[-sel, , drop = FALSE] 

    if ( stopping == "adjrsq" ) {
      mi = MASS::rlm( target ~., data = dataset[, c(da, sel) ], maxit = 2000, method = "MM")
      r2 = cor( target, fitted(mi) )^2
      tool[k] = 1 - (1 - r2) * (n - 1) / ( n - length( coef(mi) ) - 1 )
      
    } else if ( stopping  == "BIC" ) { 
      mi = MASS::rlm( target ~., data = dataset[, c(da, sel) ], maxit = 2000, method = "MM")
      tool[k] <-  BIC(mi)
    }
    moda[[ k ]] <- mi
    
  }
  ######
  ###### k equal to 2
  ######
  if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
    
    k <- k + 1
    pn <- p - k + 1   
    
    if ( ncores <= 1 ) {
      do = 2 * as.numeric( logLik( mi ) )
      fr = length( coef( mi ) )
      sta = dof = numeric(pn)
      
      for (i in 1:pn) {
        ww = MASS::rlm( target ~., data = dataset[, c(da, sel, pa + mat[i, 1]) ], maxit = 2000, method = "MM")
        sta[i] = 2 * as.numeric( logLik(ww) )
        dof[i] = length( coef(ww) )
      } 
      mat[, 3] = sta - do 
      mat[, 2] = pchisq(mat[, 3], dof - fr, lower.tail = FALSE, log.p = TRUE)      
    } else {
      do <- 2 * as.numeric( logLik( mi ) )
      fr <- length( coef( mi ) )
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mod <- foreach( i = 1:pn, .combine = rbind, .export = c("rlm"), .packages = "MASS" ) %dopar% {
        ww <- MASS::rlm( target ~., data = dataset[, c(da, sel, pa + mat[i, 1]) ], maxit = 2000, method = "MM")
        return( c( 2 * as.numeric( logLik(ww) ), length( coef(ww) ) ) )
      }
      stopCluster(cl)
      stat <- mod[, 1] - do
      pval <- pchisq(stat, mod[, 2] - fr, lower.tail = FALSE, log.p = TRUE)
      mod <- cbind( pval, stat)
      mat <- cbind(mat[, 1], mod)   
      
    }
  }
  
  ina <- which.min(mat[, 2])
  sel <- pa + mat[ina, 1]     
  
  if ( stopping == "adjrsq" ) {
  
    if ( mat[ina, 2] < threshold ) {
      ma <- MASS::rlm( target ~., data = dataset[, c(da, sela, sel) ], maxit = 2000, method = "MM")
      r2 <- cor( target, fitted(ma) )^2
      tool[k] <- 1 - (1 - r2) * (n - 1) / ( n - length( coef(ma) ) - 1)
      
      if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
        info <- info  
      } else {  
        info <- rbind(info, mat[ina, ] )
        sela <- c(sela, sel)
        mat <- mat[-ina, , drop = FALSE] 
        moda[[ k ]] <- ma
      }
      
    } else  info <- info
    
  } else if ( stopping == "BIC" ) {
  
    if ( mat[ina, 2] < threshold ) {
      ma <- MASS::rlm(target ~ target ~., data = dataset[, c(da, sela, sel) ], maxit = 2000, method = "MM")
      tool[k] <- BIC(ma)
      if ( tool[ k - 1] - tool[ k ] <= tol ) {
        info <- info
      } else {  
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina, , drop = FALSE]
        moda[[ k ]] <- ma
      }
    } else  info <- info
    
  }
  ###########
  ######   k greater than 2
  ###########
  if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
    
    while ( info[k, 2] < threshold &  k < n - 15  &  abs( tool[ k ] - tool[ k - 1 ] ) > tol  &  nrow(mat) > 0 )  {
      
      k <- k + 1   
      pn <- p - k + 1 
      
      if ( ncores <= 1 ) {
        do = 2 * as.numeric( logLik( moda[[ k - 1 ]] ) )
        fr = length( coef( moda[[ k - 1 ]] ) )
        sta = dof = numeric(pn)
        for (i in 1:pn) {
          ww = MASS::rlm( target ~., data = dataset[, c(da, sela, pa + mat[i, 1]) ], maxit = 2000, method = "MM")
          sta[i] = 2 * as.numeric( logLik(ww) )
          dof[i] = length( coef(ww) )
        } 
        mat[, 3] = sta - do
        mat[, 2] = pchisq(mat[, 3], dof - fr, lower.tail = FALSE, log.p = TRUE)

      } else {
        do = 2 * as.numeric( logLik( ma) )
        fr = length( coef( ma ) )
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach( i = 1:pn, .combine = rbind, .export = "rlm", .packages = "MASS" ) %dopar% {
          ww <- MASS::rlm( target ~., data = dataset[, c(da, sela, pa + mat[i, 1]) ], maxit = 2000, method = "MM")
          return( c( 2 * as.numeric( logLik(ww) ), length( coef(ww) ) ) )
        }
        stopCluster(cl)  
        stat = mod[, 1] - do
        pval = pchisq(stat, mod[, 2] - fr, lower.tail = FALSE, log.p = TRUE)
        mod = cbind( pval, stat)
        mat <- cbind( mat[, 1], mod )   
      }
      
      ina <- which.min(mat[, 2])
      sel <- pa + mat[ina, 1]   
      
      if ( stopping == "BIC" ) {
	  
        if ( mat[ina, 2] < threshold ) {
           ma <- MASS::rlm( target ~., data = dataset[, c(da, sela, sel)], maxit = 2000, method = "MM")
           tool[k] <-  BIC(ma)

          if ( tool[ k - 1] - tool[ k ] <= tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )
            
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- c(sela, sel)
            mat <- mat[-ina, , drop = FALSE] 
            moda[[ k ]] <- ma
          }
          
        } else  info <- rbind(info, c( 1e300, 0, 0 ) )
        
      } else if ( stopping == "adjrsq" ) {
	    if ( mat[ina, 2] < threshold ) {
          ma <- MASS::rlm( target ~., data = dataset[, c(da, sela, sel)], maxit = 2000, method = "MM")
          r2 <- cor(target, fitted(ma) )^2
          tool[k] = 1 - (1 - r2) * (n - 1) / ( n - length( coef(ma) ) - 1)
             
          if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )          
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- info[, 1]
            mat <- mat[-ina, , drop = FALSE]
            moda[[ k ]] <- ma
          }
          
        } else   info <- rbind(info, c( 1e300, 0, 0 ) )
        
      }
    }
    
  }
  
  runtime <- proc.time() - runtime
  
  d <- length(moda)
  
  if (d == 0) {
    final <- MASS::rlm( target ~., data.frame( iniset ), maxit = 2000, method = "MM")    
  } else {
    final <- MASS::rlm( target ~., dataset[, c(da, sela) ], maxit = 2000, method = "MM")
    info <- info[1:d, , drop = FALSE ]
    info <- cbind( info, tool[ 1:d ] ) 
    colnames(info) <- c( "variables", "log.p-value", "stat", stopping )
    rownames(info) <- info[, 1]
  }
  
  list( runtime = runtime, mat = t(mat), info = info, ci_test = "testIndMMReg", final = final ) 
}
