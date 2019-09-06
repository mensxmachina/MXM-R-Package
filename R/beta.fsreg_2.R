beta.fsreg_2 <- function(target, dataset, iniset = NULL, threshold = 0.05, wei = NULL, tol = 2, ncores = 1) {

  dm <- dim(dataset) 
  if ( is.null(dm) ) {
    n <- length(target)
    p <- 1
  } else {
    n <- dm[1]  ## sample size 
    p <- dm[2]  ## number of variables
  }  
  devi <- dof <- numeric( p )  
  moda <- list()
  k <- 1   ## counter
  tool <- numeric( min(n, p) )
  
  pa <- NCOL(iniset)
  da <- 1:pa
  dataset <- cbind(iniset, dataset)

  threshold <- log(threshold)  ## log of the significance level
  moda <- list()
  k <- 1   ## counter
  con <- log(n)
  sela <- NULL
  #check for NA values in the dataset and replace them with the variable median or the mode
  if( any(is.na( dataset) ) ) {
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
  
  runtime <- proc.time()
  
  devi = dof = numeric(p)
  ini <- beta.reg(target, iniset, wei = wei )
  do <- length(mi$be)
  ini <-   - 2 * mi$loglik

  if (ncores <= 1) {
    for (i in 1:p) {
      mi <- beta.reg( target, dataset[, c(da, pa + i)], wei = wei )
      devi[i] <- 2 * mi$loglik
      dof[i] = length(mi$be) 
    }
    
    stat = devi - ini
    pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
    
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach( i = 1:p, .combine = rbind, .export = "beta.reg" ) %dopar% {
      ww <- beta.reg( target, dataset[, c(da, pa + i)], wei = wei )
      return( c( 2* ww$loglik, length( ww$be ) ) )
    }
    
    stopCluster(cl)
    stat <- mod[, 1] - ini
    pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
  }
  
  mat <- cbind(1:p, pval, stat) 
  colnames(mat) <- c( "variables", "log.p-value", "stat" )
  sel <- which.min(mat[, 2])
  info <- matrix( numeric(3), ncol = 3 )
  info <- matrix( c( 1e300, 0, 0 ), ncol = 3 )
  sela <- pa + sel
  
  if ( mat[sel, 2] < threshold ) {
    info[k, ] <- mat[sel, ]
    mat <- mat[-sel, , drop = FALSE] 
    mi <- beta.reg( target, dataset[, c(da, sel) ], wei = wei )
    tool[k] <-  - 2 * mi$loglik + ( length(mi$be) + 1 ) * con
    moda[[ k ]] <- mi
  }
  
  ############
  ###       k equals 2
  ############   
  if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
    
    k <- 2
    pn <- p - k + 1   
    ini <- 2 * mi$loglik
    do <- length(mi$be)
    if ( ncores <= 1 ) {
      devi <- dof <- numeric(pn)
      for ( i in 1:pn ) {
        ww <- beta.reg( target, dataset[, c(da, sela, pa + mat[i, 1]) ], wei = wei )
        devi[i] <- 2 * ww$loglik
        dof[i] <- length( ww$be )          
      }     
      stat <- devi - ini
      pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
      
    } else {   
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mod <- foreach( i = 1:pn, .combine = rbind, .export = "beta.reg" ) %dopar% {
        ww <- beta.reg( target, dataset[, c(da, sela, pa + mat[i, 1]) ], wei = wei )
        return( c( 2 * ww$loglik, length( ww$be ) ) )
      }     
      stopCluster(cl)  
      stat <- mod[, 1] - ini
      pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )    
    }
    mat[, 2:3] <- cbind(pval, stat)
    ina <- which.min(mat[, 2])
    sel <- pa + mat[ina, 1]    
    
    if ( mat[ina, 2] < threshold ) {
      ma <- beta.reg( target, dataset[, c(da, sela, sel) ], wei = wei )
      tool[2] <-  - 2 * ma$loglik + ( length(ma$be) + 1 ) * con
      if ( tool[ 1 ] - tool[ 2 ] <= tol ) {
        info <- info    
      } else { 
        info <- rbind(info, c( mat[ina, ] ) )
        sela <- c(sela, sel)
        mat <- mat[-ina , , drop = FALSE] 
        moda[[ 2 ]] <- ma
      }
    } else  info <- info
  }  
  ############
  ##       k greater than 2
  ############ 
  if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
    while ( info[k, 2] < threshold &  k < n - 15  &  tool[ k - 1 ] - tool[ k ] > tol  &  nrow(mat) > 0 )  {
      
      ini <-  2 * moda[[ k ]]$loglik 
      do <- length( moda[[ k ]]$be )   
      k <- k + 1   
      pn <- p - k  + 1
      if (ncores <= 1) {  
        devi <- dof <- numeric(pn) 
        for ( i in 1:pn ) {
          ma <- beta.reg( target, dataset[, c(da, sela, pa + mat[i, 1] ) ], wei = wei )
          devi[i] <-  2 * ma$loglik
          dof[i] <- length( ma$be ) 
        }
        stat <- devi - ini 
        pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach( i = 1:pn, .combine = rbind, .export = "beta.reg" ) %dopar% {
          ww <- beta.reg( target, dataset[, c(da, sela, pa + mat[i, 1]) ], wei = wei )
          return( c( 2 * ww$loglik, length( ww$be ) ) )
        }
        stopCluster(cl)
        stat <- mod[, 1] - ini
        pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
      }
      mat[, 2:3] <- cbind(pval, stat)
      ina <- which.min(mat[, 2])
      sel <- pa + mat[ina, 1]    
      
      if ( mat[ina, 2] < threshold ) {
        ma <- beta.reg( target, dataset[, c(da, sela, sel) ], wei = wei )
        tool[k] <-  - 2 * ma$loglik + ( length(ma$be) + 1 ) * con  
        if ( tool[ k - 1 ] - tool[ k  ] <= tol ) {
          info <- rbind(info, c( 1e300, 0, 0 ) )  
        } else { 
          info <- rbind( info, mat[ina, ] )
          sela <- c(sela, sel)
          mat <- mat[-ina , , drop = FALSE]
          moda[[ k ]] <- ma
        }  
      } else   info <- rbind(info, c( 1e300, 0, 0 ) )
    }   
  }   
  runtime <- proc.time() - runtime 
  d <- length(sela)
  final <- NULL
  
  if ( d == 0 ) {
    final <- beta.reg( target, iniset, wei = wei )
    info <- NULL
    
  } else {
    
    if ( d >= 1 ) {
      final <- beta.reg( target, dataset[, c(da, sela) ], wei = wei )
      info <- info[1:d, , drop = FALSE]
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "log.p-value", "stat", "BIC" )
      rownames(info) <- info[, 1]
    }
    
  }    
  
  list(runtime = runtime, mat = t(mat), info = info, ci_test = "testIndBeta", final = final ) 
}         












