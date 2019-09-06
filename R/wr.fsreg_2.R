wr.fsreg_2 <- function(target, dataset, iniset = NULL, wei = NULL, threshold = 0.05, tol = 2, ncores = 1) {
  
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
  threshold <- log(threshold)
  
  if ( is.null(iniset) ) {
    da <- 0
    pa <- 0
  } else {
    pa <- NCOL(iniset)
    da <- 1:pa
    dataset <- cbind(iniset, dataset)
  }  
  dataset <- as.data.frame(dataset)
  
  runtime <- proc.time()

  devi = dof = numeric(p)
  if ( pa == 0 ) {
    mi <- survival::survreg( target ~ 1, weights = wei )
    ini <- 2 * logLik(mi)
  } else  {
    mi<- survival::survreg(target ~., data = as.data.frame( iniset ), weights = wei )
    ini <- 2* logLik(mi)
  }  
  do <- length(mi$coefficients)
  if (ncores <= 1) {
    for (i in 1:p) {
      mi <- survival::survreg( target ~ . , as.data.frame( dataset[, c(da, pa + i)] ), weights = wei )
      devi[i] <- 2 * logLik(mi)
      dof[i] <- length( mi$coefficients ) 
    }
    
    stat = devi - ini
    pval = pchisq( stat, dof - 1, lower.tail = FALSE, log.p = TRUE )
    
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach( i = 1:p, .combine = rbind, .export = "survival", .packages = "survreg") %dopar% {
      ww <- survival::survreg( target ~., data = as.data.frame( dataset[, c(da, pa + i)] ), weights = wei )
      return( c( 2 * logLik(ww), length( ww$coefficients )  ) )
    }
    
    stopCluster(cl)
    stat <- mod[, 1] - ini
    pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
  }
  
  mat <- cbind(1:p, pval, stat) 
  
  colnames(mat)[1] <- "variables"
  rownames(mat) <- 1:p
  
  sel <- which.min(pval)
  info <- matrix( c( 1e300, 0, 0 ), ncol = 3 )
  sela <- sel
  
  if ( mat[sel, 2] < threshold ) {
    info[k, ] <- mat[sel, ]
    mat <- mat[-sel, , drop = FALSE] 
    mi <- survival::survreg( target ~., data = as.data.frame( dataset[, c(da, sel) ] ), weights = wei )
    la <- logLik(mi)
    tool[k] <-  - 2 * as.numeric( la ) + attr(la, "df") * con
    moda[[ k ]] <- mi
  }
  ############
  ###       k equals 2
  ############ 
  if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
    
    k <- k + 1
    pn <- p - k + 1   
    ini <- moda[[ 1 ]]$deviance  ## residual deviance
    do <- length( coef( moda[[ 1 ]]  ) ) 
    devi <- dof <- numeric( pn )  
    
    if ( ncores <= 1 ) {
      for ( i in 1:pn ) {
        ww <- survival::survreg( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), weights = wei )
        devi[i] <- 2 * logLik(ww)
        dof[i] <- length( ww$coefficients )        
      }
      
      stat <- devi - ini
      pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
      
    } else {
      
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mod <- foreach( i = 1:pn, .combine = rbind, .export = "survival", .packages = "survreg") %dopar% {
        ww <- survival::survreg( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), weights = wei )
        return( c( 2 * logLik(ww), length( ww$coefficients ) ) )
      }
      stopCluster(cl)
      stat <- mod[, 1] - ini
      pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
    }
    
    mat[, 2:3] <- cbind(pval, stat)
    ina <- which.min(mat[, 2])
    sel <- mat[ina, 1]    
    
    if ( mat[ina, 2] < threshold ) {
      ma <- survival::survreg( target ~., data=as.data.frame( dataset[, c(da, sela, sel) ] ), weights = wei )
      la <- logLik(ma)
      tool[k] <-  - 2 * as.numeric( la ) + attr(la, "df") * con
      
      if ( tool[ k - 1 ] - tool[ k ] <= tol ) {
        info <- info
        
      } else { 
        info <- rbind(info, c( mat[ina, ] ) )
        sela <- info[, 1]
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
      
      ini <- moda[[ k ]]$deviance  ## residual deviance
      do <- length( coef( moda[[ k ]]  ) ) 
      k <- k + 1   
      pn <- p - k  + 1
      devi <- dof <- numeric( pn )  
      
      if (ncores <= 1) {  
        for ( i in 1:pn ) {
          ma <- survival::survreg( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1] ) ] ), weights = wei )
          devi[i] <- 2 * logLik(ma)
          dof[i] <- length( ma$coefficients ) 
        }
        stat <- devi - ini
        pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach( i = 1:pn, .combine = rbind, .export = "survival", .packages = "survreg") %dopar% {
          ww <- survival::survreg( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), weights = wei )
          return( c( 2 * logLik(ww), length( ww$coefficients ) ) )
        }
        stopCluster(cl)
        stat <- mod[, 1] - ini
        pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE ) 
      }
      
      mat[, 2:3] <- cbind(pval, stat)
      
      ina <- which.min(mat[, 2])
      sel <- mat[ina, 1]    
      
      if ( mat[ina, 2] < threshold ) {
        ma <- survival::survreg( target ~., data = as.data.frame( dataset[, c(da, sela, sel) ] ), weights = wei )
        la <- logLik(ma)
        tool[k] <-  - 2 * as.numeric( la ) + attr(la, "df") * con
        if ( tool[ k - 1 ] - tool[ k  ] < tol ) {
          info <- rbind(info, c( 1e300, 0, 0 ) )
          
        } else { 
          info <- rbind( info, mat[ina, ] )
          sela <- info[, 1]
          mat <- mat[-ina, , drop = FALSE]
          moda[[ k ]] <- ma
        }  ## end if ( tool[ k - 1 ] - tool[ k  ] < tol ) 
        
      } else   info <- rbind(info, c( 1e300, 0, 0 ) )
      
    }
    
  } 
  
  runtime <- proc.time() - runtime
  
  d <- length(moda)
  
  if ( d == 0 ) {
    final <- survival::survreg( target ~., data = as.data.frame( iniset ), weights = wei )
    info <- NULL
    
  } else {
      final <- survival::survreg( target ~., data = as.data.frame( dataset[, c(da, sela) ] ), weights = wei )
      info <- info[1:d, , drop = FALSE]
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "log.p-value", "stat", "BIC" )
      rownames(info) <- info[, 1]
  }    
  
  list(runtime = runtime, mat = t(mat), info = info, ci_test = "testIndTobit", final = final ) 
}    










