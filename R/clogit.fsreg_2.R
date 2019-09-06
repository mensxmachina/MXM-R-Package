clogit.fsreg_2 <- function(target, dataset, iniset = NULL, wei = NULL, threshold = 0.05, tol = 2, ncores = 1) {

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
  id <- target[, 2] #the patient id
  case <- as.logical(target[, 1]);  #
  pa <- NCOL(iniset)
  da <- 1:pa
  dataset <- cbind(iniset, dataset)
  dataset <- as.data.frame(dataset) 
    
  runtime <- proc.time()
  
  devi = dof = numeric(p)
  mi <- survival::clogit( case ~ . + strata(id), data = dataset[, 1:da, drop = FALSE] )
  ini <- 2 * logLik(mi)
  tool[k] <-  BIC(mi)
  moda[[ k ]] <- mi
  oop <- options(warn = -1) 
  on.exit( options(oop) )
	
  if (ncores <= 1) {
    for (i in 1:p) {
      mi <- survival::clogit( case ~ dataset[, pa + i] + strata(id) )
      devi[i] <- 2 * mi$loglik
      dof[i] <- length( mi$coefficients ) 
    }
    stat = devi - ini
    pval = pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE )
    
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach( i = 1:p, .combine = rbind, .export = "survival", .packages = "clogit") %dopar% {
      ww <- survival::clogit( case ~ dataset[, pa + i] + strata(id) )
      return( c( 2 * ww$loglik, length( ww$coefficients ) ) )
    }
    stopCluster(cl)
    stat <- mod[, 1] - ini
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
    mi <- survival::clogit( case ~. + strata(id), data = dataset[, c(da, sela), drop = FALSE ] )
    tool[k] <-  BIC(mi)
    moda[[ k ]] <- mi
  }
  ############
  ###       k equals 2
  ############ 
  if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
    
    k <- k + 1
    pn <- p - k + 1   
    ini <- 2 * logLik(moda[[ 1 ]])  ## residual deviance
    do <- length( coef( moda[[ 1 ]]  ) ) 
    
    if ( ncores <= 1 ) {
      devi <- dof <- numeric(pn)
      for ( i in 1:pn ) {
        ww <- survival::clogit( case ~. + strata(id), data = dataset[, c(da, sela, pa + mat[i, 1]), drop = FALSE ] )
        devi[i] <- 2 * logLik(ww)
        dof[i] <- length( ww$coefficients )        
      }
      
      stat <- devi - ini
      pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
      
    } else {
      
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mod <- foreach( i = 1:pn, .combine = rbind, .export = "survival", .packages = "clogit") %dopar% {
        ww <- survival::clogit( case ~. + strata(id), data = dataset[, c(da, sela, pa + mat[i, 1]), drop = FALSE ] )
        return( c( 2 * logLik(ww), length( ww$coefficients ) ) )
      }
      stopCluster(cl)
      
      stat <- mod[, 1] - ini
      pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
    }
    
    mat[, 2:3] <- cbind(pval, stat)
    ina <- which.min(mat[, 2])
    sel <- pa + mat[ina, 1]    
    
    if ( mat[ina, 2] < threshold ) {
      ma <- survival::clogit( case ~. + strata(id), data = dataset[, c(da, sela, sel), drop = FALSE ] )
      tool[k] <- BIC(ma)
      
      if ( tool[ k - 1 ] - tool[ k ] <= tol ) {
        info <- info
        
      } else { 
        info <- rbind(info, c( mat[ina, ] ) )
        sela <- c(sela, sel)
        mat <- mat[-ina , , drop = FALSE] 
        moda[[ k ]] <- ma
      }
      
    } else   info <- info

  ############
  ###       k greater than 2
  ############ 
  
  if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
    while ( info[k, 2] < threshold &  k < n - 15  &  tool[ k - 1 ] - tool[ k ] > tol  &  nrow(mat) > 0  )  {
      
      ini <- 2 * logLik(moda[[ k ]]) ## residual deviance
      do <- length( coef( moda[[ k ]]  ) ) 
      k <- k + 1   
      pn <- p - k  + 1
      
      if (ncores <= 1) {  
        devi = dof = numeric(pn) 
        for ( i in 1:pn ) {
          ma <- survival::clogit( case ~. + strata(id), data = dataset[, c(da, sela, pa + mat[i, 1] ), drop = FALSE ] )
          devi[i] <- 2 * logLik(ma)
          dof[i] <- length( ma$coefficients ) 
        }
        
        stat <- devi - ini
        pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach( i = 1:pn, .combine = rbind, .export = "survival", .packages = "clogit") %dopar% {
          ww <- survival::clogit( case ~. + strata(id), data = dataset[, c(da, sela, pa + mat[i, 1]), drop = FALSE ] )
          return( c( 2 * logLik(ww), length( ww$coefficients ) ) )
        }
        stopCluster(cl)
        
        stat <- mod[, 1] - ini
        pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE ) 
      }
      
      mat[, 2:3] <- cbind(pval, stat)
      ina <- which.min(mat[, 2])
      sel <- pa + mat[ina, 1]    
      
      if ( mat[ina, 2] < threshold ) {
        ma <- survival::clogit( case ~. + strata(id), data = dataset[, c(da, sela, sel), drop = FALSE ] )
        tool[k] <-  BIC(ma)
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
  
  }  ## end if ( info[k, 2] < threshold  &  nrow(mat) > 0 )
  
  runtime <- proc.time() - runtime
  d <- length(moda)
  
  if ( d == 1 ) {
    final <- survival::clogit( case ~. + strata(id), data = as.data.frame( iniset ) )
    info <- NULL
    
  } else {
    final <- survival::clogit( case ~. + strata(id), data = dataset[, c(da, sela), drop = FALSE]  )
    info <- info[1:d, , drop = FALSE]
    info <- cbind( info, tool[ 1:d ] ) 
    colnames(info) <- c( "variables", "log.p-value", "stat", "BIC" )
    rownames(info) <- info[, 1]
  }    
  
  list(runtime = runtime, mat = t(mat), info = info, ci_test = "testIndClogit", final = final ) 
}    










