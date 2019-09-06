clogit.fsreg <- function(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, tol = 2, ncores = 1) {
  
  if ( !is.null(ini) ) {
    result <- clogit.fsreg_2(target, dataset, iniset = ini, threshold = threshold, tol = tol, ncores = ncores) 
    
  } else {  ## else do the classical forward regression

    threshold <- log(threshold)  ## log of the significance level
    dataset <- as.data.frame(dataset)
    case <- as.logical(target[, 1]);  #
    id <- target[, 2] #the patient id
  
    dm <- dim(dataset)
    p <- dm[2]  ## number of variables
    devi <- dof <- numeric( p )  
    moda <- list()
    k <- 1   ## counter
    n <- dm[1]  ## sample size
    tool <- numeric( min(n, p) )
    oop <- options(warn = -1) 
    on.exit( options(oop) )
	
    runtime <- proc.time()
	  
    if (ncores <= 1) {

      for (i in 1:p) {
        mi <- survival::clogit( case ~ dataset[, i] + strata(id) )
        devi[i] <- 2 * abs( diff(mi$loglik) )
        dof[i] <- length( mi$coefficients ) 
      }
      stat <- devi 
      pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE )
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach( i = 1:p, .combine = rbind, .export = "survival", .packages = "clogit") %dopar% {
        ww <- survival::clogit( case ~ dataset[, i] + strata(id) )
        return( c( 2 * abs( diff(ww$loglik) ), length( ww$coefficients ) ) )
      }
      stopCluster(cl)
      
      stat <- mod[, 1] 
      pval <- pchisq( stat, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE )
    }
    mat <- cbind(1:p, pval, stat) 
    colnames(mat) <- c( "variables", "log.p-value", "stat" )
    rownames(mat) <- 1:p
    sel <- which.min(pval)
    info <- matrix( numeric(3), ncol = 3 )
    sela <- sel
    
    if ( mat[sel, 2] < threshold ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel,  , drop = FALSE] 
      mi <- survival::clogit( case ~ dataset[, sel] + strata(id) )
      tool[1] <-  BIC(mi)	
      moda[[ 1 ]] <- mi
    } else  {
      info <- info  
      sela <- NULL
    }
    ##########
    #####   k equals 2
    ########## 
    if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
      
      k <- 2
      pn <- p - k + 1   
      ini <- 2 * logLik(mi)   ## residual deviance
      do <- length( mi$coefficients ) 
      devi <- dof <- numeric( pn )  
      
      if ( ncores <= 1 ) {
        devi <- dof <- numeric(pn)
        for ( i in 1:pn ) {
          ww <- try( survival::clogit( case ~ dataset[, sela ] + dataset[, mat[i, 1] ] + strata(id) ), silent = TRUE )
          if ( identical( class(ww), "try-error" ) ) {
            devi[i] <- ini
            dof[i] <- do + 1
          } else {
            devi[i] <- 2 * logLik(ww)  
            dof[i] <- length( ww$coefficients )          
          }  
        }  ## end for ( i in 1:pn ) 
        stat <- devi - ini
        pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mat <- matrix(0, pn, 2)
        mod <- foreach( i = 1:pn, .combine = rbind, .export = "clogit", .packages = "survival") %dopar% {
          ww <- try( survival::clogit( case ~ dataset[, sela ] + dataset[, mat[i, 1] ] + strata(id) ), silent = TRUE )
          if ( identical( class(ww), "try-error" ) ) {
            mat[i, ] <- c(ini, do)
          } else {
            mat[i, ] <- c( 2 * logLik(ww), length( ww$coefficients ) )
          }  
        }
        stopCluster(cl)		  
        
        stat <- mod[, 1] - ini
        pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
      }
      
      mat[, 2:3] <- cbind(pval, stat)
      ina <- which.min(mat[, 2])
      sel <- mat[ina, 1]    
      
      if ( mat[ina, 2] < threshold ) {
        ma <- survival::clogit( case ~ dataset[, sela] + dataset[, sel] + strata(id) )
        tool[k] <- BIC(ma)
        
        if ( tool[ k - 1 ] - tool[ k ] <= tol ) {
          info <- rbind(info, c( 1e300, 0, 0 ) )
          
        } else { 
          info <- rbind(info, c( mat[ina, ] ) )
          sela <- info[, 1]
          mat <- mat[-ina , , drop = FALSE] 
          moda[[ k ]] <- ma
        }
        
      } else   info <- rbind(info, c( 1e300, 0, 0 ) )
    }  ## end if ( info[k, 2] < threshold  &  nrow(mat) > 0 )
    #######
    ####   k greater than 2
    ####### 
    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      
      while (  info[k, 2] < threshold  &   k < n  &  tool[ k - 1 ] - tool[ k ] > tol  &   nrow(mat) > 0 )  {
        
        ini = 2 * logLik( moda[[ k ]] )
        do = length( coef( moda[[ k ]]  ) ) 
        k <- k + 1   
        pn <- p - k  + 1
        devi <- dof <- numeric( pn )  
        
        if (ncores <= 1) {  
          devi = dof = numeric(pn) 
          for ( i in 1:pn ) {
            ma <- try( survival::clogit(  case ~ . + strata(id), data = dataset[, c(sela, mat[i, 1])] ), silent = TRUE)
          if ( identical( class(ma), "try-error" ) ) {
              devi[i] <- ini
              dof[i] <- do + 1
            } else {
              devi[i] <- 2 * logLik(ma)  
              dof[i] <- length( ma$coefficients )          
            }  
          }
          stat = devi - ini
          pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
          
        } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mat <- matrix(0, pn, 2)
          mod <- foreach( i = 1:pn, .combine = rbind, .export = "clogit", .packages = "survival") %dopar% {
            ww <- try( survival::clogit(  case ~ . + strata(id), data = dataset[, c(sela, mat[i, 1])] ), silent = TRUE)
            if ( identical( class(ww), "try-error" ) ) {
              mat[i, ] <- c(ini, do)
            } else {
              mat[i, ] <- c( 2 * logLik(ww), length( ww$coefficients ) )
            }  
          }
          stopCluster(cl)
          
          stat <- mod[, 1] - ini
          pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
        }
        
        mat[, 2:3] <- cbind(pval, stat)
        ina <- which.min(mat[, 2])
        sel <- mat[ina, 1]    
        
        if ( mat[ina, 2] < threshold ) {
          sela <- c(sela, sel)
          ma <- survival::clogit(  case ~ . + strata(id), data = dataset[, sela] )
          tool[k] <-  BIC(ma)	
          
          if ( tool[ k - 1 ] - tool[ k  ] <= tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )
            
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- info[, 1]
            mat <- mat[-ina , , drop = FALSE]
            moda[[ k ]] <- ma
          } 
          
        } else  info <- rbind(info, c( 1e300, 0, 0 ) )
      }  ## end while (  info[k, 2] < threshold  &   k < n  &  tool[ k - 1 ] - tool[ k ] > tol  &   nrow(mat) > 0 ) 
    }  ## end if ( nrow(info) > 1  &  nrow(mat) > 0 )
    
    runtime <- proc.time() - runtime
    d <- p - dim(mat)[1]
    final <- NULL
    
    if ( d >= 1 ) {
      final <- survival::clogit(  case ~ . + strata(id), data = dataset[, sela, drop = FALSE] )
      info <- info[1:d, , drop = FALSE]
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "log.p-values", "stat", "BIC" )
      rownames(info) <- info[, 1]
      mat <- cbind( mat, exp(mat[, 2]) )
      colnames(mat)[4] <- "p-value"
      
    }
    result = list( runtime = runtime, mat = mat, info = info, ci_test = "testIndClogit", final = final ) 
    
  }
  
  result
}    










