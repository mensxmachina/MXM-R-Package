tobit.fsreg <- function(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, tol = 2, ncores = 1) {
  
  if ( !is.null(ini) ) {
    result <- tobit.fsreg_2(target, dataset, iniset = ini, wei = wei, threshold = threshold, tol = tol, ncores = ncores) 
    
  } else {  ## else do the classical forward regression
    
    threshold <- log(threshold)  ## log of the significance level
    dm <- dim(dataset)
    p <- dm[2]  ## number of variables
    devi <- dof <- numeric( p )  
    moda <- list()
    k <- 1   ## counter
    n <- dm[1]  ## sample size
    con <- log(n)
    tool <- numeric( min(n, p) )
    dataset <- as.data.frame(dataset)
    
    runtime <- proc.time()
    
      mi <- survival::survreg( target ~ 1, weights = wei, dist = "gaussian" )
      ini <-  2 * as.numeric( logLik(mi) ) 

    if (ncores <= 1) {
      
      for (i in 1:p) {
         mi <- survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
         devi[i] <- 2 * as.numeric( logLik(mi) )
         dof[i] <- length( mi$coefficients ) 
      }
      stat = devi - ini
      pval = pchisq( stat, dof - 1, lower.tail = FALSE, log.p = TRUE )
	
    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mod <- foreach( i = 1:p, .combine = rbind, .export = "survival", .packages = "survreg") %dopar% {
        ww <- survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
        return( c( 2 * as.numeric( logLik(mi) ), length( ww$coefficients ) ) )
      }
      stopCluster(cl)

      stat = mod[, 1] - ini
      pval = pchisq( stat, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE )
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
        mi <- survival::survreg( target ~., data = dataset[, sel, drop = FALSE], weights = wei, dist = "gaussian" )
        tool[1] <-  - 2 * logLik(mi) + ( length( mi$coefficients ) + 1 ) * con 		
        moda[[ 1 ]] <- mi
      }  else  {
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
            ww <- survival::survreg( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei, dist ="gaussian" )
            devi[i] <- 2 * logLik(ww)  
            dof[i] <- length( ww$coefficients )          
          }
          stat <- devi - ini
          pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
          
        } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mod <- foreach( i = 1:pn, .combine = rbind, .export = "survival", .packages = "survreg") %dopar% {
            ww <- survival::survreg( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei, dist = "gaussian" )
            mata[i, ] <- c( 2 * logLik(ww), length( ww$coefficients ) )
          }
          stopCluster(cl)		  
  
          stat <- mod[, 1] - ini
          pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
        }
        
        mat[, 2:3] <- cbind(pval, stat)
        ina <- which.min(mat[, 2])
        sel <- mat[ina, 1]    
        
        if ( mat[ina, 2] < threshold ) {
            ma <- survival::survreg( target ~ dataset[, sela] + dataset[, sel], weights = wei, dist = "gaussian" )
            tool[k] <-  - 2 * logLik(ma) + ( length( ma$coefficients) + 1 ) * con		  
          
          if ( tool[ k - 1 ] - tool[ k ] <= tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )
            
          } else { 
            info <- rbind(info, c( mat[ina, ] ) )
            sela <- info[, 1]
            mat <- mat[-ina , , drop = FALSE] 
            moda[[ k ]] <- ma
          }
          
        } else   info <- rbind(info, c( 1e300, 0, 0 ) )
      }
      #######
      ####   k greater than 2
      ####### 
      if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
        
        while ( ( info[k, 2] < threshold ) &  ( k < n ) & ( tool[ k - 1 ] - tool[ k ] > tol ) & ( nrow(mat) > 0 ) )  {
          
          ini = 2 * logLik( moda[[ k ]] )
          do = length( coef( moda[[ k ]]  ) ) 
          k <- k + 1   
          pn <- p - k  + 1
          
          if (ncores <= 1) {  
            devi = dof = numeric(pn) 
            for ( i in 1:pn ) {
              ma <- survival::survreg( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1] ) ] ), weights = wei, dist = "gaussian")
              devi[i] <- 2 * logLik(ma)
              dof[i] = length( ma$coefficients ) 
            }
            stat = devi - ini
            pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
            
          } else {
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            devi = dof = numeric(pn)
            mod <- foreach( i = 1:pn, .combine = rbind, .export = "survival", .packages = "survreg") %dopar% {
              ww <- survival::survreg( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei, dist = "gaussian")
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
              ma <- survival::survreg( target ~., data = dataset[, c(sela, sel) ], weights = wei, dist = "gaussian" )
              tool[k] <- - 2 * logLik(ma) + ( length( ma$coefficients) + 1 ) * con	
 
            if ( tool[ k - 1 ] - tool[ k  ] <= tol ) {
              info <- rbind(info, c( 1e300, 0, 0 ) )
              
            } else { 
              info <- rbind( info, mat[ina, ] )
              sela <- info[, 1]
              mat <- mat[-ina , , drop = FALSE]
              if ( !is.matrix(mat) )  mat <- matrix(mat, ncol = 3) 
              moda[[ k ]] <- ma
            } 
            
          } else  info <- rbind(info, c( 1e300, 0, 0 ) )
        }
      } 
      
      runtime <- proc.time() - runtime
      d <- p - dim(mat)[1]
      final <- NULL

    if ( d >= 1 ) {
      final <- survival::survreg( target ~., data = dataset[, sela, drop = FALSE], weights = wei, dist = "gaussian" )
      info <- info[1:d, , drop = FALSE]
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "log.p-values", "stat", "BIC" )
      rownames(info) <- info[, 1]
      mat <- cbind( mat, exp(mat[, 2]) )
      colnames(mat)[4] <- "p-value"
      
    }
    result = list( runtime = runtime, mat = mat, info = info, ci_test = "testIndTobit", final = final ) 
    
  }
  
  result
}    










