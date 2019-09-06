normlog.fsreg <- function(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, tol = 2, ncores = 1) {
  
  if ( !is.null(ini) ) {
    result <- normlog.fsreg_2(target, dataset, iniset = ini, wei = wei, threshold = threshold, tol = tol, ncores = ncores) 
    
  } else {  ## else do the classical forward regression
    
    threshold <- log(threshold)  ## log of the significance level
    p <- dim(dataset)[2]  ## number of variables
    devi <- dof <- phi <- numeric( p )  
    moda <- list()
    k <- 1   ## counter
    n <- length(target)  ## sample size
    con <- log(n)
    tool <- numeric( min(n, p) )
    dataset <- as.data.frame(dataset)
    
    runtime <- proc.time()
    ini <- glm( target ~ 1, family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )$deviance
    
    if (ncores <= 1) {
  		
        for (i in 1:p) {
          mi <- glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
          devi[i] <- mi$deviance
          dof[i] <- length( mi$coefficients ) 
          phi <- summary(mi)[[ 14 ]]
        }

      stat <- (ini - devi)/(dof - 1) /phi
      pval <- pf( stat, dof - 1, n - dof , lower.tail = FALSE, log.p = TRUE )
    } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
          ww <- glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
          return( c( ww$deviance, length( ww$coefficients ), summary(ww)[[ 14 ]] ) )
        }
        stopCluster(cl)

    }
    mat <- cbind(1:p, pval, stat) 
    colnames(mat) <- c( "variables", "log.p-value", "stat" )
    rownames(mat) <- 1:p
    sel <- which.min(pval)
    info <- matrix( numeric(3), ncol = 3 )
    sela <- sel
    
    if ( mat[sel, 2] < threshold ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE] 
        mi <- glm( target ~ dataset[, sel], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
        tool[1] <- BIC( mi )

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
      ini <- mi$deviance
      do <- length( mi$coefficients ) 
      devi <- dof <- phi <- numeric( pn )  
      
      if ( ncores <= 1 ) {
        devi <- dof <- numeric(pn)
          for ( i in 1:pn ) {
            ww <- glm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
            devi[i] <- ww$deviance
            dof[i] <- length( ww$coefficients )
            phi[i] <- summary(ww)[[ 14 ]]
          }
          
        stat <- (ini - devi) /(dof - do) / phi
        pval <- pf( stat, dof - do, n - dof, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- glm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
            return( c( ww$deviance, length( ww$coefficients ), summary(ww)[[ 14 ]] ) )
          }
          stopCluster(cl)

        stat <- (ini - mod[, 1]) / (mod[, 2] - do)/ mod[, 3]
        pval <- pf( stat, mod[, 2] - do, n - mod[, 2], lower.tail = FALSE, log.p = TRUE )
      }
      
      mat[, 2:3] <- cbind(pval, stat)
      ina <- which.min(mat[, 2])
      sel <- mat[ina, 1]    
      
      if ( mat[ina, 2] < threshold ) {
          ma <- glm( target ~ dataset[, sela] + dataset[, sel], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
          tool[k] <- BIC( ma )

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
    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      
      while ( ( info[k, 2] < threshold ) &  ( k < n ) & ( tool[ k - 1 ] - tool[ k ] > tol ) & ( nrow(mat) > 0 ) )  {
        
        ini <- moda[[ k ]]$deviance
        do <- length( coef( moda[[ k ]]  ) ) 
        k <- k + 1   
        pn <- p - k  + 1
        devi <- dof <- phi <- numeric( pn )  
        
        if (ncores <= 1) {  
          devi = dof = numeric(pn) 
            for ( i in 1:pn ) {
              ma <- glm( target ~., data = dataset[, c(sela, mat[i, 1] ) ], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
              devi[i] <- ma$deviance
              dof[i] <- length( ma$coefficients ) 
              phi[i] <- summary(ma)[[ 14 ]]
            }
          stat <- (ini - devi)/(dof - do)/phi
          pval <- pf( stat, dof - do, n - dof, lower.tail = FALSE, log.p = TRUE )
          
        } else {
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
              ww <- glm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
              return( c( ww$deviance, length( ww$coefficients ), summary(ww)[[ 14 ]] ) )
            }
            stopCluster(cl)

          stat <- (ini - mod[, 1])/(mod[, 2] - do)/mod[, 3]
          pval <- pf( stat, mod[, 2] - do, n - mod[, 2], lower.tail = FALSE, log.p = TRUE )
        }
        
        mat[, 2:3] <- cbind(pval, stat)
        ina <- which.min(mat[, 2])
        sel <- mat[ina, 1]    
        
        if ( mat[ina, 2] < threshold ) {
            ma <- glm( target ~., data = as.data.frame( dataset[, c(sela, sel) ] ), family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
            tool[k] <- BIC( ma )

          if ( tool[ k - 1 ] - tool[ k  ] <= tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )
            
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- info[, 1]
            mat <- mat[-ina, , drop = FALSE]
            moda[[ k ]] <- ma
          } 
          
        } else  info <- rbind(info, c( 1e300, 0, 0 ) )
      }
    } 
    
    runtime <- proc.time() - runtime
    final <- NULL
    
    d <- p - dim(mat)[1]
    
    if ( d >= 1 ) {
        final <- glm( target ~., data = dataset[, sela, drop = FALSE], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
      info <- info[1:d, , drop = FALSE]
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "log.p-values", "stat", "BIC" )
      rownames(info) <- info[, 1]
      mat <- cbind( mat, exp(mat[, 2]) )
      colnames(mat)[4] <- "p-value"
    }
    result = list( runtime = runtime, mat = mat, info = info, ci_test = "testIndNormLog", final = final ) 
    
  }
  
  result
}    










