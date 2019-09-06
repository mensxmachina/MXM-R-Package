mm.fsreg <- function(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, stopping = "BIC", tol = 2, ncores = 1 ) {  
  
  ###### If there is an initial set of variables do this function
  if ( !is.null(ini) ) {
      result <- mm.fsreg_2(target, dataset, iniset = ini, threshold = threshold, wei = NULL, stopping = stopping, tol = tol, ncores = ncores) 
      
  } else {  ## else do the classical forward regression
  
  threshold <- log(threshold)
  p <- dim(dataset)[2]  ## number of variables
  pval <- stat <- dof <- numeric( p )  
  moda <- list()
  k <- 1   
  n <- length(target)  ## sample size
  con <- log(n)
  tool <- numeric( min(n, p) )
  ci_test <- "testIndReg"
  runtime <- proc.time()
    
    if (ncores <= 1) {
       for (i in 1:p) {
        ww = MASS::rlm( target ~ dataset[, i ], maxit = 2000, method = "MM") 
        stat[i] = 2 * as.numeric( logLik(ww) )
        dof[i] = length( coef(ww) )
       }

       fit0 = MASS::rlm( target ~ 1, maxit = 2000, method = "MM") 
       stat0 = 2 * as.numeric( logLik(fit0) )
       difa = abs( stat - stat0 )
       pval = pchisq(difa, dof - 1, lower.tail = FALSE, log.p = TRUE)
       mat <- cbind(1:p, pval, difa)

    } else {
        fit0 = MASS::rlm( target ~ 1, maxit = 2000, weights = wei, method = "MM" ) 
        stat0 = 2 * logLik(fit0)
         cl <- makePSOCKcluster(ncores)
         registerDoParallel(cl)
         mod <- foreach( i = 1:p, .combine = rbind, .export = "rlm", .packages = "MASS" ) %dopar% {
           ww = MASS::rlm( target ~ dataset[, i], maxit = 2000, method = "MM") 
           return( c( 2 * as.numeric( logLik(ww) ), length( coef(ww) ) ) )
         }
         stopCluster(cl)
         
         difa = abs( mod[, 1] - stat0 )
         pval = pchisq(difa, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE)
         mod = cbind( pval, difa)
         mat <- cbind(1:p, mod)      
    }
    
    colnames(mat) <- c( "variables", "log.p-value", "stat" )
    rownames(mat) <- 1:p
    sel <- which.min(mat[, 2])
    info <- matrix( numeric(3), ncol = 3 )
    sela <- sel
    dataset <- as.data.frame(dataset)
    
    if ( mat[sel, 2] < threshold ) {
      info[1, ] <- mat[sel, , drop = FALSE]
      mat <- mat[-sel, , drop = FALSE] 

      if ( stopping == "adjrsq" ) {
          ma = MASS::rlm( target ~ dataset[, sel], maxit = 2000, method = "MM")
          r2 = cor( target, fitted(ma) )^2
          tool[1] = 1 - (1 - r2) * (n - 1) / ( n - length( coef(ma) ) - 1 )
                 
      } else if ( stopping  == "BIC" ) {
          ma = MASS::rlm( target ~ dataset[, sel], maxit = 2000, method = "MM")
          tool[1] <- BIC(ma)

      }
      moda[[ 1 ]] <- ma
      
    }  else  {
      info <- info  
      sela <- NULL
    }    
     ######
     ####   k equal to 2
     ######
    if ( info[k, 2] < threshold  &  nrow(mat) > 0 )  {
      
      k <- k + 1
      pn <- p - k + 1   
      
      if ( ncores <= 1 ) {
          do = 2 * as.numeric( logLik( moda[[ 1 ]] ) )
          fr = length( coef( moda[[ 1 ]] ) )
          sta = dof = numeric(pn)
          for (i in 1:pn) {
            ww = MASS::rlm( target ~ dataset[, sel] + dataset[, mat[i, 1] ], maxit = 2000, method = "MM")
            sta[i] = 2 * as.numeric( logLik(ww) )
            dof[i] = length( coef(ww) )
          } 
          mat[, 3] = abs( sta - do )
          mat[, 2] = pchisq(mat[, 3], dof - fr, lower.tail = FALSE, log.p = TRUE)
        

      } else {
          do = 2 * as.numeric( logLik( moda[[ 1 ]] ) )
          fr = length( coef( moda[[ 1 ]] ) )
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mod <- foreach( i = 1:pn, .combine = rbind, .export = c("rlm"), .packages = "MASS" ) %dopar% {
            ww <- MASS::rlm( target ~ dataset[, sel] + dataset[, mat[i, 1] ], maxit = 2000, method = "MM")
            return( c( 2 * as.numeric( logLik(ww) ), length( coef(ww) ) ) )
          }
          stopCluster(cl)
          
          difa = abs( mod[, 1] - do )
          pval = pchisq(difa, mod[, 2] - fr, lower.tail = FALSE, log.p = TRUE)
          mod = cbind( pval, difa)
         
         mat <- cbind(mat[, 1], mod)   
      }

      ina <- which.min(mat[, 2])
      sel <- mat[ina, 1]     
      
        if ( stopping == "adjrsq" ) {
          if ( mat[ina, 2] < threshold ) {
        
              ma = MASS::rlm( target ~ dataset[, sela] + dataset[, sel], maxit = 2000, method = "MM")
              r2 = cor( target, fitted(ma) )^2
              tool[k] = 1 - (1 - r2) * (n - 1) / ( n - length( coef(ma) ) - 1)

            if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
              info <- info
            
            } else {  
              info <- rbind(info, mat[ina, ] )
              sela <- info[, 1]
              mat <- mat[-ina, , drop = FALSE ] 
              moda[[ k ]] <- ma
            }

          } else  info <- info

        } else if ( stopping == "BIC" ) {
		
          if ( mat[ina, 2] < threshold ) {
              ma = MASS::rlm( target ~ dataset[, sela] + dataset[, sel], maxit = 2000, method = "MM")
              tool[2] = BIC(ma)
            if ( tool[ k - 1] - tool[ k ] <= tol ) {
              info <- info
            
            } else {  
              info <- rbind(info, mat[ina, ] )
              sela <- info[, 1]
              mat <- mat[-ina , ]
              if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
              moda[[ k ]] <- ma
            }
          } else  info <- info
          
        }
   }
     ######
     ###### k greater than 2
     ######
    if ( nrow(info) > 1  &  nrow(mat) > 0 )  {
      while ( info[k, 2] < threshold  &  k < n - 15 & abs( tool[ k ] - tool[ k - 1 ] ) > tol & nrow(mat) > 0 )  {
         
        k <- k + 1   
        pn <- p - k + 1 
        
        if ( ncores <= 1 ) {

            do = 2 * as.numeric( logLik( moda[[ k - 1 ]] ) )
            fr = length( coef( moda[[ k - 1 ]] ) )
            sta = dof = numeric(pn)      
            for (i in 1:pn) {
              ww = MASS::rlm( target ~., data = dataset[, c(sela, mat[i, 1]) ], maxit = 2000, method = "MM")
              sta[i] = 2 * as.numeric( logLik(ww) )
              dof[i] = length( coef(ww) )
            } 
            mat[, 3] = abs( sta - do )
            mat[, 2] = pchisq(mat[, 3], dof - fr, lower.tail = FALSE, log.p = TRUE)

         } else {

            do = 2 * as.numeric( logLik( moda[[ k - 1 ]] ) )
            fr = length( coef( moda[[ k - 1 ]] ) )
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            mod <- foreach( i = 1:pn, .combine = rbind, .export = c("rlm"), .packages = "MASS" ) %dopar% {
              ww <- MASS::rlm( target ~., data = dataset[, c(sela, mat[ i, 1]) ], maxit = 2000, method = "MM")
              return( c( 2 * as.numeric( logLik(ww) ), length( coef(ww) ) ) )
            }
            stopCluster(cl)  
            difa = abs( mod[, 1] - do )
            pval = pchisq(difa, mod[, 2] - fr, lower.tail = FALSE, log.p = TRUE)
            mod = cbind( pval, difa)

          mat <- cbind( mat[, 1], mod )   

        }
       
        ina <- which.min(mat[, 2])
        sel <- mat[ina, 1]   
        
          if ( stopping == "BIC" ) {
		  
            if ( mat[ina, 2] < threshold ) {
              ma = MASS::rlm( target ~., data = dataset[, c(sela, sel)], maxit = 2000, method = "MM")
              tool[k] =  BIC(ma)
            if ( tool[ k - 1] - tool[ k ] <= tol ) {
                info <- rbind(info, c( Inf, 0, 0 ) )
              
            } else { 
              info <- rbind( info, mat[ina, ] )
              sela <- info[, 1]
              mat <- mat[-ina, , drop = FALSE] 
              moda[[ k ]] <- ma
            }
             
            } else   info <- rbind(info, c( Inf, 0, 0 ) )

          } else if ( stopping == "adjrsq" ) {
            if ( mat[ina, 2] < threshold ) {

                ma = MASS::rlm( target ~., data = dataset[, c(sela, sel)], maxit = 2000, method = "MM")
                r2 = cor(target, fitted(ma) )^2
                tool[k] <- 1 - (1 - r2) * (n - 1) / ( n - length( coef(ma) ) - 1)
             
              if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
                info <- rbind(info, c( Inf, 0, 0 ) )
              
              } else { 
                info <- rbind( info, mat[ina, ] )
                sela <- info[, 1]
                mat <- mat[-ina, , drop = FALSE]
                moda[[ k ]] <- ma
              }

            } else  info <- rbind(info, c( Inf, 0, 0 ) )
        
          }
      }

    }

    runtime <- proc.time() - runtime

    d <- length(sela)
    final <- NULL

    if ( d >= 1 ) {
      final <- MASS::rlm( target ~., data = dataset[, sela, drop = FALSE], maxit = 2000, method = "MM")  
      info <- info[1:d, , drop = FALSE]
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "log.p-value", "stat", stopping )
      rownames(info) <- info[, 1]
    }

    result <- list(runtime = runtime, mat = t(mat), info = info, ci_test = ci_test, final = final ) 
    
  }  
  
  result
}
  