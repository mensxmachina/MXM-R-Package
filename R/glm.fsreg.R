glm.fsreg <- function(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, tol = 2, ncores = 1) {
  
  ## target can be Real valued (normal), binary (binomial) or counts (poisson)
  ## dataset is a matrix or a data.frame with the predictor variables
  ## test is the test used, but not really required because depending on the data this will be decided.
  ## there is no hrm is psecifying this though
  ## threshold is the level of significance
  ## method can be either BIC or adjrsq (for non-robust linear models only). The BIC is a consistent method for selecting
  ## models, thus we also use it to avoid overfitting
  ## stopping is based on "BIC"
  ## tol is the tolerance value for the method. If BIC is used as the stopping rule, the default is 2, but usually can be 2 or 4.
  ## If BIC is used as a way to proceed, the tol is 0.
  ## ncores is for parallel processing 
  ###### If there is an initial set of variables do this function
  if ( !is.null(ini) ) {
    result <- glm.fsreg_2(target, dataset, iniset = ini, wei = wei, threshold = threshold, tol = tol, ncores = ncores) 
    
  } else {  ## else do the classical forward regression
    
    threshold <- log(threshold)  ## log of the significance level
    p <- dim(dataset)[2]  ## number of variables
    devi <- dof <- numeric( p )  
    moda <- list()
    k <- 1   ## counter
    n <- length(target)  ## sample size
    tool <- numeric( min(n, p) )
    con <- length(target)
	
    #########
    ## if it is binomial regression
    #########
    
    if ( is.matrix(target)  &  NCOL(target) == 2 )  {

      dataset <- as.data.frame(dataset)
        
        ci_test <- "testIndBinom"
        runtime <- proc.time()
        wei <- target[, 2]
        y <- target[, 1] / wei
        devi <- dof <- numeric(p)
        ini <- glm( y ~ 1, weights = wei, family = binomial, y = FALSE, model = FALSE )$deviance  ## residual deviance

        if (ncores <= 1) {
          for (i in 1:p) {
            mi <- glm( y ~ dataset[, i], weights = wei, family = binomial, y = FALSE, model = FALSE )
            devi[i] <- mi$deviance
            dof[i] <- length( mi$coefficients ) 
          }

          stat <- ini - devi
          pval <- pchisq( stat, dof - 1, lower.tail = FALSE, log.p = TRUE )
          
        } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
            ww <-  glm( y ~ dataset[, i], weights = wei, family = binomial, y = FALSE, model = FALSE )
            return( c( ww$deviance, length( ww$coefficients ) ) )
          }
          stopCluster(cl)
          stat <- ini - mod[, 1]
          pval <- pchisq( stat, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE )
        }
        
        mat <- cbind(1:p, pval, stat) 
        colnames(mat) <- c( "variables", "log.p-value", "stat" )
        rownames(mat) <- 1:p
        
        sel <- which.min(pval)
        info <- matrix( numeric(3), ncol = 3 )
        sela <- sel
        
        if ( mat[sel, 2] < threshold ) {
          info[1, ] <- mat[sel, , drop = FALSE]
          mat <- mat[-sel, ,drop = FALSE] 
          if ( !is.matrix(mat) )  mat <- matrix(mat, ncol = 3) 
          mi <- glm( y ~ dataset[, sel], weights = wei, family = binomial, y = FALSE, model = FALSE )
          tool[1] <- BIC( mi )
          moda[[ 1 ]] <- mi
        }  else  {
          info <- info  
          sela <- NULL
        }
        ########
        #####   k equals 2
        ######## 
        if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
          
          k <- 2
          pn <- p - k + 1   
          ini <- moda[[ 1 ]]$deviance  ## residual deviance
          do <- length( coef( moda[[ 1 ]]  ) ) 
          devi <- dof <- numeric(pn)
          
          if ( ncores <= 1 ) {
            for ( i in 1:pn ) {
              ww <- glm( y ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei, family = binomial, y = FALSE, model = FALSE )
              devi[i] <- ww$deviance
              dof[i] <- length( ww$coefficients )         
            }
            stat = ini - devi
            pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
            
          } else {
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
              ww <- glm( y ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei, family = binomial, y = FALSE, model = FALSE )
              return( c( ww$deviance, length( ww$coefficients ) ) )
            }
            stopCluster(cl)
            stat <- ini - mod[, 1]
            pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
            
          }
          
          mat[, 2:3] <- cbind(pval, stat)
          ina <- which.min(mat[, 2])
          sel <- mat[ina, 1]    
          
          if ( mat[ina, 2] < threshold ) {
            ma <- glm( y ~ dataset[, sela] + dataset[, sel], weights = wei, family = binomial, y = FALSE, model = FALSE )
            tool[k] <- BIC( ma )
            if ( tool[ k - 1 ] - tool[ k ] <= tol ) {
              info <- info
              
            } else { 
              info <- rbind(info, c( mat[ina, ] ) )
              sela <- info[, 1]
              mat <- mat[-ina , , drop = FALSE] 
              moda[[ k ]] <- ma
            }
            
          } else  info <- info
        }
        ########
        #####     k greater than 2
        ######## 
        if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
          
          while ( info[k, 2] < threshold &  k < n - 15 &  tool[ k - 1 ] - tool[ k ] > tol  &  nrow(mat) > 0  )  {
            
            ini <- moda[[ k ]]$deviance  ## residual deviance
            do <- length( coef( moda[[ k ]]  ) ) 
            k <- k + 1   
            pn <- p - k  + 1
            devi <- dof <- numeric( pn )  
            
            if (ncores <= 1) {  
              devi <- dof <- numeric(pn) 
              for ( i in 1:pn ) {
                ma <- glm( y ~., data = dataset[, c(sela, mat[i, 1] ) ], weights = wei, family = binomial, y = FALSE, model = FALSE )
                devi[i] <- ma$deviance
                dof[i] <- length( coef( ma ) ) 
              }

              stat <- ini - devi
              pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
              
            } else {
              cl <- makePSOCKcluster(ncores)
              registerDoParallel(cl)
              devi <- dof <- numeric(pn)
              mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
                ww <- glm( y ~., data = dataset[, c(sela, mat[i, 1] ) ], weights = wei, family = binomial, y = FALSE, model = FALSE )
                return( c( ww$deviance, length( ww$coefficients ) ) )
              }
              stopCluster(cl)
              stat <- ini - mod[, 1]
              pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
            }
            
            mat[, 2:3] <- cbind(pval, stat)
            ina <- which.min(mat[, 2])
            sel <- mat[ina, 1]    
            
            if ( mat[ina, 2] < threshold ) {
              ma <- glm( y ~., data = dataset[, c(sela, sel) ], weights = wei, family = binomial, y = FALSE, model = FALSE )
              tool[k] <- BIC( ma )

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
        
        d <- p - nrow(mat)
        final <- NULL

        if ( d >= 1 ) {
          final <- glm( y ~., data = dataset[, sela, drop = FALSE], weights = wei, family = binomial, y = FALSE, model = FALSE )
          info <- info[1:d, ]
          if ( d == 1 )  info <- matrix(info, nrow = 1)
          info <- cbind( info, tool[ 1:d ] ) 
          colnames(info) <- c( "variables", "log.p-values", "stat", "BIC" )
          rownames(info) <- info[, 1]
        }
        
        mat <- cbind( mat, exp(mat[, 2]) )
        colnames(mat)[4] <- "p-values"
        result <- list(runtime = runtime, mat = t(mat), info = info, ci_test = ci_test, final = final ) 
        
    } else  {
      
      #############################
      #####    if it is binary or poisson regression
      #############################
      
      la <- length( unique(target) )
      
      if ( is.matrix(dataset) & is.null(wei) ) {
        
        if ( la == 2 ) {
          runtime <- proc.time()
          if ( is.factor(target) )  target <- as.numeric(target) - 1
          info <- Rfast::fs.reg(target, dataset, sig = exp(threshold), tol = tol, type = "logistic")
          runtime <- proc.time() - runtime
          ci_test <- "testIndLogistic"
          d <- dim(info)[1]
          sela <- info[, 1]
          if (d == 0) {
            final <- NULL
          } else  final <- glm(target ~ dataset[, info[, 1]], binomial)
          mat <- NULL
        } else {
          runtime <- proc.time()
          info <- Rfast::fs.reg(target, dataset, sig = exp(threshold), tol = tol, type = "poisson")
          runtime <- proc.time() - runtime
          ci_test <- "testIndPois"
          d <- dim(info)[1]
          sela <- info[, 1]
          if (d == 0) {
            final <- NULL
          } else  final <- glm(target ~ dataset[, info[, 1]], poisson)
          mat <- NULL
        }
        result = list( runtime = runtime, mat = mat, info = info, ci_test = ci_test, final = final ) 
        
      } else {

     	dataset <- as.data.frame(dataset)

        if ( la == 2 ) {
          
          if ( is.factor(target) )  target <- as.numeric(target) - 1
          oiko <- binomial(logit)  ## binomial regression
          ci_test <- "testIndLogistic"
        } else {
          oiko <- poisson(log)  ## poisson regression
          ci_test <- "testIndPois"
        } 
        
        runtime <- proc.time()
        
        devi <- numeric(p)
        dof <- devi + 2
        ini <- glm( target ~ 1, family = oiko, weights = wei, y = FALSE, model = FALSE )$deviance  ## residual deviance
        
        if (ncores <= 1) {
          for (i in 1:p) {
            mi <- glm( target ~ dataset[, i], family = oiko, weights = wei, y = FALSE, model = FALSE )
            devi[i] <- mi$deviance
            dof[i] <- length( coef( mi ) ) 
          }
          stat <- ini - devi
          pval <- pchisq( stat, dof - 1, lower.tail = FALSE, log.p = TRUE )
        } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mata <- matrix(0, p, 2)
          mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
            ww <- glm( target ~ dataset[, i], family = oiko, weights = wei, y = FALSE, model = FALSE )
            mata[i, ] <- c( ww$deviance, length( coef( ww ) )  )
          }
          stopCluster(cl)
          stat = ini - mod[, 1]
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
          mat <- mat[-sel, , drop = FALSE] 
          mi <- glm( target ~ dataset[, sel], family = oiko, weights = wei, y = FALSE, model = FALSE )
          tool[1] <- BIC( mi )
          moda[[ 1 ]] <- mi
        }  else  {
          info <- info  
          sela <- NULL
        }
        ##########
        #####   k equals 2
        ########## 
        if ( info[1, 2] < threshold  &  nrow(mat) > 0 ) {
          
          k <- 2
          pn <- p - k + 1   
          devi <- dof <- numeric( pn )  
          
          ini <- mi$deviance  ## residual deviance
          do <- length( coef( mi ) ) 
          
          if ( ncores <= 1 ) {
            devi <- dof <- numeric(pn)
            for ( i in 1:pn ) {
              ww <- glm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = oiko, weights = wei, y = FALSE, model = FALSE )
              devi[i] <- ww$deviance
              dof[i] <- length( coef( ww ) )          
            }
            stat <- ini - devi
            pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
           
          } else {
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
              ww <- glm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = oiko, weights = wei, y = FALSE, model = FALSE )
              return( c( ww$deviance, length( ww$coefficients ) ) ) 
            }
            stopCluster(cl)
            stat <- ini - mod[, 1]
            pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
          }
          
          mat[, 2:3] <- cbind(pval, stat)
          ina <- which.min(mat[, 2])
          sel <- mat[ina, 1]    
          
          if ( mat[ina, 2] < threshold ) {
             ma <- glm( target ~ dataset[, sela] + dataset[, sel], family = oiko, weights = wei, y = FALSE, model = FALSE )
             tool[k] <- BIC( ma )
            if ( tool[ k - 1 ] - tool[ k ] <= tol ) {
              info <- rbind(info, c( 1e300, 0, 0 ) )
            } else { 
              info <- rbind(info, c( mat[ina, ] ) )
              sela <- info[, 1]
              mat <- mat[-ina , ,drop = FALSE ] 
              moda[[ k ]] <- ma
            }
            
          } else   info <- rbind(info, c( 1e300, 0, 0 ) )
        }
        
        #######
        ####   k greater than 2
        ####### 
        
        if ( info[k, 2] > 0  &  nrow(mat) > 0 ) {
          
          while ( ( info[k, 2] < threshold ) &  ( k < n ) & ( tool[ k - 1 ] - tool[ k ] > tol ) & ( nrow(mat) > 0 ) )  {
            
            ini = moda[[ k ]]$deviance  ## residual deviance
            do = length( coef( moda[[ k ]]  ) ) 
            k <- k + 1   
            pn <- p - k  + 1
            devi <- dof <- numeric( pn )  
            
            if (ncores <= 1) {  
              devi = dof = numeric(pn) 
              for ( i in 1:pn ) {
                ma <- glm( target ~., data = dataset[, c(sela, mat[i, 1] ) ], family = oiko, weights = wei, y = FALSE, model = FALSE )
                devi[i] <- ma$deviance
                dof[i] = length( ma$coefficients ) 
              }
              stat = ini - devi
              pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
              
            } else {
              cl <- makePSOCKcluster(ncores)
              registerDoParallel(cl)
              mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
                ww <- glm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, weights = wei, y = FALSE, model = FALSE )
                return( c( ww$deviance, length( coef( ww ) ) ) )
              }
              stopCluster(cl)
              stat <- ini - mod[, 1]
              pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
            }
            
            mat[, 2:3] <- cbind(pval, stat)
            ina <- which.min(mat[, 2])
            sel <- mat[ina, 1]    
            
            if ( mat[ina, 2] < threshold ) {
              ma <- glm( target ~., data = dataset[, c(sela, sel) ], family = oiko, weights = wei, y = FALSE, model = FALSE )
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
        d <- p - dim(mat)[1]
        final <- NULL
        if ( d >= 1 ) {
          final <- glm( target ~., data = dataset[, sela, drop = FALSE ], family = oiko, weights = wei, y = FALSE, model = FALSE )
          info <- info[1:d,  , drop = FALSE]
          info <- cbind( info, tool[ 1:d ] ) 
          colnames(info) <- c( "variables", "log.p-values", "stat", "BIC" )
          rownames(info) <- info[, 1]
          mat <- cbind( mat, exp(mat[, 2]) )
          colnames(mat)[4] <- "p-value"
        }
        result = list( runtime = runtime, mat = mat, info = info, ci_test = ci_test, final = final ) 
      }
      
    }  
    
  }
  
  result
}    










