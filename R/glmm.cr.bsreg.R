glmm.cr.bsreg <- function(target, dataset, id, threshold = 0.05, wei = NULL) {
  
    threshold <- log(threshold)
    dm <- dim(dataset)
    if ( is.null(dm) ) {
      n <- length(target)
      p <- 1
    } else {
      n <- dm[1]  ## sample size 
      p <- dm[2]  ## number of variables
    }  
    if ( p > n ) {
      res <- paste("The number of variables is higher than the sample size. No backward procedure was attempted")
      
    } else {
      
      tic <- proc.time()
      #check for NA values in the dataset and replace them with the variable median or the mode
      if( any(is.na(dataset)) ) {
        warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
        dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
      }
      ###################
      ###################
      ini <- coxme::coxme( target ~ dataset + (1 | id), weights = wei )
      likini <- ini$loglik[2] 
      stat <- numeric(p)
      if ( p == 1 )  {
        mod <- coxme::coxme( target ~ 1 + (1 | id), weights = wei )
        stat <- 2 * ( likini - mod$loglik[1] )	
      } else {
        for (j in 1:p) {
          mod <- coxme::coxme( target ~ dataset[, -j, drop = FALSE] + (1 | id), weights = wei )
          stat[j] <- 2 * ( likini - mod$loglik[2] )
        }
      }  
      mat <- cbind(1:p, pchisq( stat, 1, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 
      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10), ncol = 3 )
      
      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, final = ini ) 
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE] 
        
        i <- 1  
        
        if ( info[1, 2] > threshold & dim(mat)[1] > 0) {
          
          while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
            
            ini <- coxme::coxme( target ~ dat + (1 | id), weights = wei )
            likini <- ini$loglik[2] 
            i <- i + 1        
            k <- p - i + 1
            
            if ( k == 1 ) {
              mod <- coxme::coxme(target ~ 1 + (1 | id), REML = FALSE, weights = wei)
              stat <- 2 * ( likini - mod$loglik[1] )
              pval <- pchisq( stat, 1, lower.tail = FALSE, log.p = TRUE)
              
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], pval, stat) )
                dat <- dataset[, -info[, 1], drop = FALSE ]
                mat <- matrix(nrow = 0, ncol = 3)
              } else {
                info <- rbind(info, c(0, -10, -10)) 
                final <- ini
              }  
              
            } else { 
              stat <- numeric(k)
              
              for (j in 1:k) {
                mod <- coxme::coxme( target ~ dat[, -j, drop = FALSE] + (1 | id), weights = wei )
                stat[j] <- 2 * ( likini - mod$loglik[2] )
              }
              mat[, 2:3] <- cbind( pchisq( stat, 1, lower.tail = FALSE, log.p = TRUE), stat )
              sel <- which.max( mat[, 2] )
              
              if ( mat[sel, 2] < threshold ) {
                final <- ini
                info <- rbind(info, c(0, -10, -10) )
                
              } else {
                info <- rbind(info, mat[sel, ] )
                mat <- mat[-sel, , drop = FALSE] 
                dat <- dataset[, -info[, 1], drop = FALSE ]
              }
              
            }
          }  ## end while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 ) 
          runtime <- proc.time() - tic
          info <- info[ info[, 1] > 0, , drop = FALSE]
          colnames(mat) <- c("Variables", "log.p-values", "statistic")
          res <- list(runtime = runtime, info = info, mat = mat, final = final ) 
          
        } else {
          runtime <- proc.time() - tic
          res <- list(runtime = runtime, info = info, mat = NULL, final = mod ) 
        }  ## end if ( info[1, 2] > threshold & dim(mat)[1] > 0) 
        
      }  ## end  if ( mat[sel, 2] < threshold ) 
    }  ## end if ( p > n )
    
  res
}  


