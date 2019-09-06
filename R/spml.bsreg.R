spml.bsreg <- function(target, dataset, threshold = 0.05) {
  
  threshold <- log(threshold)
  dm <- dim(dataset)
  if ( !is.matrix(target) )  target <- cbind( cos(target), sin(target) ) 
  if ( is.null(dm) ) {
    n <- dim(target)[1]
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
      ini <- Rfast::spml.reg( target, dataset )
      likini <- ini$loglik
      stat <- numeric(p)
      if (p == 1) {
        y <- ( atan(target[, 2]/target[, 1]) + pi * I(target[, 1] < 0) ) %% (2 * pi)
        mod <- Rfast::spml.mle(y, rads = TRUE)
        stat <- 2 * ( likini - mod$loglik )
      }	else {
        for (j in 1:p) {
          mod <- Rfast::spml.reg( target, dataset[, -j] )
          stat[j] <- 2 * ( likini - mod$loglik )
        }
      }  
      mat <- cbind(1:p, pchisq( stat, 2, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 
      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      
      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = "testIndSPML", final = ini ) 
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE] 
        
        i <- 1  
        
        while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
          
          ini <- Rfast::spml.reg( target, dat )
          likini <- ini$loglik
          i <- i + 1        
          k <- p - i + 1
          
          if ( k == 1 ) {
            y <- ( atan(target[, 2]/target[, 1]) + pi * I(target[, 1] < 0) ) %% (2 * pi)
            mod <- Rfast::spml.mle(y, rads = TRUE)
            stat <- 2 * ( likini - mod$loglik )
            pval <- pchisq( stat, 2, lower.tail = FALSE, log.p = TRUE)
            
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
              mod <- Rfast::spml.reg( target,  dat[, -j] )
              stat[j] <- 2 * ( likini - mod$loglik )
            }
            mat[, 2:3] <- cbind( pchisq( stat, 2, lower.tail = FALSE, log.p = TRUE), stat )
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
        }
        runtime <- proc.time() - tic
        info <- info[ info[, 1] > 0, , drop = FALSE]
        colnames(mat) <- c("Variables", "log.p-values", "statistic")
        res <- list(runtime = runtime, info = info, mat = mat, ci_test = "testIndSPML", final = final ) 
        
      }  

  }
  
  res
}  

