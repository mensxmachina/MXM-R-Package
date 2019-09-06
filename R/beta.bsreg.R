beta.bsreg <- function(target, dataset, threshold = 0.05, wei = NULL) {
  
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
    res <- paste("The number of variables is hiher than the sample size. No backward procedure was attempted")
    
  } else {
    
    tic <- proc.time()
    ##################################
    # target checking and initialize #
    ################################## 
      ini <- beta.reg( target,  dataset, wei = wei )
      dofini <- length( ini$be )
      likini <- ini$loglik 
      stat <- dof <- numeric(p)
      if (p == 1) {
	      if ( is.null(wei) ) {
          mod <- Rfast::beta.mle(target)
        } else mod <- betamle.wei(target, wei)
        stat <- 2 * ( likini - mod$loglik )
        dof <- dofini - length( mod$be ) 
        pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
	    } else {
        for (j in 1:p) {
          mod <- beta.reg( target, dataset[, -j], wei = wei )
          stat[j] <- 2 * ( likini - mod$loglik )
          dof[j] <- dofini - length( mod$be )
        }
      }
      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-value", "statistic" )
      rownames(mat) <- 1:p 
      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )

      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = "testIndBeta", final = ini ) 
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE ]  

      i <- 1  
      
        while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
          
          ini <- beta.reg( target, dat, wei = wei )
          likini <- ini$loglik
          dofini <- length(ini$be)
          i <- i + 1        
          k <- p - i + 1
          
          if ( k == 1 ) {
            if ( is.null(wei) ) {
              mod <- Rfast::beta.mle(target)
            } else mod <- betamle.wei(target, wei)
            
            stat <- 2 * ( likini - mod$loglik )
            dof <- dofini - length( mod$be ) 
            pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
            
            if (pval > threshold ) {
              final <- "No variables were selected"
              info <- rbind(info, c(mat[, 1], pval, stat) )
              dat <- dataset[, -info[, 1], drop = FALSE ] 
              mat <- NULL
            } else {
              info <- rbind(info, c(0, -10, -10)) 
              final <- ini
              mat[, 2:3] <- c(pval, stat)
            }
            
          } else { 
            
            stat <- dof <- numeric(k)
            for (j in 1:k) {
              mod <- beta.reg( target,  dat[, -j], wei = wei )
              stat[j] <- 2 * ( likini - mod$loglik )
              dof[j] <- dofini - length( mod$be ) 
            }
            
            mat[, 2:3] <- cbind( pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
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
        res <- list(runtime = runtime, info = info, mat = mat, ci_test = "testIndBeta", final = final ) 

    }
  }
  res
} 
  






