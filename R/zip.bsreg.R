zip.bsreg <- function(target, dataset, threshold = 0.05, wei = NULL) {
  
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
    
    runtime <- proc.time()
    ##################################
    # target checking and initialize #
    ################################## 
    ini <- zip.mod( target,  dataset, wei = wei )
    dofini <- length( ini$be[, 1] )
    likini <- ini$loglik 
    stat <- dof <- numeric(p)
    if ( is.null(wei) ) {
      lgy <- sum( lgamma(target + 1) )  
    } else  lgy <- sum( wei * gamma(target + 1) )  

    if ( p == 1) {
      if ( is.null(wei) ) {
        mod <- Rfast::zip.mle(target)
        lgy <- sum( lgamma(target + 1) )  
      } else {
        mod <- zipmle.wei(target, wei)
        lgy <- sum( wei * gamma(target + 1) )  
      }  
      stat <- 2 * ( likini - mod$loglik )
      dof <- dofini - 1 
      pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)      
    } else {
      for (j in 1:p) {
        mod <- zip.reg( target, dataset[, -j], wei = wei, lgy = lgy )
        stat[j] <- 2 * ( likini - mod$loglik )
        dof[j] <- dofini - length( mod$be )
      }  
    }
    mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
    colnames(mat) <- c("variable", "p-value", "statistic" )
    rownames(mat) <- 1:p 
    sel <- which.max( mat[, 2] )
    info <- matrix( c(0, -10, -10) , ncol = 3 )
    sela <- sel 
    
    if ( mat[sel, 2] < threshold ) {
      runtime <- runtime - proc.time()
      res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = "testIndZIP", final = ini ) 
      
    } else {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE] 
      dat <- dataset[, -sel, drop = FALSE] 
    
    i <- 1  
    
      while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
        
        ini <- zip.mod( target, dat, wei = wei )
        likini <- ini$loglik
        dofini <- length(ini$be[, 1])
        i <- i + 1        
        k <- p - i + 1
        
        if ( k == 1 ) {
          if ( is.null(wei) ) {
            mod <- Rfast::zip.mle(target)
          } else mod <- zipmle.wei(target, wei)
          stat <- 2 * ( likini - mod$loglik )
          dof <- dofini - length( mod$be ) 
          pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
          
          if (pval > threshold ) {
            final <- "No variables were selected"
            info <- rbind(info, c(mat[, 1], pval, stat) )
            dat <- as.data.frame( dataset[, -info[, 1] ] )
            mat <- NULL
          } else {
            info <- rbind(info, c(0, -10, -10)) 
            final <- ini
            mat[, 2:3] <- c(pval, stat)
          }     
          
        } else { 
          
          stat <- dof <- numeric(k)
          for (j in 1:k) {
            mod <- zip.reg( target,  dat[, -j], wei = wei, lgy = lgy )
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
          } ## end if ( mat[sel, 2] < threshold )
          
        }
      }
  
      runtime <- proc.time() - runtime	
      info <- info[ info[, 1] > 0, , drop = FALSE]
      res <- list(runtime = runtime, info = info, mat = mat, ci_test = "testIndzIP", final = final ) 

    }  
  }  
  
  res
} 







