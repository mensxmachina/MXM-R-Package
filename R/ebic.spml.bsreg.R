ebic.spml.bsreg <- function(target, dataset, gam = NULL) {
  
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
    
    logn <- log(n)
    if ( is.null(gam) ) {
      con <- 2 - log(p) / logn
      if ( (con) < 0 )  con <- 0
    } else con <- 2 * gam
    tool <- numeric(p + 1)
    
    ini <- Rfast::spml.reg( target, dataset )
    bic0 <-  - 2 * ini$loglik + length(ini$be) * logn    ## initial BIC  
    tool[1] <- bic0
    bic <- numeric(p)
    M <- dim(dataset)[2] - 1
    
    if ( M == 0 ) {
        if ( is.matrix(target) ) y <- ( atan(y[, 2]/y[, 1]) + pi * I(y[, 1] < 0) ) %% (2 * pi)
        bic <-  - 2 * Rfast::spml.mle(y)$loglik + 2 * logn
      
      if (bic0 - bic < 0 ) {
        info <- matrix( 0, nrow = 0, ncol = 2 )
        mat <- matrix( c(1, bic), ncol = 2 )
      } else { 
        info <- matrix( c(1, bic), ncol = 2 )
        mat <- matrix(0, nrow = 0, ncol = 2 )
      }
      runtime <- proc.time() - tic
      colnames(info) <- c("Variables", "eBIC")
      colnames(mat) <- c("Variables", "eBIC")
      res <- list(runtime = runtime, info = info, mat = mat )
      
    } else { 
      ###########  
      for (j in 1:p) {
        mod <- Rfast::spml.reg( target, dataset[, -j, drop = FALSE] )
        bic[j] <-  - 2 * mod$loglik + length(mod$be) * logn + con * lchoose(p, M)
      }
      
      mat <- cbind(1:p, bic )
      sel <- which.min( mat[, 2] )
      info <- matrix( c(0, 0), ncol = 2 )
      colnames(info) <- c("Variables", "eBIC")
      colnames(mat) <- c("Variables", "eBIC")
      
      if ( bic0 - mat[sel, 2] < 0  ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = info, mat = mat )
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE] 
        tool[2] <- info[1, 2]
        
        i <- 2  
        
        if ( tool[2] != 0 ) {
          
          while ( tool[i - 1] - tool[i ] > 0  &  NCOL(dat) > 0 )  {   
            
            ini <- Rfast::spml.reg( target, dat )
            M <- dim(dat)[2]
            bic0 <-  - 2 * ini$loglik + length(ini$be) * logn + con * lchoose(p, M)
            i <- i + 1        
            
            if ( M == 1 ) {
              
              if ( is.matrix(target) ) y <- ( atan(y[, 2]/y[, 1]) + pi * I(y[, 1] < 0) ) %% (2 * pi)
              bic <-  - 2 * Rfast::spml.mle(y)$loglik + 2 * logn
              tool[i] <- bic
              if (bic0 - bic < 0 ) {
                runtime <- proc.time() - tic
                res <- list(runtime = runtime, info = info, mat = mat )
              } else {
                runtime <- proc.time() - tic		
                info <- rbind(info, c(mat[, 1], bic) )
                mat <- mat[-1, , drop = FALSE]
                res <- list(runtime = runtime, info = info, mat = mat )
                dat <- dataset[, -info[, 1], drop = FALSE ]
              }  
              
            } else { 
              bic <- numeric(M)
              M <- dim(dat)[2] - 1
              for ( j in 1:(M + 1) ) {
                mod <- Rfast::spml.reg( target, dat[, -j, drop = FALSE] )
                bic[j] <-  - 2 * mod$loglik + length(mod$be) * logn + con * lchoose(p, M)
              }  
              mat[, 2] <- bic
              sel <- which.min( mat[, 2] )
              tool[i] <- mat[sel, 2]
              if ( bic0 - mat[sel, 2] < 0 ) {
                runtime <- proc.time() - tic
                res <- list(runtime = runtime, info = info, mat = mat )
              } else {
                info <- rbind(info, mat[sel, ] )
                mat <- mat[-sel, , drop = FALSE] 
                dat <- dataset[, -info[, 1], drop = FALSE ]
              }  ## end if ( bic0 - mat[sel, 2] < 0 )
            }  ## end if (M == 1)
          }  ## end while
          
        }  ## end if ( tool[2] > 0 )
      }  ## end if ( bic0 - mat[sel, 2] < 0  ) 
    }  ## end  if (M == 0)
  }  ##  end  if ( p > n )
  
  res
}  





