ebic.glmm.cr.bsreg <- function(target, dataset, id, wei = NULL, gam = NULL) {

    dm <- dim(dataset)
    n <- dm[1]  ## sample size 
    p <- dm[2]  ## number of variables
    if ( p > n ) {
      res <- paste("The number of variables is higher than the sample size. No backward procedure was attempted")
      
    } else {
      
      tic <- proc.time()
      
      logn <- log(n)
      if ( is.null(gam) ) {
        con <- 2 - log(p) / logn
      } else con <- 2 * gam
      if ( con < 0 )  con <- 0
      tool <- numeric(p + 1)
      
      ini <- coxme::coxme( target ~ dataset + (1|id), weights = wei )
      bic0 <-  BIC(ini)   ## initial eBIC  
      tool[1] <- bic0
      bic <- numeric(p)
      M <- dim(dataset)[2] - 1
      
      if ( M == 0 ) {
        mod <- coxme::coxme( target ~ 1 + (1|id), weights = wei)
        bic <- BIC(mod)      
        if (bic0 - bic < 0 ) {
          info <- matrix( 0, nrow = 0, ncol = 2 )
          mat <- matrix( c(1, bic - bic0), ncol = 2 )
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
          mod <- coxme::coxme( target ~ dataset[, -j, drop = FALSE] + (1|id), weights = wei)
          bic[j] <- BIC(mod) + con * lchoose(p, M) 
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
              
              ini <- coxme::coxme( target ~ dat + (1|id), weights = wei )
              M <- dim(dat)[2]
              bic0 <-  BIC(mod) + con * lchoose(p, M)
              i <- i + 1        
              
              if ( M == 1 ) {
                mod <- coxme::coxme(target ~ 1 + (1|id), weights = wei )
                bic <-  BIC(mod)
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
                  mod <- coxme::coxme( target ~ dat[, -j, drop = FALSE] + (1|id), weights = wei )
                  bic[j] <- BIC(mod) + con * lchoose(p, M)
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
      runtime <- proc.time() - tic
      res <- list(runtime = runtime, info = info, mat = mat )
    }  ##  end  if ( p > n )
    
  res
}  





