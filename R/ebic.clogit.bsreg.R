ebic.clogit.bsreg <- function(target, dataset, wei = NULL, gam = NULL) {
  
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
    case <- target[, 1]
    id <- target[, 2]
     
    ini <- survival::clogit( case ~ . + strata(id),  data = dataset )
    bic0 <-  BIC(ini)    ## initial bic  
    tool[1] <- bic0
    bic <- numeric(p)
    M <- dim(dataset)[2] - 1
    
    if ( M == 0 ) {
      info <- matrix( c(1, bic0), ncol = 2 )
      mat <- matrix(0, nrow = 0, ncol = 2 )
      runtime <- proc.time() - tic
      colnames(info) <- c("Variables", "eBIC")
      colnames(mat) <- c("Variables", "eBIC")
      res <- list(runtime = runtime, info = info, mat = mat )
      
    } else { 
      ########### 
      for (j in 1:p) {
        mod <- survival::clogit( case ~ . + strata(id),  data = dataset[, -j, drop = FALSE] )
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
            
            ini <- survival::clogit( case ~ . + strata(id), data = dat )
            M <- dim(dat)[2]
            bic0 <-  BIC(ini) + con * lchoose(p, M)
            i <- i + 1        
            
            if ( M == 1 ) {
              tool[i] <- BIC(ini)
              runtime <- proc.time() - tic		
              info <- rbind(info, c(mat[, 1], bic) )
              mat <- mat[-1, , drop = FALSE]
              res <- list(runtime = runtime, info = info, mat = mat )
              dat <- dataset[, -info[, 1], drop = FALSE ]
              
            } else { 
              bic <- numeric(M)
              M <- dim(dat)[2] - 1
              for ( j in 1:(M + 1) ) {
                mod <- survival::clogit( case ~ . + strata(id), data = dat[, -j, drop = FALSE] )
                bic[j] <-  BIC(mod) + con * lchoose(p, M)
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





