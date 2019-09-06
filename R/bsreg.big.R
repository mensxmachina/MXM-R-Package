bsreg.big <- function(target, dataset, threshold = 0.01, test = "testIndLogistic") {
  
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
    ind <- 1:p
    ini <- big.model(y = target, x = dataset, test = test)
    inidevi <- ini$devi
    phi <- ini$phi
    if ( is.na(phi) )  phi <- 1
    if (p == 1) {
      devi <- big.model( y = target, x = NULL, test = test)$devi
      stat <- (devi - inidevi)/phi
    }	else {
      devi <- numeric(p)
      for (j in 1:p)   devi[j] <- big.model( y = target,  x = dataset[, -j], test = test )$devi
      stat <- (devi - inidevi)/phi
    }  
    mat <- cbind(1:p, pf( stat, 1, n - p, lower.tail = FALSE, log.p = TRUE), stat )
    colnames(mat) <- c("variable", "log.p-values", "statistic" )
    rownames(mat) <- 1:p 
    sel <- which.max( mat[, 2] )
    info <- matrix( c(0, -10, -10) , ncol = 3 )
      
      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = test, final = ini ) 
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        ind[sel] <- 0
        final <- "No variables were selected"
 
        i <- 1  
        
          while ( info[i, 2] > threshold  &  sum(ind) > 0 )  {   
            
            ini <- big.model( y = target, x = dataset[, ind[ ind > 0] ], test = test )
            inidevi <- ini$devi
            phi <- ini$phi
            if ( is.na(phi) )  phi <- 1
            i <- i + 1        
            k <- p - i + 1
            
            if ( k == 1 ) {
              mod <- big.model(y = target, x = NULL, test = test)
              stat <- (mod$devi - inidevi)/phi
              pval <- pf( stat, 1, n - 2, lower.tail = FALSE, log.p = TRUE)
              
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
              devi <- numeric(k)
              ela <- ind[ ind > 0]
              ela <- matrix( rep(ela, length(ela) ), ncol = length(ela) )
              diag(ela) <- 0
              
              for (j in ela ) {
                b <- ela[, j]
                devi[j] <- big.model( y = target, x = dataset[, b[b > 0] ], test = test )$devi
              }
              stat <- (devi - inidevi)/phi
              mat[, 2:3] <- cbind( pf( stat, 1, n - k - 1, lower.tail = FALSE, log.p = TRUE), stat )
              sel <- which.max( mat[, 2] )
              
              if ( mat[sel, 2] < threshold ) {
                final <- ini
                info <- rbind(info, c(0, -10, -10) )
                
              } else {
                info <- rbind(info, mat[sel, ] )
                mat <- mat[-sel, , drop = FALSE] 
                ind[sel] <- 0
              }
              
            }  ##  end  if ( k == 1 )
          }  ##  end  while ( info[i, 2] > threshold  &  sum(ind) > 0 )  { 
          runtime <- proc.time() - tic
          info <- info[ info[, 1] > 0, , drop = FALSE]
          colnames(mat) <- c("Variables", "log.p-values", "statistic")
          res <- list(runtime = runtime, info = info, mat = mat, ci_test = test, final = final ) 
        
      }  ##  end  if ( mat[sel, 2] < threshold ) {
    }  ##  end  if ( p > n ) {

  res
}  





