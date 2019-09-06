llr.bsreg <- function(target, dataset, threshold = 0.05, wei = NULL) {
  
  threshold <- log(threshold)
  dm <- dim(dataset)
  if ( is.null(dm) ) {
    n <- length(dataset)
    p <- 1
  } else {
    n <- dm[1]  ## sample size 
    p <- dm[2]  ## number of variables
  }  
  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. 
                 No backward procedure was attempted")
  } else {
    #check for NA values in the dataset and replace them with the variable median or the mode
    if ( any(is.na(dataset)) ) {
      #dataset = as.matrix(dataset);
      warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
      if ( is.matrix(dataset) )  {
        dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
      } else {
        poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
        for( i in poia )  {
          xi <- dataset[, i]
          if( is.numeric(xi) )   {                    
            xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
          } else if ( is.factor( xi ) ) {
            xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
          }
          dataset[, i] <- xi
        }
      }
    }
    runtime <- proc.time()
    dataset <- as.data.frame(dataset)
    
    if (p == 1) {
      ini <- survival::survreg( target ~.,  data = dataset, weights = wei, dist = "loglogistic" )
      mod <- anova(ini)
      stat <- mod[2, 2]
      dof <- mod[2, 1]
    } else {
      ini <- survival::survreg( target ~.,  data = dataset, weights = wei, dist = "loglogistic" )
      dofini <- length( ini$coefficients )
      stat <- dof <- numeric(p)
      for (i in 1:p) {
        mod <- survival::survreg( target ~.,  data = dataset[, -i ,drop = FALSE], weights = wei, dist = "loglogistic" )
        stat[i] <- 2 * abs(logLik(mod) - logLik(ini) )
        dof[i] <- dofini - length( mod$coefficients ) 
      }
    }
    mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
    colnames(mat) <- c("variable", "log.p-values", "statistic" )
    rownames(mat) <- 1:p 
    
    sel <- which.max( mat[, 2] )
    info <- matrix( c(0, -10, -10) , ncol = 3 )
    
    if ( mat[sel, 2] < threshold ) {
      runtime <- proc.time() - runtime
      res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = "censIndLLR", final = ini ) 
      
    } else {
      
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE] 
      dat <- dataset[, -sel ,drop = FALSE]
    
    i <- 1  
    
    while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   
      i <- i + 1
      k <- p - i + 1
      ini <- survival::survreg( target ~., data = dat, weights = wei, dist = "loglogistic" )
      
      if ( k == 1 ) {
        mod <- survival::survreg(target ~ 1, data = dat, weights = wei, dist = "loglogistic")
        stat <- 2 * abs( logLik(ini) - logLik(mod) )
        dof <- length( ini$coefficients ) - length( mod$coefficients ) 
        pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
        
        if (pval > threshold ) {
          final <- "No variables were selected"
          info <- rbind(info, c(mat[, 1], pval, stat) )
          dat <- dataset[, -info[, 1], drop = FALSE ]
          mat <- matrix(nrow = 0, ncol = 3)
        } else {
          info <- rbind(info, c(0, -10, -10)) 
          final <- ini
          mat[, 2:3] <- c(pval, stat)
        }
        
      } else {
        
        stat <- dof <- numeric(k)
        
        for (j in 1:k) {
          mod <- survival::survreg( target ~.,  data = dat[, -j, drop = FALSE], weights = wei, dist = "loglogistic" )
          stat[j] <- 2 * abs( logLik(mod) - logLik(ini) )
          dof[j] <- length( ini$coefficients ) - length( mod$coefficients ) 
        }
        mat[, 2:3] <- cbind( pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
        sel <- which.max( mat[, 2] )
        
        if ( mat[sel, 2] < threshold ) {
          final <- ini
          info <- rbind(info, c(0, -10, -10) )
          
        } else {
          info <- rbind(info, mat[sel, ] )
          mat <- mat[-sel, ,drop = FALSE] 
          dat <- dataset[, -info[, 1], drop = FALSE ] 
        }
      }  
    }  ## end while

    info <- info[ info[, 1] > 0, , drop = FALSE ]
    res <- list(runtime = runtime, info = info, mat = mat, ci_test = "censIndLLR", final = final ) 
    
    }  ## end  if ( mat[sel, 2] < threshold ) else 
  }  ##  end if (p > n)  else 

}    






