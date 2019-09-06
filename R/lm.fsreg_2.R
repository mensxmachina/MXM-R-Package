lm.fsreg_2 <- function(target, dataset, iniset = NULL, threshold = 0.05, wei = NULL, stopping = "BIC", tol = 2, ncores = 1 ) {
  
  threshold = log(threshold)
  p <- dim(dataset)[2]  ## number of variables
  pval <- stat <- dof <- numeric( p )  
  moda <- list()
  ## percentages
  if ( min( target ) > 0  &  max( target ) < 1 )   target <- log( target / (1 - target) ) 
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any(is.na(dataset)) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
      for ( i in poia )  {
        xi <- dataset[, i]
        if ( is.numeric(xi) )
        {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) )   xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        dataset[, i] <- xi
      }
    }
  }
  
  ## is there already an initial set of variables to start with?
  pa <- NCOL(iniset)
  da <- 1:pa
  dataset <- cbind(iniset, dataset)
  dataset <- as.data.frame(dataset)
  
  n <- length(target)  ## sample size
  tool <- numeric( length( min(n, p) ) )
  info <- matrix( c( 1e300, 0, 0 ), ncol = 3 )
  k <- 1   ## counter k is 1, step 1
  mi <- lm(target ~. , data = as.data.frame(iniset), weights = wei, y = FALSE, model = FALSE)
  runtime <- proc.time()
  
  if (ncores <= 1) {
    
      for (i in 1:p) {
        ww <- lm( target ~., data = dataset[, c(da, pa + i)], weights = wei, y = FALSE, model = FALSE )
        tab <- anova(mi, ww)
        stat[i] <- tab[2, 5] 
        df1 <- tab[2, 3]   ;  df2 = tab[2, 1]
        pval[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
      } 
      mat <- cbind(1:p, pval, stat)
   
  } else {
    
    cl <- makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach( i = 1:p, .combine = rbind ) %dopar% {
       ww <- lm( target ~., data = dataset[, c(da, pa + i)], weights = wei, y = FALSE, model = FALSE )
       tab <- anova( mi, ww )
       stat <- tab[2, 5] 
       df1 <- tab[2, 3]   ;  df2 = tab[2, 1]
       pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE )
       return( c(pval, stat) )
    }
    stopCluster(cl)
    mat <- cbind(1:p, mod)      
  }
  
  colnames(mat) <- c( "variables", "log.p-value", "stat" )
  rownames(mat) <- 1:p
  sel <- which.min(mat[, 2])
  sela <- pa + sel
  
  if ( mat[sel, 2] < threshold ) {
    
    info[k, ] <- mat[sel, ]
    mat <- mat[-sel, , drop = FALSE] 

    if ( stopping == "adjrsq" ) {
      mi = lm( target ~., data = dataset[, c(da, sel) ], weights = wei, y = FALSE, model = FALSE )
      tool[k] <- as.numeric( summary( mi )[[ 9 ]] )
    } else if ( stopping  == "BIC" ) { 
      mi = lm( target ~., data = dataset[, c(da, sel) ], y = FALSE, model = FALSE )
      tool[k] <- BIC( mi )
    }
    moda[[ k ]] <- mi
    
  }
  ######
  ###### k equal to 2
  ######
  if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
    
    k <- k + 1
    pn <- p - k + 1   
    
    if ( ncores <= 1 ) {
      for (i in 1:pn) {
        ww = lm( target ~., data = dataset[, c(da, sel, pa + mat[i, 1]) ], weights = wei, y = FALSE, model = FALSE )
        tab = anova( mi, ww )
        mat[i, 3] = tab[2, 5] 
        df1 = tab[2, 3]   ;  df2 = tab[2, 1]
        mat[i, 2] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
      }      
    } else {
      cl <- makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
        ww <- lm( target ~., data = dataset[, c(da, sel, pa + mat[i, 1]) ], weights = wei, y = FALSE, model = FALSE )
        tab <- anova( mi, ww )
        stat <- tab[2, 4] 
        df1 <- tab[2, 3]   ;  df2 = tab[2, 1]
        pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE )
        return( c(pval, stat) )
      }
      stopCluster(cl)
      mat <- cbind(mat[, 1], mod)   
      
    }
  }
  
  ina <- which.min(mat[, 2])
  sel <- pa + mat[ina, 1]     
  
  if ( stopping == "adjrsq" ) {
  
    if ( mat[ina, 2] < threshold ) {
      ma <- lm( target ~., data = dataset[, c(da, sela, sel) ], weights = wei, y = FALSE, model = FALSE )
      tool[k] <- as.numeric( summary( ma )[[ 9 ]] )
      if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
        info <- info       
      } else {  
        info <- rbind(info, mat[ina, ] )
        sela <- c(sela, sel)
        mat <- mat[-ina, , drop = FALSE] 
        moda[[ k ]] <- ma
      }
      
    } else  info <- info
    
  } else if ( stopping == "BIC" ) {
  
    if ( mat[ina, 2] < threshold ) {
       ma <- lm( target ~., data = dataset[, c(da, sela, sel) ], weights = wei, y = FALSE, model = FALSE )
       tool[k] <- BIC( ma )
       if ( tool[ k - 1] - tool[ k ] <= tol ) {
         info <- info
      } else {  
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina, , drop = FALSE]
        moda[[ k ]] <- ma
      }
    } else  info <- info
    
  }
  ###########
  ######   k greater than 2
  ###########
  if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
    
    while ( info[k, 2] < threshold &  k < n - 15  &  abs( tool[ k ] - tool[ k - 1 ] ) > tol  &  nrow(mat) > 0 )  {
      
      k <- k + 1   
      pn <- p - k + 1 
      
      if ( ncores <= 1 ) {
        for ( i in 1:pn ) {
          ww <- lm( target ~., data = dataset[, c(da, sela, pa + mat[i, 1]) ], weights = wei, y = FALSE, model = FALSE )
          tab <- anova( ww )
          mat[i, 3] <- tab[pa + k, 4] 
          df1 <- tab[pa + k, 1]   ;  df2 = tab[pa + k + 1, 1]
          mat[i, 2] <- pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
        }
      } else {
        cl <- makePSOCKcluster(ncores)
        doParallel::registerDoParallel(cl)
        mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
          ww <- lm( target ~., data = dataset[, c(da, sela, pa + mat[i, 1]) ], weights = wei, y = FALSE, model = FALSE )
          tab <- anova(ma, ww)
          stat <- tab[2, 5] 
          df1 <- tab[2, 3]   ;  df2 = tab[2, 1]
          pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE )
          return( c( pval, stat ) )
        }
        stopCluster(cl)
        mat <- cbind( mat[, 1], mod )   
      }
      
      ina <- which.min(mat[, 2])
      sel <- pa + mat[ina, 1]   
      
      if ( stopping == "BIC" ) {
	  
        if ( mat[ina, 2] < threshold ) {	
          ma <- lm( target ~., data = dataset[, c(da, sela, sel)], weights = wei, y = FALSE, model = FALSE )
          tool[k] <- BIC( ma )
          
          if ( tool[ k - 1] - tool[ k ] <= tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- c(sela, sel)
            mat <- mat[-ina, , drop = FALSE] 
            moda[[ k ]] <- ma
          }
          
        } else  info <- rbind(info, c( 1e300, 0, 0 ) )
        
      } else if ( stopping == "adjrsq" ) {
	    if ( mat[ina, 2] < threshold ) {
          
          ma <- lm( target ~., data = dataset[, c(da, sela, sel)], weights = wei, y = FALSE, model = FALSE )
          tool[k] <- as.numeric( summary(ma)[[ 9 ]] )           
          if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )
            
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- info[, 1]
            mat <- mat[-ina, , drop = FALSE]
            moda[[ k ]] <- ma
          }
          
        } else   info <- rbind(info, c( 1e300, 0, 0 ) )
        
      }
    }
    
  }
  
  runtime <- proc.time() - runtime
  
  d <- length(moda)
  
  if (d == 0) {
    final <- lm( target ~., data = as.data.frame( iniset ), weights = wei, y = FALSE, model = FALSE )
  } else {
    final <- lm( target ~., data =  dataset[, c(da, sela) ], weights = wei, y = FALSE, model = FALSE )
    info <- info[1:d, , drop = FALSE ]
    info <- cbind( info, tool[ 1:d ] ) 
    colnames(info) <- c( "variables", "log.p-value", "stat", stopping )
    rownames(info) <- info[, 1]
  }
  
  list( runtime = runtime, mat = t(mat), info = info, ci_test = "testIndReg", final = final ) 
}
