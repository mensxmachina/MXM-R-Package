bic.betafsreg <- function( target, dataset, wei = NULL, tol = 2, ncores = 1 ) {
  
  p <- ncol(dataset)  ## number of variables
  moda <- list()
  k <- 1   ## counter
  n <- length(target)  ## sample size
  con <- log(n)
  tool <- NULL
  info <- matrix( 0, ncol = 2 )
  sela <- NULL
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any( is.na(dataset) ) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
      for( i in poia )  {
        xi <- dataset[, i]
        if ( is.numeric(xi) ) {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) ) {
          xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        }
        dataset[, i] <- xi
      }
    }
  }
  ##################################
  # target checking and initialize #
  ##################################
  runtime <- proc.time()

    if ( is.null(wei) ) {
      ini <-  - 2 * Rfast::beta.mle(target)$loglik + 2 * con 
    } else ini <-  - 2 * betamle.wei(target, wei)$loglik + 2 * con 
    bico <- beta.regs(target, dataset, wei, logged = TRUE, ncores = ncores)[, 3]
    mat <- cbind(1:p, bico)
    bico <- NULL
    colnames(mat) <- c("variable", "BIC")
    rownames(mat) <- 1:p
    sel <- which.min( mat[, 2] )
    
    if ( ini - mat[sel, 2] > tol ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE]
      sela <- sel
      mi <- beta.reg( target, dataset[, sel], wei = wei )
      tool[1] <-  - 2 * mi$loglik + ( length(mi$be) + 1 ) * con
      moda[[ 1 ]] <- mi
    }  else  {
      info <- info  
      sela <- NULL
    }
    ######
    ###     k equals 2
    ######
    if ( length(moda) > 0  &  nrow(mat) > 0 ) {
      k <- 2
      pn <- p - k  + 1
      mod <- list()
      
      if ( ncores <= 1 ) {
        bico <- numeric( pn )
        for ( i in 1:pn ) {
          ma <- beta.reg( target, dataset[, c(sel, mat[i, 1]) ], wei = wei )
          bico[i] <-  - 2 * ma$loglik + ( length(ma$be) + 1 ) * con
        }
        mat[, 2] <- bico 
      } else {
        cl <- makePSOCKcluster(ncores)
        doParallel::registerDoParallel(cl)
        mod <- foreach( i = 1:pn, .combine = rbind, .export = "beta.reg") %dopar% {
          ww <- beta.reg( target, dataset[, c(sel, mat[i, 1]) ], wei = wei )
          return( - 2 * ww$loglik + ( length(ma$be) + 1 ) * con )
        }
        stopCluster(cl)
        mat[, 2] <- mod
      }
      ina <- which.min( mat[, 2] )
      sel <- mat[ina, 1]
      if ( tool[1] - mat[ina, 2] <= tol ) {
        info <- info
        sela <- NULL
      } else {
        tool[2] <- mat[ina, 2]
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina, , drop = FALSE]
        mi <- beta.reg( target, dataset[, sela], wei = wei )
        tool[2] <-  - 2 * mi$loglik + ( length(mi$be) + 1 ) * con
        moda[[ 2 ]] <- mi
      }
    }
    #########
    ####      k is greater than 2
    #########
    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while ( k < n - 15 & tool[ k - 1 ] - tool[ k ] > tol & nrow(mat) > 0 ) {  
        k <- k + 1
        pn <- p - k + 1
        if (ncores <= 1) {
          for ( i in 1:pn ) {
            ma <- beta.reg( target, dataset[, c(sela, mat[i, 1]) ], wei = wei )
            mat[i, 2] <-  - 2 * ma$loglik + ( length(ma$be) + 1 ) * con
          } 
        } else {
          cl <- makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach( i = 1:pn, .combine = rbind, .export = "beta.reg") %dopar% {
            ww <- beta.reg( target, dataset[, c(sela, mat[i, 1]) ], wei = wei )
            return(  - 2 * ww$loglik + ( length(ww$be) + 1 ) * con )
          }
          stopCluster(cl)
          mat[, 2] <- mod
        }
        ina <- which.min( mat[, 2] )
        sel <- mat[ina, 1]
        if ( tool[k - 1] - mat[ina, 2]  <= tol ) {
          info <- rbind( info,  c( -10, Inf ) )
          tool[k] <- Inf
        } else {
          tool[k] <- mat[ina, 2]
          info <- rbind(info, mat[ina, ] )
          sela <- info[, 1]
          mat <- mat[-ina, , drop = FALSE]
          ma <- beta.reg( target, dataset[, sela], wei = wei )
          tool[k] <-  - 2 * ma$loglik + ( length(ma$be) + 1 ) * con         
          moda[[ k ]] <- ma
        }
      }
    }
    runtime <- proc.time() - runtime

    d <- length(sela)
    final <- NULL
    if ( d >= 1 ) {
      final <- beta.reg( target, dataset[, sela], wei = wei )
      info <- info[1:d, , drop = FALSE]
      colnames(info) <- c( "variables", "BIC" )
      rownames(info) <- info[, 1]
    }
    list(runtime = runtime, mat = t(mat), info = info, ci_test = "testIndBeta", final = final )
  } 
  
