bic.fsreg <- function( target, dataset, test = NULL, wei = NULL, tol = 2, ncores = 1 ) {

  p <- ncol(dataset)  ## number of variables
  bico <- numeric( p )
  moda <- list()
  k <- 1   ## counter
  n <- length(target)  ## sample size
  con <- log(n)
  tool <- NULL
  info <- matrix( 0, ncol = 2 )
  result <- NULL
  sela <- NULL
  #check for NA values in the dataset and replace them with the variable median or the mode
  if( any( is.na(dataset) ) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    }else{
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
  if ( is.null(test) ) {
    ## linear regression 
    if ( is.numeric(target)  ||  is.vector(target)  ) { 
	  la <- length( unique(target) )
	  if ( la > 2 ) {
	    if ( sum( round(target) - target ) == 0 )   test <- "testIndPois" 
	   } else if ( la == 2 ) {
	     test <- "testIndLogistic" 
	   }	 
       ## surival data
    } else if ( sum( class(target) == "Surv" ) == 1 ) {
      test <- "censIndCR"
    }  ## end if ( is.null(test) )
    #available conditional independence tests
  }	 ## end if ( is.null(test) )
    av_models = c("testIndReg", "testIndRQ", "testIndBeta", "censIndCR", "censIndWR", "censIndLLR", 
                  "testIndLogistic", "testIndPois", "testIndNB", "testIndZIP", "testIndGamma", 
                  "testIndNormLog", "testIndTobit", "testIndMMReg"); 
  
  dataset <- as.data.frame(dataset)
    
  if ( test == "testIndLogistic"  |  test == "testIndPois"  ||  test == "testIndReg" ) {
    result <- bic.glm.fsreg( target, dataset, wei = wei, tol = tol, ncores = ncores ) 
  
  } else if ( test == "testIndBeta" ) {
    result <- bic.betafsreg(target, dataset, wei = wei, tol = tol, ncores = ncores )

  } else if ( test == "testIndMMReg" ) {
    result <- bic.mm.fsreg(target, dataset, wei = wei, tol = tol, ncores = ncores )
	
  } else if ( test == "testIndZip" ) {
    result <- bic.zipfsreg(target, dataset, wei = wei, tol = tol, ncores = ncores )
    
  } else if ( test == "testIndGamma" ) {
    result <- bic.gammafsreg(target, dataset, wei = wei, tol = tol, ncores = ncores )
    
  } else if ( test == "testIndNormLog" ) {
    result <- bic.normlog.fsreg(target, dataset, wei = wei, tol = tol, ncores = ncores )

  } else if ( test == "testIndTobit" ) {
    result <- bic.tobit.fsreg(target, dataset, wei = wei, tol = tol, ncores = ncores )
	
  } else if ( test == "censIndWR" ) {
    result <- bic.wr.fsreg(target, dataset, wei = wei, tol = tol, ncores = ncores )
    
  } else if ( test == "censIndLLR" ) {
    result <- bic.llr.fsreg(target, dataset, wei = wei, tol = tol, ncores = ncores )
    
  } else if ( test == "testIndClogit" ) {
    result <- bic.clogit.fsreg(target, dataset, wei = wei, tol = tol, ncores = ncores )
    
  } else {
 
    ci_test <- test <- match.arg(test, av_models ,TRUE);
    #convert to closure type
    if ( test == "censIndCR" ) {
      test <- survival::coxph 

	  } else if ( test == "censIndWR" ) {
      test <- survival::survreg 
      
    } else if ( test == "testIndOrdinal" ) {  
      test <- ordinal::clm
  
    } else if (test == "testIndMultinom") {
      test <- nnet::multinom 
      
    } else if ( test == "testIndNB" ) {
      test <- MASS::glm.nb

    } else if ( test == "testIndRQ" ) {
      test <- quantreg::rq 
    } 
    runtime <- proc.time()
      
    ini <- test( target ~ 1 )
    if ( ci_test == "censIndCR" ) {
      ini <-  - 2 * ini$loglik
    } else {  
      la <- logLik(ini)
      ini <-  - 2 * as.numeric( la ) + attr(la, "df") * con   ## initial BIC  
    }   
    if (ncores <= 1) {
        for (i in 1:p) {
          mi <- test( target ~ dataset[, i], weights = wei )
		      la <- logLik(mi)
          bico[i] <-  - 2 * as.numeric( la ) + attr(la, "df") * con
        }
      mat <- cbind(1:p, bico)
      if( any( is.na(mat) ) )     mat[ which( is.na(mat) ) ] <- ini

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
        ww <- test( target ~ dataset[, i], weights = wei )
		    la <- logLik(ww)
        return( - 2 * as.numeric( la ) + attr(la, "df") * con )
      }
      stopCluster(cl)
      mat <- cbind(1:p, mod)
      if ( any( is.na(mat) ) )  mat[ which( is.na(mat) ) ] <- ini 
    }
    colnames(mat) <- c("variable", "BIC")
    rownames(mat) <- 1:p
    sel <- which.min( mat[, 2] )
    
    if ( ini - mat[sel, 2] > tol ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE]
      sela <- sel
      mi <- test( target ~ dataset[, sel], weights = wei )
      la <- logLik(mi)
      tool[1] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
      moda[[ 1 ]] <- mi
    } else  {
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
          ma <- test( target ~., data = dataset[, c(sel, mat[i, 1]) ], weights = wei )
		      la <- logLik(ma)
          bico[i] <-  - 2 * as.numeric( la ) + attr(la, "df") * con
        }
        mat[, 2] <- bico

      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
          ww <- test( target ~., data = dataset[, c(sel, mat[i, 1]) ], weights = wei )
          la <- logLik(ww)
          return( - 2 * as.numeric( la ) +  attr(la, "df") * con )
        }
        stopCluster(cl)
        mat[, 2] <- mod
      }
      ina <- which.min( mat[, 2] )
      sel <- mat[ina, 1]

      if ( tool[1] - mat[ina, 2] <= tol ) {
        info <- info
      } else {
        tool[2] <- mat[ina, 2]
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina, , drop = FALSE]
        mi <- test( target ~., data = dataset[, sela], weights = wei )
	      la <- logLik(mi)
        tool[2] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
        moda[[ 2 ]] <- mi
      }
   }
    #########
    ####      k is greater than 2
    #########
    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while ( k < n - 15 & tool[ k - 1 ] - tool[ k ] > tol &  nrow(mat) > 0 ) {

        k <- k + 1
        pn <- p - k + 1
        if (ncores <= 1) {
          for ( i in 1:pn ) {
            ma <- test( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei )
		        la <- logLik(ma)
            mat[i, 2] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
          }

        } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- test( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei )
		      	la <- logLik(ww)
            return( - 2 * as.numeric( la ) +  attr(la, "df") * con )
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
          ma <- test( target ~., data = dataset[, sela], weights = wei )
		      la <- logLik(ma)
          tool[k] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
          moda[[ k ]] <- ma
        }
      }
    }

    runtime <- proc.time() - runtime

    d <- length(sela)
    final <- NULL
    if ( d >= 1 ) {
      final <- test( target ~., data = dataset[, sela, drop = FALSE], weights = wei )
      info <- info[1:d, , drop = FALSE]
      colnames(info) <- c( "variables", "BIC" )
      rownames(info) <- info[, 1]
    }
    result = list( runtime = runtime, mat = t(mat), info = info, final = final, ci_test = ci_test)
  }  
  result
}
