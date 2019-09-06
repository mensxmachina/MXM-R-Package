bs.reg <- function(target, dataset, threshold = 0.05, wei = NULL, test = NULL, user_test = NULL) {
  
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
   #check for NA values in the dataset and replace them with the variable median or the mode
   if( any(is.na(dataset)) ) {
     warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
     if ( is.matrix(dataset) )  {
       dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
     } else {
       poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
       for ( i in poia )  {
         xi <- dataset[, i]
         if ( is.numeric(xi) )  {                    
           xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
         } else if ( is.factor( xi ) ) {
           xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
         }
         dataset[, i] <- xi
       }
     }
   }
  ## dependent (target) variable checking if no test was given, 
  ## but other arguments are given. For some cases, these are default cases
  if ( is.null(test)  &  is.null(user_test) ) {
    ## surival data
    if ( identical( class(target), "Surv" ) ) {
      ci_test <- test <- "censIndCR"
      ## ordinal, multinomial or perhaps binary data
    } else if ( is.factor(target)  &  length( unique(target) ) == 2 ) {
      ci_test <- test <- "testIndLogistic"
      ## count data
    } else if ( length( unique(target) ) > 2  &  !is.factor(target) ) {
      if ( sum( round(target) - target ) == 0 ) {
        ci_test <- test <- "testIndPois"
      } else  ci_test <- test <- "testIndReg"  
    }
  }
   
  av_models = c("testIndReg", "testIndMMReg", "testIndBeta", "censIndCR", "censIndWR", "censIndLLR", "testIndRQ",
                "testIndLogistic", "testIndPois", "testIndNB", "testIndZIP", "testIndBinom", "testIndGamma", 
                "testIndNormLog", "testIndTobit", "testIndClogit", "testIndFisher", "testIndQPois", 
                "testIndQBinom", "testIndMultinom", "testIndOrdinal", "testIndSPML")
  
  ci_test <- test
  test <- match.arg(test, av_models, TRUE);
   ############ 
   ###  GLMs 
   ############
    if ( test == "testIndPois"  |  test == "testIndReg"  | test == "testIndLogistic"  | test == "testIndBinom" | test == "testIndMMReg" ) {
     res <- glm.bsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold ) ) 
      
    } else  if ( test == "testIndFisher" ) {
      res <- cor.bsreg(target, dataset, threshold = exp( threshold ) ) 
        
    } else if ( test == "testIndBeta" ) {
      res <- beta.bsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold ) ) 
	  
    } else if ( test == "testIndZIP" ) {
      res <- zip.bsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold ) ) 	
      
    } else if ( test == "testIndGamma" ) {
      res <- gammabsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold )) 	
      
    } else if ( test == "testIndNormLog" ) {
      res <- normlog.bsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold ) ) 	

    } else if ( test == "testIndTobit" ) {
      res <- tobit.bsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold ) ) 	
      
    } else if ( test == "testIndClogit" ) {
      res <- clogit.bsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold ) ) 	
      
    } else if ( test == "testIndQPois" ) {
      res <- quasipois.bsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold ) ) 
      
    } else if ( test == "testIndQBinom" ) {
      res <- quasibinom.bsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold ) ) 
      
	  } else if ( test == "gSquare" ) {
	    res <- bs.g2(target, dataset, threshold = exp( threshold ) )
	    
	  } else if ( test == "testIndSPML" ) {
	    res <- spml.bsreg(target = target, dataset = dataset, threshold = exp( threshold ) ) 	
	    
	  } else if ( test == "censIndLLR" ) {
	    res <- llr.bsreg(target = target, dataset = dataset, threshold = exp( threshold ) ) 	
	    
    } else {
	   
     if ( test == "censIndCR" ) {
        test <- survival::coxph 

	  } else if ( test == "censIndWR" ) {
      test <- survival::survreg 
      
    } else if ( test == "testIndOrdinal" ) {
      test <- ordinal::clm

    } else if ( test == "testIndMultinom" ) {
      test <- nnet::multinom
	  
    } else if ( test == "testIndNB" ) {
      test <- MASS::glm.nb

    } else if ( test == "testIndRQ" ) {
      test <- quantreg::rq
	
	  } else if (test == "testIndRQ") {
	    test <- quantreg::rq
	  }
    ###################
    ###################
    dataset <- as.data.frame(dataset)
    ini <- test( target ~.,  data = dataset, weights = wei )
	  dofini <- length( coef(ini) )
	  likini <- logLik(ini) 
	  stat <- dof <- numeric(p)
	  if (p == 1) {
		   mod <- test( target ~ 1, weights = wei )
       if ( ci_test == "censIndCR")  {
		     stat <- 2 * likini - 2 * mod$loglik
		     dof <- dofini
		   } else  stat <- 2 * ( likini - logLik(mod) )
		   dof <- dofini - length( coef(mod) ) 	  
	  }	else {
	    for (j in 1:p) {
		    mod <- test( target ~.,  data = dataset[, -j, drop = FALSE], weights = wei )
		    stat[j] <- 2 * ( likini - logLik(mod) )
		    dof[j] <- dofini - length( coef(mod) ) 
	    }
	  }  
      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 
      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )

      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = ci_test, final = ini ) 
    
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE] 
        final <- "No variables were selected"
		
        i <- 1  
         
         while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   

           ini <- test( target ~., data = dat, weights = wei )
           likini <- logLik(ini) 
           dofini <- length( coef(ini) )
           i <- i + 1        
           k <- p - i + 1
		      
		      if ( k == 1 ) {
		        mod <- test(target ~ 1, weights = wei)
		        if ( ci_test == "censIndCR")  {
		          dof <- dofini
		          stat <- 2 * ( likini - mod$loglik )
		        } else {
		          stat <- 2 * ( likini - logLik(mod) )
		          dof <- dofini - length( coef(mod) ) 
		        }  
		        pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
		      
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
            stat <- dof <- numeric(k)
             
	          for (j in 1:k) {
		          mod <- test( target ~.,  data = dat[, -j, drop = FALSE], weights = wei )
		          stat[j] <- 2 * ( likini - logLik(mod) )
		          dof[j] <- dofini - length( coef(mod) ) 
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
 
           }  ##  end if ( k == 1 ) 
         }  ##  end while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 ) 
		 
        runtime <- proc.time() - tic
        info <- info[ info[, 1] > 0, , drop = FALSE]
        colnames(mat) <- c("Variables", "log.p-values", "statistic")
        res <- list(runtime = runtime, info = info, mat = mat, ci_test = ci_test, final = final ) 

      }  ##  end if ( mat[sel, 2] < threshold )  else 
    }  ## end if test ==
    
  }  ## end if ( p > n ) else
  
  res
}  
   
    