iamb.bs <- function(target, dataset, threshold = 0.05, wei = NULL, test = NULL, user_test = NULL) {
  
  threshold <- log(threshold)
  dm <- dim(dataset)
  n <- dm[1]  ## sample size 
  p <- dm[2]  ## number of variables
  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. No backward procedure was attempted")
    
  } else {
    
    runtime <- proc.time()
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
          if ( is.numeric(xi) )  {                    
            xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
          } else if ( is.factor( xi ) ) {
            xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
          }
          dataset[, i] <- xi
        }
      }
    }
    dataset <- as.data.frame(dataset)
    la <- length( unique(target) )
    ##################################
    # target checking and initialize #
    ################################## 
    if ( is.null(test)  &  is.null(user_test) ) {
      ## surival data
      if ( sum( class(target) == "Surv" ) == 1 ) {
        ci_test <- test <- "censIndCR"
        ## ordinal, multinomial or perhaps binary data
      } else if ( is.factor(target) ||  is.ordered(target) || la== 2 ) {
        ci_test <- test <- "testIndLogistic"
        ## count data
      } else if ( la > 2  &  !is.factor(target) ) {
        if ( sum( round(target) - target ) == 0 ) {
          ci_test <- test <- "testIndPois"
        } else  ci_test <- test <- "testIndReg"  
      }
    }
    
    av_models <- c("testIndReg", "testIndBeta", "censIndCR", "testIndRQ", "censIndWR", "testIndLogistic", "testIndPois", 
                   "testIndNB", "testIndZIP", "testIndGamma", "testIndNormLog", "testIndTobit", "testIndMMReg", 
				   "testIndMultinom", "testIndOrdinal");
    
    ci_test <- test
    test <- match.arg(test, av_models, TRUE);
    ############ 
    ###  GLMs 
    ############
    
    if ( test == "testIndPois"  ||  test == "testIndReg"  ||  ( test == "testIndLogistic"  &  la == 2 ) ) {
      res <- iamb.glmbs(target = target, dataset = dataset, wei = wei, threshold = exp(threshold) ) 
      
    } else if ( test == "testIndBeta" ) {
      res <- iamb.betabs(target = target, dataset = dataset, wei = wei, threshold = exp(threshold) ) 
      
    } else if ( test == "testIndZIP" ) {
      res <- iamb.zipbs(target = target, dataset = dataset, wei = wei, threshold = exp(threshold) ) 	
      
    } else if ( test == "testIndGamma" ) {
      res <- iamb.gammabs(target = target, dataset = dataset, wei = wei, threshold = exp(threshold) ) 	
      
    } else if ( test == "testIndNormLog" ) {
      res <- iamb.normlogbs(target = target, dataset = dataset, wei = wei, threshold = exp(threshold) ) 	
      
    } else if ( test == "testIndTobit" ) {
      res <- iamb.tobitbs(target = target, dataset = dataset, wei = wei, threshold = exp(threshold) ) 	
      
    } else {
      
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

      } else if ( !is.null(user_test) ) {
        test <- user_test
      }
      ###################
      ###################
      a1 <- internaliamb.bs( target = target, dataset = dataset, threshold = threshold, test = test, wei = wei, p = p ) 
      ind <- 1:p
      a2 <- list()
      poies <- a1$mat[, 1]
      if ( length(poies) > 0 ) {
        ind[-poies] <- 0
        ind <- ind[ind > 0]
        dat <- dataset[, poies, drop = FALSE ]
        a2 <- internaliamb.bs(target = target, dataset = dat, threshold = threshold, test = test, wei = wei, p = length(ind))  
        poies <- a2$mat[, 1]
        ind[-poies] <- 0
        ind <- ind[ind > 0]
        dat <- dat[, poies, drop = FALSE]
        i <- 2
      } else {
        ind <- NULL
        a2$mat <- NULL  
        a2$final <- paste("No variables were selected")
      } 
      while ( length(a1$mat[, 1]) - length(a2$mat[, 1]) != 0 ) {
        i <- i + 1
        a1 <- a2
        a2 <- internaliamb.bs( target = target, dataset = dat, threshold = threshold, test = test, wei = wei, p = length(ind) ) 
        poies <- a2$mat[, 1]
        if ( length(poies) > 0 ) {
          ind[-poies] <- 0
          ind <- ind[ind > 0]
          dat <- dat[, poies, drop = FALSE]
        } else  {
          dat <- NULL  
          ind <- NULL
        }  
      }
      
      runtime <- proc.time() - runtime
    # if ( !is.null(a2$mat) ) {
    #   final <- test( target ~.,  data = as.data.frame(dataset[, ind]), weights = wei )
    #   a2$mat[, 1] <- ind 
    # } else final <- "No variables were selected"
    res <- list(runtime = runtime, ci_test = ci_test, vars = ind, mat = a2$mat, final = a2$final ) 
  
    }
  }  
  res
}
      


      
      
         
