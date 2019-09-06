fbedreg.bic <- function(target, dataset, wei = NULL, fbedreg.object, test = NULL, graph = TRUE) {
  
  if ( fbedreg.object$info[1, 1] == 0 ) {
    res <- NULL
    mod <- paste("No associations were found, hence no model is produced.")
  } else {
  
  if ( any(is.na(dataset) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- which( is.na(dataset), arr.ind = TRUE )[2]
      for ( i in poia )  {
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
  
  ci_test <- test 
  p <- dim(fbedreg.object$res)[1]
  bic <- numeric(p)
  dataset <- as.data.frame(dataset[, fbedreg.object$res[, 1], drop = FALSE])

  if ( ci_test == "testIndFisher" | ci_test == "testIndReg" ) {
    
    for (i in 1:p) {
      mod <- lm( target ~ ., data = dataset[, 1:i, drop = FALSE ], weights = wei )
      bic[i] <- BIC(mod)      
    }
    
  } else if (test == "testIndMMReg") {
    for (i in 1:p) {
      mod <- MASS::rlm(target ~., data = dataset[, 1:i, drop = FALSE ], maxit = 2000, weights = wei )  
      bic[i] <- BIC( mod )    
    }
  
  } else if ( ci_test == "testIndRQ" ) {
    for (i in 1:p) {
      mod <- quantreg::rq( target ~., data = dataset[, 1:i, drop = FALSE ], weights = wei  )
      la <- logLik(mod)
      bic[i] <-  - 2 * as.numeric( la ) + attr(la, "df") * log( length(target) )   
    }
    
  } else if ( ci_test == "testIndBeta" ) {
    for (i in 1:p) {
      mod <- beta.mod( target, dataset[, 1:i, drop = FALSE ], wei = wei )
      bic[i] <-  - 2 * mod$loglik + ( length( mod$be ) + 1 ) * log( length(target) )
    }

  } else if ( ci_test == "testIndPois") {
    for (i in 1:p) {
      mod <- glm( target ~ ., data = dataset[, 1:i, drop = FALSE ], family = poisson, weights = wei )
      bic[i] <- BIC( mod )
    }  
    
  } else if ( ci_test == "testIndNB" ) {
    for (i in 1:p) {
      mod <- MASS::glm.nb( target ~ ., data = dataset[, 1:i, drop = FALSE ], weights = wei )
      bic[i] <- BIC(mod)
    }
    
  } else if ( ci_test == "testIndZIP" ) {
    for (i in 1:p) {
      mod <- zip.mod( target, dataset[, 1:i, drop = FALSE ], wei = wei )
      bic[i] <-  -2 * mod$loglik + ( length( coef(mod$be) ) + 1) * log( length(target) )
    }
    
  } else if (ci_test == "testIndIGreg") {
    for (i in 1:p) {
      mod <- glm(target ~., data = dataset[, 1:i, drop = FALSE ], family = inverse.gaussian(log), weights = wei)
      bic[i] <- BIC(mod)
    }
    
  } else if ( is.matrix(target)  & ci_test == "testIndMVreg" ) {
    if ( all(target > 0 & target < 1)  &  Rfast::Var( Rfast::rowsums(target) ) == 0 )   target = log( target[, -1]/(target[, 1]) ) 
    for (i in 1:p) {
      mod <- lm( target ~., data = dataset[, 1:i, drop = FALSE ], weights = wei )
      bic[i] <- NULL 
    }
    
  } else if ( is.matrix(target)   &  ci_test == "testIndBinom" ) {
    y <- target[, 1] /target[, 2]
    for (i in 1:p) {
      mod <- glm( y ~., data = dataset[, 1:i, drop = FALSE ], weights = target[, 2], family = binomial )
      bic[i] <- BIC(mod)
    }
    
  } else if ( ci_test == "testIndGamma" ) {
    for (i in 1:p) {
      mod <- glm( target ~ ., data = dataset[, 1:i, drop = FALSE ], weights = wei, family = Gamma(link = log) )
      bic[i] <- BIC(mod)
    }
    
  } else if ( ci_test == "testIndNormLog" ) {
    for (i in 1:p) {
      mod <- glm( target ~ ., data = dataset[, 1:i, drop = FALSE ], weights = wei, family = gaussian(link = log) )
      bic[i] <- BIC(mod)
    }
    
  } else if ( ci_test == "testIndTobit" ) {
    for (i in 1:p) {
      mod <- survival::survreg( target ~ ., data = dataset[, 1:i, drop = FALSE ], weights = wei, dist = "gaussian" )
      bic[i] <-  - 2 * as.numeric( logLik(mod) ) + ( length( mod$coefficients ) + 1 ) * log( NROW(dataset) )
    }
    
  } else if (ci_test == "censIndCR") {
    for (i in 1:p) {
      mod <- survival::coxph( target ~ ., data = dataset[, 1:i, drop = FALSE ], weights = wei )
      bic[i] <- BIC(mod)
    }
    
  } else if (ci_test == "censIndWR") {
    for (i in 1:p) {
      mod <- survival::survreg( target ~ ., data = dataset[, 1:i, drop = FALSE ], weights = wei )
      bic[i] <-  - 2 * as.numeric( logLik(mod) ) + ( length( mod$coefficients ) + 1 ) * log( NROW(dataset) )
    }
    
  } else if (ci_test == "testIndClogit") {
    case <- as.logical(target[, 1]);  
    id <- target[, 2]
    for (i in 1:p) {
      mod <- survival::clogit(case ~ . + strata(id), data = dataset[, 1:i, drop = FALSE ] )  ## wieghts are ignored here anyway
      bic[i] <- BIC(mod)
    }
    
  } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
    
    if ( length( unique(target)) == 2 ) {
      for (i in 1:p) {
        mod <- glm( target ~., data = dataset[, 1:i, drop = FALSE ], family = binomial, weights = wei ) 
        bic[i] <- BIC(mod)
      }
    } else if ( !is.ordered(target) ) { 
      target <- as.factor( as.numeric( as.vector(target) ) )
      for (i in 1:p) {
        mod <- nnet::multinom( target ~., data = dataset[, 1:i, drop = FALSE ], trace = FALSE, weights = wei )
        bic[i] <- BIC(mod)
      }
    } else if ( is.ordered(target) ) {
      for (i in 1:p) {
        mod <- ordinal::clm( target ~., data = dataset[, 1:i, drop = FALSE ], weights = wei  )
        bic[i] <- BIC(mod)
      }  
    }    
    
  }
  res <- cbind(fbedreg.object$res, bic)
  colnames(res)[dim(res)[2]] <- "BIC"
  if  (graph)   plot(1:p, bic, type = "b", xlab = "Number of selected variables", ylab = "BIC values")
  
  }  ## end if ( fbedreg.object$info[1, 1] == 0 )
  
  list(res = res, mod = mod) 
  
}

