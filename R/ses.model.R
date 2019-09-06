# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of signatures and generated models. It could be numeric from 1 to total number of signatures or "all" for all the 
## signatures. Default is 1.
ses.model <- function(target, dataset, wei = NULL, sesObject, nsignat = 1, test = NULL) {
  
  signature <- sesObject@selectedVars  
  if ( sum( is.na(signature) ) > 0 | length(signature) == 0 ) {
    mod <- paste("No associations were found, hence no model is produced.")
    signature = NULL
    res <- list(mod = mod, signature = signature)  
    
  } else {    
  
  if ( any(is.na(dataset) ) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if (class(dataset) == "matrix")  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- which( is.na(dataset), arr.ind = TRUE )[2]
      for( i in poia )  {
        xi <- dataset[, i]
        if(class(xi) == "numeric") {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) ) {
          xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        }
        dataset[, i] <- xi
      }
    }
  }
  
  if ( is.null(test) ) {  
    ci_test <- sesObject@test
  } else ci_test <- test 

  if ( nsignat == 1 || ( nsignat > 1 & nrow(sesObject@signatures) == 1 ) ) {
    p <- length(signature)
    
   if ( ci_test == "testIndFisher"  |  ci_test == "testIndReg" ) {
      mod = lm( target ~ ., data = data.frame(dataset[, signature ]), weights = wei )
      bic = BIC(mod)     
 
   } else if ( ci_test == "testIndMMReg" ) {
      mod = MASS::rlm(target ~., data = data.frame(dataset[, signature ]), maxit = 2000, weights = wei )
      bic = BIC( mod )  
     
   } else if ( ci_test == "testIndSpearman"  |  ci_test == "testIndRQ" ) {
      mod <- quantreg::rq( target ~., data = data.frame(dataset[, signature ]), weights = wei )
	    la <- logLik(mod)
      bic <-  - 2 * as.numeric( la ) + attr(la, "df") * log( length(target) )
      
    } else if ( ci_test == "testIndBeta" ) {     
      mod <- beta.mod( target, dataset[, signature ], wei = wei )
      bic <- - 2 * mod$loglik + ( length( mod$be ) + 1 ) * log( length(target) )

    } else if ( ci_test == "testIndPois") {
        mod <- glm( target ~ ., data = data.frame(dataset[, signature ]) , family = poisson, weights = wei )
        bic <- BIC( mod )
        
    } else if ( ci_test == "testIndQPois") {
      mod <- glm( target ~ ., data = data.frame(dataset[, signature ]), family = quasipoisson, weights = wei )
      bic <- NA

    } else if ( ci_test == "testIndQBinom") {
      mod <- glm( target ~ ., data = data.frame(dataset[, signature ]), family = quasibinomial, weights = wei )
      bic <- NA
      
    } else if ( ci_test == "testIndNB" ) {
      mod = MASS::glm.nb( target ~ ., data = data.frame(dataset[, signature ]), weights = wei )
      bic = BIC(mod)
      
    } else if ( ci_test == "testIndZIP" ) {
      mod <- zip.mod( target, dataset[, signature], wei = wei )
      bic <-  -2 * mod$loglik + ( length( coef(mod$be) ) + 1) * log( length(target) )
      ## bic = BIC(mod)
      
    } else if ( is.matrix(target) & ci_test == "testIndMVreg" ) {
      if ( all(target > 0 & target < 1)  &  Rfast::Var( Rfast::rowsums(target) ) == 0 )   target = log( target[, -1]/(target[, 1]) ) 
       mod = lm( target ~., data = data.frame(dataset[, signature ]), weights = wei )
       bic = NULL
      
    } else if ( is.matrix(target) &  ci_test == "testIndBinom" ) {
      mod = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, signature ]), weights = target[, 2], family = binomial )
      bic = BIC(mod)
      
    } else if ( ci_test == "testIndGamma" ) {
      mod <- glm( target ~ ., data = as.data.frame(dataset[, signature ]), weights = wei, family = Gamma(link = log) )
      bic <- BIC(mod)
      
    } else if ( ci_test == "testIndNormLog" ) {
      mod <- glm( target ~ ., data = as.data.frame(dataset[, signature ]), weights = wei, family = gaussian(link = log) )
      bic <- BIC(mod)
      
    } else if ( ci_test == "testIndTobit" ) {
      mod <- survival::survreg( target ~ ., data = as.data.frame(dataset[, signature ]), weights = wei, dist = "gaussian" )
      bic <-  - 2 * logLik(mod) + (length(mod$coefficients) + 1) * log( NROW(dataset)  )
      
    } else if (ci_test == "censIndCR") {
      mod = survival::coxph( target ~ ., data = data.frame(dataset[, signature ]), weights = wei )
      bic = BIC(mod)
      
    } else if (ci_test == "censIndWR") {
      mod = survival::survreg( target ~ ., data = as.data.frame(dataset[, signature ]), weights = wei )
      bic =  - 2 * logLik(mod) + (length(mod$coefficients) + 1) * log( NROW(dataset)  )

    } else if (ci_test == "censIndER") {
      mod = survival::survreg( target ~ ., data = data.frame(dataset[, signature ]), weights = wei, dist = "exponential" )
      bic = - 2 * logLik(mod) + (length(mod$coefficients) + 1) * log( NROW(dataset)  )
      
    } else if (ci_test == "censIndLLR") {
      mod <- survival::survreg( target ~ ., data = data.frame(dataset[, signature ]), weights = wei, dist = "loglogistic" )
      bic <-  - 2 * as.numeric( logLik(mod) ) + ( length( mod$coefficients ) + 1 ) * log( NROW(dataset) )
      
    } else if (ci_test == "testIndClogit") {
      case = as.logical(target[, 1]);  
      id = target[, 2]
      mod = survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , signature] ) )  ## weights are ignored here anyway
      bic = BIC(mod)
      
    } else if ( ci_test == "testIndLogistic" ) {
      mod = glm( target ~., data = data.frame(dataset[, signature ]) , family = binomial, weights = wei ) 
      bic = BIC(mod)

    } else if ( ci_test == "testIndMultinom" || ci_test == "gSquare" ) {
      mod = nnet::multinom( target ~., data = data.frame(dataset[, signature ]) , trace = FALSE, weights = wei )
      bic = BIC(mod)
        
    } else if ( ci_test == "testIndOrdinal" ) {
        mod = ordinal::clm( target ~., data = data.frame(dataset[, signature ]), weights = wei )
        bic = BIC(mod)
    }
    
    if ( is.null( colnames(dataset) ) ) {
      names(signature) = paste("Var", signature, sep = " ")
    } else    names(signature) = colnames(dataset)[signature]
    
    signature = c(signature, bic)
    names(signature)[length(signature)] = "bic"
  } 

  #############  more than one signatures
 
   if ( nsignat > 1 & nrow(sesObject@signatures) > 1 ) {
      
     if ( nsignat > nrow(sesObject@signatures) )  nsignat = nrow(sesObject@signatures)

     con <- log( NROW(dataset) )
     bic <- numeric(nsignat)
     signature <- sesObject@signatures[1:nsignat, , drop = FALSE] 
     mod <- list()
    
    for ( i in 1:nsignat ) {
	
     if ( ci_test == "testIndFisher" | ci_test == "testIndReg" ) {
        mod[[ i ]] = lm( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
        bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndMMReg" ) {
        mod[[ i ]] = MASS::rlm(target~., data = data.frame( dataset[, signature[i, ] ] ), maxit = 2000, weights = wei )
        bic[i] = BIC( mod[[ i ]] ) 
        
     } else if ( ci_test == "testIndSpearman"  |  ci_test == "testIndRQ" ) {
       mod[[ i ]] = quantreg::rq( target ~., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
	     la <- logLik( mod[[ i ]] ) 
       bic[i] <-  - 2 * as.numeric( la ) + attr(la, "df") * con
	  
     } else if ( ci_test == "testIndBeta" ) {
       mod[[ i ]] <- beta.mod( target, dataset[, signature[i, ] ], wei = wei )
       bic[i] <-  - 2 * mod[[ i ]]$loglik + ( length( mod[[ i ]]$be ) + 1 ) * con
   
     } else if ( ci_test == "testIndPois ") {
       mod[[ i ]] = glm( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), family = poisson, weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndQPois") {
       mod[[ i ]] <- glm( target ~ ., data = data.frame(dataset[, signature[i, ] ]), family = quasipoisson, weights = wei )
       bic[i] <- NA
       
     } else if ( ci_test == "testIndQBinom") {
       mod[[ i ]] <- glm( target ~ ., data = data.frame(dataset[, signature[i, ] ]), family = quasibinomial, weights = wei )
       bic[i] <- NA
       
     } else if ( ci_test == "testIndNB" ) {
       mod[[ i ]] = MASS::glm.nb( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndZIP" ) {
       mod[[ i ]] <- zip.mod( target, dataset[, signature[i, ] ], wei = wei )
       bic[i] <-  -2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1) * con

     } else if ( is.matrix(target) & ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1)  &  Rfast::Var( Rfast::rowsums(target) ) == 0 )   target = log( target[, -1]/(target[, 1]) ) 
       mod[[ i ]] = lm( target ~.,  data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
       bic[i] = NULL
       
     } else if ( is.matrix(target)  &  ci_test == "testIndBinom" ) {
       mod[[ i ]] = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, signature[i, ] ]), weights = target[, 2], family = binomial )
       bic[i] = BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndGamma" ) {
       mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, signature[i, ] ]), weights = wei, family = Gamma(link = log) )
       bic[i] <- BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndNormLog" ) {
       mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, signature[i, ] ]), weights = wei, family = gaussian(link = log) )
       bic[i] <- BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndTobit" ) {
       mod[[ i ]] <- survival::survreg( target ~ ., data = as.data.frame(dataset[, signature[i, ] ]), weights = wei, dist = "gaussian" )
       bic[i] <-  - 2 * as.numeric( logLik(mod[[ i ]]) ) + ( length( coef(mod[[ i ]]) ) + 1 ) * con

     } else if (ci_test == "censIndCR") {
       mod[[ i ]] = survival::coxph( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if (ci_test == "censIndWR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
       bic[i] = - 2 * logLik(mod[[ i ]]) + (length(mod[[ i ]]$coefficients) + 1) * con

     } else if (ci_test == "censIndER") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei, dist = "exponential" )
       bic[i] = - 2 * logLik(mod[[ i ]]) + (length(mod[[ i ]]$coefficients) + 1) * con
       
     } else if (ci_test == "censIndLLR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei, dist = "loglogistic" )
       bic[i] = - 2 * logLik(mod[[ i ]]) + (length(mod[[ i ]]$coefficients) + 1) * con
       
     } else if (ci_test == "testIndClogit") {
       case = as.logical(target[, 1]);  
       id = target[, 2]
       mod[[ i ]] = survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , signature[i, ] ] ) )  ## weights are ignored here anyway
       bic[i] = BIC( mod[[ i ]] )
       
     } else if ( ci_test == "testIndLogistic") {
       mod[[ i ]] = glm( target ~., data = data.frame(dataset[, signature[i, ] ]) , family = binomial, weights = wei ) 
       bic[[ i ]] = BIC(mod[[ i ]])

     } else if ( ci_test == "testIndMultinom" || ci_test == "gSquare" ) { 
       mod[[ i ]] = nnet::multinom( target ~., data = data.frame( dataset[, signature[i, ] ] ), trace = FALSE, weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndOrdinal" ) {
       mod[[ i ]] = ordinal::clm( target ~., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )
      
     }  
  
    } ## end for ( i in 1:nsignat ) 
     signature = cbind(signature, bic)   

  }  ## end if ( nsignat > 1 & nrow(sesObject@signatures) > 1 )
  ####### all signatures

  if ( nsignat == "all" ) { 
    signature = sesObject@signatures
    bic = numeric( nrow(signature) )
    mod = list()
    con <- log( NROW(dataset) )
	
    for ( i in 1:nrow(signature) ) {
	
     if ( ci_test == "testIndFisher" | ci_test == "testIndReg" ) {
        mod[[ i ]] = lm( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
        bic[i] = BIC( mod[[ i ]] )
      
     } else if ( ci_test == "testIndMMReg" ) {
       mod[[ i ]] = MASS::rlm(target~., data = data.frame( dataset[, signature[i, ] ] ), maxit = 2000, weights = wei)
       bic[[ i ]] = BIC( mod[[ i ]] )
       
     } else if ( ci_test == "testIndSpearman"  |  ci_test == "testIndRQ" ) {
       mod[[ i ]] = quantreg::rq( target ~., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
	     la <- logLik( mod[[ i ]] )
       bic[i] <-  - 2 * as.numeric( la ) + attr(la, "df") * con

     } else if ( ci_test == "testIndBeta" ) {
       mod[[ i ]] <- beta.mod( target, dataset[, signature[i, ] ], wei = wei )
       bic[i] <-  - 2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1 ) * con

     } else if ( ci_test == "testIndPois ") {
       mod[[ i ]] = glm( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), poisson, weights = wei )
       bic[i] = BIC( mod[[ i ]] )
       
     } else if ( ci_test == "testIndQPois") {
       mod[[ i ]] <- glm( target ~ ., data = data.frame(dataset[, signature[i, ] ]), family = quasipoisson, weights = wei )
       bic[i] <- NA
       
     } else if ( ci_test == "testIndQBinom") {
       mod[[ i ]] <- glm( target ~ ., data = data.frame(dataset[, signature[i, ] ]), family = quasibinomial, weights = wei )
       bic[i] <- NA

     } else if ( ci_test == "testIndNB" ) {
       mod[[ i ]] = MASS::glm.nb( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndZIP" ) {
       mod[[ i ]] <- zip.mod( target, dataset[, signature[i, ] ], wei = wei )
       bic[i] <-  -2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1) * log( length(target) )

     } else if ( is.matrix(target) & ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1)  &  Rfast::Var( Rfast::rowsums(target) ) == 0 )   target = log( target[, -1]/(target[, 1]) ) 
       mod[[ i ]] = lm( target ~., data = dataset[, signature[i, ]], weights = wei )
       bic[i] = NULL
       
     } else if ( is.matrix(target)  &  ci_test == "testIndBinom" ) {      
       mod[[ i ]] = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, signature[i, ] ]), weights = target[, 2], family = binomial )
       bic[i] = BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndGamma" ) {
       mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, signature[i, ] ]), weights = wei, family = Gamma(link = log) )
       bic[i] <- BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndNormLog" ) {
       mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, signature[i, ] ]), weights = wei, family = gaussian(link = log) )
       bic[i] <- BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndTobit" ) {
       mod[[ i ]] <- survival::survreg( target ~ ., data = as.data.frame(dataset[, signature[i, ] ]), weights = wei, dist = "gaussian" )
       bic[i] <-  - 2 * as.numeric( logLik(mod[[ i ]]) ) + ( length( coef(mod[[ i ]]) ) + 1 ) * con
       
     } else if (ci_test == "censIndCR") {
       mod[[ i ]] = survival::coxph( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if (ci_test == "censIndWR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
       bic[i] = - 2 * logLik(mod[[ i ]]) + (length(mod[[ i ]]$coefficients) + 1) * con
	   
     } else if (ci_test == "censIndER") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei, dist = "exponential" )
       bic[i] = BIC( mod[[ i ]] )
       
     } else if (ci_test == "censIndLLR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, signature[i, ] ] ), weights = wei, dist = "loglogistic" )
       bic[i] = - 2 * logLik(mod[[ i ]]) + (length(mod[[ i ]]$coefficients) + 1) * con
       
     } else if (ci_test == "testIndClogit") {
       case = as.logical(target[, 1]);  
       id = target[, 2]
       mod[[ i ]] = survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , signature[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndLogistic") {
       mod[[ i ]] = glm( target ~., data = data.frame(dataset[, signature[i, ] ]), family = binomial, weights = wei ) 
       bic[[ i ]] = BIC(mod[[ i ]])
  
     } else if ( ci_test == "testIndMultinom" || ci_test == "gSquare"  ) { 
       mod[[ i ]] = nnet::multinom( target ~., data = data.frame( dataset[, signature[i, ] ] ), trace = FALSE, weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( is.ordered(target) ) {
       mod[[ i ]] = ordinal::clm( target ~., data = data.frame( dataset[, signature[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )      
     }
    
  }  ## end for ( i in 1:nrow(signature) )
    
    signature = cbind(signature, bic)
	
  }  ## end if ( nsignat == "all" )

   res <- list(mod = mod, signature = signature)  
  
  }  ## end if ( sum( is.na(sesObject@selectedVars) ) > 0 ) { 

  res
}
  
 