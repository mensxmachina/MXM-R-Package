# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of signatures and generated models. It could be numeric from 1 to total number of signatures or "all" for all the 
## signatures. Default is 1.
waldmmpc.model = function(target, dataset, wei = NULL, wald.mmpcObject, test = NULL) {
  
  if ( sum( is.na(wald.mmpcObject@selectedVars) ) > 0 ) {
    mod = paste("No associations were found, hence no model is produced.")
    signature = NULL
    bic = NULL
  }
  
  if ( any(is.na(dataset) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    dataset = apply(dataset, 2, function(x) { x[ which(is.na(x)) ] = median(x,na.rm = TRUE) } );
  }
  
  if ( is.null(test) ) {  
    ci_test = wald.mmpcObject@test
  } else ci_test = test 
  
  rob <- FALSE
  signature <- wald.mmpcObject@selectedVars  
  p <- length(signature)

  if ( ci_test == "waldMMReg" ) {
    if ( min(target) > 0 & max(target) < 1 )   target <- log( target/(1 - target) ) 
      mod = MASS::rlm(target ~., data = data.frame(dataset[, signature ]), maxit = 2000, weights = wei, method = "MM" )
      bic = BIC( mod )        

  } else if ( ci_test == "waldBeta" ) {
    mod <- beta.mod( target, dataset[, signature ], wei = wei )
    bic <-  - 2 * mod$loglik + ( length( coef(mod$be) ) + 1 ) * log( length(target) )

  } else if ( ci_test == "waldPois") {
    mod <- glm( target ~ ., data = as.data.frame(dataset[, signature ]), family = poisson, weights = wei )
    bic <- BIC( mod )
    
  } else if ( ci_test == "waldNB" ) {
    mod <- MASS::glm.nb( target ~ ., data = as.data.frame(dataset[, signature ]), weights = wei )
    bic <- BIC(mod)
    
  } else if ( ci_test == "waldZIP" ) {
    mod <- zip.mod( target, dataset[, signature], wei = wei )
    bic <-  -2 * mod$loglik + ( p + 1) * log( length(target) )
    
  } else if ( ci_test == "waldBinom" ) {
    mod <- glm( target[, 1] /target[, 2] ~., data = as.data.frame(dataset[, signature ]), weights = target[, 2], family = binomial )
    bic <- BIC(mod)
    
  } else if ( ci_test == "waldGamma" ) {
    mod <- glm( target ~ ., data = as.data.frame(dataset[, signature ]), weights = wei, family = Gamma(link = log) )
    bic <- BIC(mod)
    
  } else if ( ci_test == "waldNormLog" ) {
    mod <- glm( target ~ ., data = as.data.frame(dataset[, signature ]), weights = wei, family = gaussian(link = log) )
    bic <- BIC(mod)
    
  } else if ( ci_test == "waldTobit" ) {
    mod <- survival::survreg( target ~ ., data = as.data.frame(dataset[, signature ]), weights = wei, dist = "gaussian" )
    bic <-  - 2 * as.numeric( logLik(mod) ) + ( p + 1 ) * log( dim(dataset)[1] )
    
  } else if (ci_test == "waldCR") {
    mod <- survival::coxph( target ~ ., data = as.data.frame(dataset[, signature ]), weights = wei )
    bic <- BIC(mod)
    
  } else if (ci_test == "waldWR") {
    mod <- survival::survreg( target ~ ., data = as.data.frame(dataset[, signature ]), weights = wei )
    bic <- BIC(mod)
    
  } else if (ci_test == "waldER") {
    mod <- survival::survreg( target ~ ., data = as.data.frame(dataset[, signature ]), weights = wei, dist = "exponential" )
    bic <- BIC(mod)

  } else if ( ci_test == "waldBinary" ) {
    mod <- glm( target ~., data = as.data.frame(dataset[, signature ]), family = binomial, weights = wei ) 
    bic <- BIC(mod)
    
  } else if ( ci_test == "waldOrdinal" ) {
    mod <- ordinal::clm( target ~., data = as.data.frame(dataset[, signature ]), weights = wei ) 
    bic <- BIC(mod)
    
  } else if (ci_test == "waldIGreg") {
    mod <- glm( target ~ ., data = as.data.frame(dataset[, signature ]), family = inverse.gaussian(log), weights = wei )
    bic <- BIC(mod)
    
  }
  
  if ( is.null( colnames(dataset) ) ) {
    names(signature) = paste("Var", signature, sep = " ")
  } else  names(signature) = colnames(dataset)[signature]
  
  signature <- c(signature, bic)
  names(signature)[length(signature)] = "bic"
  
  list(mod = mod, signature = signature)  
  
}

