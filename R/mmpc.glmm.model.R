# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of signatures and generated models. It could be numeric from 1 to total number of signatures or "all" for all the 
## signatures. Default is 1.
mmpc.glmm.model = function(target, dataset, reps = NULL, group, slopes = FALSE, wei = NULL, mmpcglmm.Object, test = NULL) {
  
  if ( sum( is.na(mmpcglmm.Object@selectedVars) ) > 0 ) {
    mod <- paste("No associations were found, hence no model is produced.")
    signature <- NULL
    res <- list(mod = mod, signature = signature)  
    
  } else {
  
  if ( any(is.na(dataset) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
  }
  
  if ( is.null(test) ) {  
    ci_test <- mmpcglmm.Object@test
    slopes <- mmpcglmm.Object@slope
  } else ci_test <- test 
  
  signature <- mmpcglmm.Object@selectedVars  

    if ( ci_test == "testIndGLMMLogistic" ) {
      if ( is.null(reps) ) {
        mod <- lme4::glmer( target ~ dataset[, signature] + (1|group), weights = wei, family = binomial ) 
      } else {
        reps <- reps 
        if (slopes ) {
          mod <- lme4::glmer( target ~ reps + dataset[, signature] + (reps|group), weights = wei, family = binomial )
        } else  mod <- lme4::glmer( target ~ reps + dataset[, signature] + (1|group), weights = wei, family = binomial ) 
      }
      
    } else if ( ci_test == "testIndGLMMPois" )  {  
      if ( is.null(reps) ) {
        mod <- lme4::glmer( target ~ dataset[, signature] + (1|group), weights = wei, family = poisson ) 
      } else {
        reps <- reps 
        if (slopes ) {
          mod <- lme4::glmer( target ~ reps + dataset[, signature] + (reps|group), weights = wei, family = poisson )
        } else  mod <- lme4::glmer( target ~ reps + dataset[, signature] + (1|group), weights = wei, family = poisson ) 
      }
      
    } else if ( ci_test == "testIndGLMMGamma" )  {  
      if ( is.null(reps) ) {
        mod <- lme4::glmer( target ~ dataset[, signature] + (1|group), weights = wei, family = Gamma(log) ) 
      } else {
        reps <- reps 
        if (slopes ) {
          mod <- lme4::glmer( target ~ reps + dataset[, signature] + (reps|group), weights = wei, family = Gamma(log) )
        } else  mod <- lme4::glmer( target ~ reps + dataset[, signature] + (1|group), weights = wei, family = Gamma(log) ) 
      }
      
    } else if ( ci_test == "testIndGLMMNormLog" )  {  
      if ( is.null(reps) ) {
        mod <- lme4::glmer( target ~ dataset[, signature] + (1|group), weights = wei, family = gaussian(log) ) 
      } else {
        reps <- reps 
        if (slopes ) {
          mod <- lme4::glmer( target ~ reps + dataset[, signature] + (reps|group), weights = wei, family = gaussian(log) )
        } else  mod <- lme4::glmer( target ~ reps + dataset[, signature] + (1|group), weights = wei, family = gaussian(log) ) 
      }
      
    } else if ( ci_test == "testIndGLMMOrdinal" )  {  
       mod <- ordinal::clmm( target ~ dataset[, signature] + (1|group), weights = wei )
       
    } else if ( ci_test == "testIndGLMMCR" )  {  
      mod <- coxme::coxme( target ~ dataset[, signature] + (1|group), weights = wei )
       
    } else if ( ci_test == "testIndGLMMReg"  |  ci_test == "testIndLMM" ) {
      if ( is.null(reps) ) {
        mod <- lme4::lmer( target ~ dataset[, signature] + (1|group), weights = wei, REML = FALSE )
      } else {
        reps <- reps 
        if ( slopes ) {
          mod <- lme4::lmer( target ~ reps + dataset[, signature] + (reps|group), weights = wei, REML = FALSE ) 
        } else {
          reps <- reps 
          mod <- lme4::lmer( target ~ reps + dataset[, signature] + (1|group), weights = wei, REML = FALSE )        
        }
      }   
    }
    
    bic <- BIC(mod)
    if ( is.null( colnames(dataset) ) ) {
      names(signature) = paste("Var", signature, sep = " ")
    } else  names(signature) = colnames(dataset)[signature]
    signature <- c(signature, bic)
    names(signature)[length(signature)] = "bic"
    
    res <- list(mod = mod, signature = signature)  
    
  } ## if ( sum( is.na(mmpcglmm.Object@selectedVars) ) > 0 ) {
    
}

