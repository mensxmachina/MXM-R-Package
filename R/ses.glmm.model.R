# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of signatures and generated models. It could be numeric from 1 to total number of signatures or "all" for all the 
## signatures. Default is 1.
ses.glmm.model = function(target, dataset, reps = NULL, group, slopes = FALSE, wei = NULL, sesglmm.Object, nsignat = 1, test = NULL) {
  
  if ( sum( is.na(sesglmm.Object@selectedVars) ) > 0 ) {
    mod <- paste("No associations were found, hence no model is produced.")
    signature <- NULL
    res <- list(mod = mod, signature = signature)  
    
  } else {
  
  if ( any(is.na(dataset) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
  }
  
  if ( is.null(test) ) {  
    ci_test <- sesglmm.Object@test
    slopes <- sesglmm.Object@slope
  } else ci_test = test 
  
  if ( nsignat == 1 || ( nsignat > 1 & nrow(sesglmm.Object@signatures) == 1 ) ) {
    signature <- sesglmm.Object@selectedVars  

  if ( test == "testIndGLMMLogistic" ) {
    if ( is.null(reps) ) {
      mod <- lme4::glmer( target ~ dataset[, signature] + (1|group), weights = wei, family = binomial ) 
    } else {
      reps <- reps 
      if (slopes ) {
        mod <- lme4::glmer( target ~ reps + dataset[, signature] + (reps|group), weights = wei, family = binomial )
      } else  mod <- lme4::glmer( target ~ reps + dataset[, signature] + (1|group), weights = wei, family = binomial ) 
    }
    
  } else if ( test == "testIndGLMMPois" )  {  
    if ( is.null(reps) ) {
      mod <- lme4::glmer( target ~ dataset[, signature] + (1|group), weights = wei, family = poisson ) 
    } else {
      reps <- reps 
      if (slopes ) {
        mod <- lme4::glmer( target ~ reps + dataset[, signature] + (reps|group), weights = wei, family = poisson )
      } else  mod <- lme4::glmer( target ~ reps + dataset[, signature] + (1|group), weights = wei, family = poisson ) 
    }
    
  } else if ( test == "testIndGLMMOrdinal" )  {  
    mod <- ordinal::clmm( target ~ dataset[, signature] + (1|group), weights = wei )
    
  } else if ( test == "testIndGLMMCR" )  {  
    mod <- coxme::coxme( target ~ dataset[, signature] + (1|group), weights = wei )
    
  } else if ( test == "testIndGLMMReg"  |  test == "testIndLMM" ) {
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
  
  if ( is.null( colnames(dataset) ) ) {
    names(signature) = paste("Var", signature, sep = " ")
  } else  names(signature) = colnames(dataset)[signature]
  signature <- c( signature, BIC(mod) )
  names(signature)[length(signature)] <- "bic"
  
  res <- list(mod = mod, signature = signature)  
  
  }  ## end if ( nsignat == 1 || ( nsignat > 1 & nrow(sesglmm.Object@signatures) == 1 ) ) 
  
  #############  more than one signatures
  if ( nsignat > 1 & nrow(sesglmm.Object@signatures) > 1 ) {
    
    if ( nsignat > nrow(sesglmm.Object@signatures) )  nsignat <- nrow(sesglmm.Object@signatures)
    
    bic <- numeric(nsignat)
    signature <- sesglmm.Object@signatures[1:nsignat, , drop = FALSE] 
    mod <- list()
    
    for ( i in 1:nsignat ) {
      if ( test == "testIndGLMMLogistic" ) {
        if ( is.null(reps) ) {
          mod[[ i ]] = lme4::glmer( target ~ dataset[, signature[i, ]] + (1|group), weights = wei, family = binomial ) 
          bic[i] = BIC( mod[[ i ]] )
        } else {
          reps = reps 
          if (slopes ) {
            mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (reps|group), weights = wei,family = binomial )
          } else  mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (1|group), weights = wei, family = binomial ) 
          bic[i] = BIC( mod[[ i ]] )
        }
        
      } else if ( test == "testIndGLMMPois" )  {  
        if ( is.null(reps) ) {
          mod[[ i ]] = lme4::glmer( target ~ dataset[, signature[i, ]] + (1|group), weights = wei, family = poisson ) 
          bic[i] = BIC( mod[[ i ]] )
        } else {
          reps = reps 
          if (slopes ) {
            mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (reps|group), weights = wei, family = poisson )
          } else  mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (1|group), weights = wei, family = poisson ) 
            bic[i] = BIC( mod[[ i ]] )
        }
        
      } else if ( test == "testIndGLMMGamma" )  {  
        if ( is.null(reps) ) {
          mod[[ i ]] = lme4::glmer( target ~ dataset[, signature[i, ]] + (1|group), weights = wei, family = Gamma(log) ) 
          bic[i] = BIC( mod[[ i ]] )
        } else {
          reps = reps 
          if (slopes ) {
            mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (reps|group), weights = wei, family = Gamma(log) )
          } else  mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (1|group), weights = wei, family = Gamma(log) ) 
            bic[i] = BIC( mod[[ i ]] )
        }
        
      } else if ( test == "testIndGLMMNormLog" )  {  
        if ( is.null(reps) ) {
          mod[[ i ]] = lme4::glmer( target ~ dataset[, signature[i, ]] + (1|group), weights = wei, family = gaussian(log) ) 
          bic[i] = BIC( mod[[ i ]] )
        } else {
          reps = reps 
          if (slopes ) {
            mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (reps|group), weights = wei, family = gaussian(log) )
          } else  mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (1|group), weights = wei, family = gaussian(log) ) 
            bic[i] = BIC( mod[[ i ]] )
        }
        
      } else if ( test == "testIndGLMMOrdinal" )  {  
        mod[[ i ]] <- ordinal::clmm( target ~ dataset[, signature] + (1|group), weights = wei )
        bic[i] <- BIC( mod[[ i ]] )
        
      } else if ( test == "testIndGLMMCR" )  {  
        mod[[ i ]] <- coxme::coxme( target ~ dataset[, signature] + (1|group), weights = wei )
        bic[i] <- BIC( mod[[ i ]] )
        
      } else if ( test == "testIndGLMMReg"  |  test == "testIndLMM" ) {
        if ( is.null(reps) ) {
          mod[[ i ]] = lme4::lmer( target ~ dataset[, signature[i, ]] + (1|group), weights = wei, REML = FALSE )
          bic[i] = BIC( mod[[ i ]] )
        } else {
          reps = reps 
          if ( slopes ) {
            mod[[ i ]] = lme4::lmer( target ~ reps + dataset[, signature[i, ]] + (reps|group), weights = wei, REML = FALSE ) 
            bic[i] = BIC( mod[[ i ]] )
          } else {
            reps = reps 
            mod[[ i ]] = lme4::lmer( target ~ reps + dataset[, signature[i, ]] + (1|group), weights = wei, REML = FALSE )        
            bic[i] = BIC( mod[[ i ]] )
          }
        }   
      }
      
    }
    signature = cbind(signature, bic)
    res <- list(mod = mod, signature = signature)  
    
  }  ## end if ( nsignat > 1 & nrow(sesglmm.Object@signatures) > 1 )
  
  if ( nsignat == "all" ) { 
    signature = sesglmm.Object@signatures
    bic = numeric( NROW(signature) )
    mod = list()

      for ( i in 1:nsignat ) {
        if ( test == "testIndGLMMLogistic" ) {
          if ( is.null(reps) ) {
            mod[[ i ]] = lme4::glmer( target ~ dataset[, signature[i, ]] + (1|group), weights = wei, family = binomial ) 
            bic[i] = BIC( mod[[ i ]] )
          } else {
            reps = reps 
            if (slopes ) {
              mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (reps|group), weights = wei, family = binomial )
            } else  mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (1|group), weights = wei, family = binomial ) 
              bic[i] = BIC( mod[[ i ]] )
          }
          
        } else if ( test == "testIndGLMMPois" )  {  
          if ( is.null(reps) ) {
            mod[[ i ]] = lme4::glmer( target ~ dataset[, signature[i, ]] + (1|group), weights = wei, family = poisson ) 
            bic[i] = BIC( mod[[ i ]] )
          } else {
            reps = reps 
            if (slopes ) {
              mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (reps|group), weights = wei, family = poisson )
            } else  mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature] + (1|group), weights = wei, family = poisson ) 
              bic[i] = BIC( mod[[ i ]] )
          }
          
        } else if ( test == "testIndGLMMGamma" )  {  
          if ( is.null(reps) ) {
            mod[[ i ]] = lme4::glmer( target ~ dataset[, signature[i, ]] + (1|group), weights = wei, family = Gamma(log) ) 
            bic[i] = BIC( mod[[ i ]] )
          } else {
            reps = reps 
            if (slopes ) {
              mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (reps|group), weights = wei, family = Gamma(log) )
            } else  mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature] + (1|group), weights = wei, family = Gamma(log) ) 
              bic[i] = BIC( mod[[ i ]] )
          }
          
        } else if ( test == "testIndGLMMNormLog" )  {  
          if ( is.null(reps) ) {
            mod[[ i ]] = lme4::glmer( target ~ dataset[, signature[i, ]] + (1|group), weights = wei, family = gaussian(log) ) 
            bic[i] = BIC( mod[[ i ]] )
          } else {
            reps = reps 
            if (slopes ) {
              mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature[i, ]] + (reps|group), weights = wei, family = gaussian(log) )
            } else  mod[[ i ]] = lme4::glmer( target ~ reps + dataset[, signature] + (1|group), weights = wei, family = gaussian(log) ) 
              bic[i] = BIC( mod[[ i ]] )
          }
          
        } else if ( test == "testIndGLMMOrdinal" )  {  
          mod[[ i ]] <- ordinal::clmm( target ~ dataset[, signature] + (1|group), weights = wei )
          bic[i] <- BIC( mod[[ i ]] )
          
        } else if ( test == "testIndGLMMCR" )  {  
          mod[[ i ]] <- coxme::coxme( target ~ dataset[, signature] + (1|group), weights = wei )
          bic[i] <- BIC( mod[[ i ]] )
          
        } else if ( test == "testIndGLMMReg"  |  test == "testIndLMM" ) {
          if ( is.null(reps) ) {
            mod[[ i ]] = lme4::lmer( target ~ dataset[, signature[i, ]] + (1|group), weights = wei, REML = FALSE )
            bic[i] = BIC( mod[[ i ]] )
          } else {
            reps = reps 
            if ( slopes ) {
              mod[[ i ]] = lme4::lmer( target ~ reps + dataset[, signature[i, ]] + (reps|group), weights = wei, REML = FALSE ) 
              bic[i] = BIC( mod[[ i ]] )
            } else {
              reps = reps 
              mod[[ i ]] = lme4::lmer( target ~ reps + dataset[, signature[i, ]] + (1|group), weights = wei, REML = FALSE )        
              bic[i] = BIC( mod[[ i ]] )
            }
          }   
        }
        
      }
    
      signature = cbind(signature, bic)
      res <- list(mod = mod, signature = signature)  
      
  }  ## end if ( nsignat == "all") 
  
  } ## if ( sum( is.na(sesglmm.Object@selectedVars) ) > 0 ) { 
  
  res
}

