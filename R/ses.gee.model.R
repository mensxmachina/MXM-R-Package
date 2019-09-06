# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of signatures and generated models. It could be numeric from 1 to total number of signatures or "all" for all the 
## signatures. Default is 1.
ses.gee.model = function(target, dataset, reps = NULL, group, correl = "exchangeable", se = "jack", 
                         wei = NULL, sesgee.Object, nsignat = 1, test = NULL) {
  
  if ( sum( is.na(sesgee.Object@selectedVars) ) > 0 ) {
    mod <- paste("No associations were found, hence no model is produced.")
    signature <- NULL
    res <- list(mod = mod, signature = signature)  
    
  } else {
    
    if ( any(is.na(dataset) ) ) {
      warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    }
    
    if ( is.null(test) ) {  
      ci_test <- sesgee.Object@test
    } else ci_test <- test 
    
    if ( nsignat == 1 || ( nsignat > 1 & nrow(sesgee.Object@signatures) == 1 ) ) {
      signature <- sesgee.Object@selectedVars  
      
     
      if (test == "testIndGEEReg") {
          oiko <- gaussian
        } else if (test == "testIndGEELogistic") {
          oiko <- binomial(logit)
        } else if (test == "testIndGEEPois") {
          oiko <- poisson(log)
        } else if (test == "testIndGEEGamma") {
          oiko <- Gamma(log)
        } else if (test == "testIndGEENormLog") {
          oiko <- gaussian(log)
        } 
        
        if ( is.null(reps) ) {
          mod <- try( geepack::geeglm( target ~ dataset[, signature], family = oiko, id = group, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
        } else {
          mod <- try( geepack::geeglm( target ~ reps + dataset[, signature], family = oiko, id = group, weights = wei, corstr = correl, std.err = se ), silent = TRUE) 
        }       
        
      }
      
      if ( is.null( colnames(dataset) ) ) {
        names(signature) = paste("Var", signature, sep = " ")
      } else  names(signature) = colnames(dataset)[signature]
      signature <- c(signature, bic)
      names(signature)[length(signature)] = "bic"
      
      res <- list(mod = mod, signature = signature)  
    
    #############  more than one signatures
    if ( nsignat > 1 & nrow(sesgee.Object@signatures) > 1 ) {
      
      if ( nsignat > nrow(sesgee.Object@signatures) )  nsignat = nrow(sesgee.Object@signatures)
      
      bic <- numeric(nsignat)
      signature <- sesgee.Object@signatures[1:nsignat, , drop = FALSE] 
      mod <- list()
      
      for ( i in 1:nsignat ) {
         
          if (test == "testIndGEEReg") {
            oiko <- gaussian
          } else if(test == "testIndGEELogistic") {
            oiko <- binomial(logit)
          } else if (test == "testIndGEEPois") {
            oiko <- poisson(log)
          } else if (test == "testIndGEEGamma") {
            oiko <- Gamma(log)
          } else if (test == "testIndGEENormLog") {
            oiko <- gaussian(log)
          } 
          
          if ( is.null(reps) ) {
            mod[[ i ]] <- try( geepack::geeglm( target ~ dataset[, signature], family = oiko, id = group, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
          } else {
            mod[[ i ]] <- try( geepack::geeglm( target ~ reps + dataset[, signature], family = oiko, id = group, weights = wei, corstr = correl, std.err = se ), silent = TRUE) 
          }            
        
      }
      signature <- cbind(signature, bic)
      
    }
    
    if ( nsignat == "all" ) { 
      signature <- sesgee.Object@signatures
      bic <- numeric( nrow(signature) )
      mod <- list()
      
      for ( i in 1:nrow(signature) ) {
           
          if (test == "testIndGEEReg") {
              oiko <- gaussian
            } else if(test == "testIndGEELogistic") {
              oiko <- binomial(logit)
            } else if (test == "testIndGEEPois") {
              oiko <- poisson(log)
            } else if (test == "testIndGEEGamma") {
              oiko <- Gamma(log)
            } else if (test == "testIndGEENormLog") {
              oiko <- gaussian(log)
            } 
            
            if ( is.null(reps) ) {
              mod[[ i ]] <- try( geepack::geeglm( target ~ dataset[, signature], family = oiko, id = group, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
            } else {
              mod[[ i ]] <- try( geepack::geeglm( target ~ reps + dataset[, signature], family = oiko, id = group, weights = wei, corstr = correl, std.err = se ), silent = TRUE) 
            }       
            
      }       
      signature <- cbind(signature, bic)
      
    }
    
  } ## if ( sum( is.na(sesgee.Object@selectedVars) ) > 0 ) { 
  
  res
}
