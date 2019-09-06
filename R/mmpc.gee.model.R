# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of signatures and generated models. It could be numeric from 1 to total number of signatures or "all" for all the 
## signatures. Default is 1.
mmpc.gee.model = function(target, dataset, reps = NULL, group, correl = "exchangeable", se = "jack", wei = NULL, 
                          mmpcgee.Object, test = NULL) {
  
  if ( sum( is.na(mmpcgee.Object@selectedVars) ) > 0 ) {
    mod <- paste("No associations were found, hence no model is produced.")
    signature = NULL
    res <- list(mod = mod, signature = signature)  
    
  } else {

    if ( any(is.na(dataset) ) ) {
      warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    }
    
    if ( is.null(test) ) {  
      ci_test <- mmpcgee.Object@test
    } else ci_test <- test 
    
    signature <- mmpcgee.Object@selectedVars  
      
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
      
    if ( is.null( colnames(dataset) ) ) {
      names(signature) = paste("Var", signature, sep = " ")
    } else  names(signature) = colnames(dataset)[signature]
    
    res <- list(mod = mod, signature = signature)  
    
  } ## end if ( sum( is.na(mmpcgee.Object@selectedVars) ) > 0 ) 
  
  res
  
}

