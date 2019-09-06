ses.timeclass.model <- function(target, dataset, id, reps, wei = NULL, sestimeclass.Object, nsignat = 1) {
  
  signature <- sestimeclass.Object@signatures
  
  if ( sum( is.na(signature) ) > 0 ) {
    mod <- paste("No associations were found, hence no model is produced.")
    signature = NULL
    
  } else {
    
    if ( any(is.na(dataset) ) ) {
      warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    }

    la <- length( unique(id) )
    tar <- numeric(la)
    for(i in 1:la)   tar[i] <- unique( target[id == i] )
    target <- tar
    tar <- NULL
    
    if ( nsignat == 1 || ( nsignat > 1 & nrow(sestimeclass.Object@signatures) == 1 ) ) {
      signature <- sestimeclass.Object@selectedVars  
      
      dat <- dataset[, signature, drop = FALSE]
      x <- group.mvbetas(dat, id, reps)
      n <- 0.5 * dim(x)[1]
      xx <- cbind(x[1:n, ], x[(n + 1):(2 * n), ]) 

      if ( sestimeclass.Object@test == "testIndTimeLogistic" ) {
        mod <- glm(target ~ xx, binomial, weights = wei)
      } else  mod <- nnet::multinom(target ~ xx, weights = wei, trace = FALSE)
      
      if ( is.null( colnames(dataset) ) ) {  
        names(signature) = paste("Var", signature, sep = " ")
      } else  names(signature) = colnames(dataset)[signature]
    
      signature <- c( signature, BIC(mod) )
      names(signature)[length(signature)] = "bic" 
    
    }  ## end if ( nsignat == 1 || ( nsignat > 1 & nrow(sesglmm.Object@signatures) == 1 ) ) 
    
    if ( nsignat > 1 & nrow(sestimeclass.Object@signatures) > 1 ) {
      
      if ( nsignat > nrow(sestimeclass.Object@signatures) )  nsignat <- nrow(sestimeclass.Object@signatures)
      
      bic <- numeric(nsignat)
      signature <- sestimeclass.Object@signatures[1:nsignat, , drop = FALSE] 
      mod <- list()
     
      for ( i in 1:nsignat ) {
        
        dat <- dataset[, signature[i, ], drop = FALSE]
        x <- group.mvbetas(dat, id, reps)
        n <- 0.5 * dim(x)[1]
        xx <- cbind(x[1:n, ], x[(n + 1):(2 * n), ])
        
        if ( sestimeclass.Object@test == "testIndTimeLogistic" ) {
          mod[[ i ]] = glm( target ~ xx[, signature[i, ]] + (1|group), weights = wei, family = binomial ) 
          bic[i] = BIC( mod[[ i ]] )
        } else {
          mod[[ i ]] = nnet::multinom( target ~ xx[, signature[i, ]], weights = wei, trace = FALSE ) 
          bic[i] = BIC( mod[[ i ]] )
        }
      }  ## end for ( i in 1:nsignat )
      signature <- cbind(signature, bic)
    }  ##  end if ( nsignat > 1 & nrow(sestimeclass.Object@signatures) > 1 ) 
    
    if ( nsignat == "all" ) { 
      signature <- sestimeclass.Object@signatures
      bic <- numeric( NROW(signature) )
      mod <- list()
      
      for ( i in 1:nsignat ) {
        
        dat <- dataset[, signature[i, ], drop = FALSE]
        x <- group.mvbetas(dat, id, reps)
        n <- 0.5 * dim(x)[1]
        xx <- cbind(x[1:n, ], x[(n + 1):(2 * n), ])
        
        if ( sestimeclass.Object@test == "testIndTimeLogistic" ) {
          mod[[ i ]] = glm( target ~ xx[, signature[i, ]] + (1|group), weights = wei, family = binomial ) 
          bic[i] = BIC( mod[[ i ]] )
        } else {
          mod[[ i ]] = nnet::multinom( target ~ xx[, signature[i, ]], weights = wei, trace = FALSE ) 
          bic[i] = BIC( mod[[ i ]] )
        }
      }  ## end for ( i in 1:nsignat )
      signature <- cbind(signature, bic)
    }
  
  }  ## end if ( sum( is.na(signature) ) > 0 )
  
  list(mod = mod, signature = signature)  
}