fbed.glmm.reg <- function(target, dataset, id, reps = NULL, ini = NULL, threshold = 0.05, wei = NULL, K = 0, method = "LR", gam = NULL, backward = TRUE, test = "testIndGLMMReg") {
  
  if ( length(K) > 1 ) {
    
    result <- kfbed.glmm.reg(y = target, x = dataset, id = id, univ = ini, alpha = threshold, wei = wei, K = K, method = method, gam = gam, backward = backward, test = test)
    
  } else {
    
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any( is.na(dataset) ) ) {  
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
  }
  
  if ( method =="LR" ) {
    if ( test == "testIndGLMMOrdinal" ) {
      if ( !is.null(reps) ) {
        result <- fbed.glmm.ordinal.reps(y = target, x = dataset, id = id, reps = reps, univ = ini, alpha = threshold, wei = wei, K = K)
      } else  result <- fbed.glmm.ordinal(y = target, x = dataset, id = id, univ = ini, alpha = threshold, wei = wei, K = K)       
    } else if ( test == "testIndGLMMCR" ) {
      result <- fbed.glmm.cr(y = target, x = dataset, id, univ = ini, alpha = threshold, wei = wei, K = K)  
    } else if ( test == "testIndGLMMReg" ) {
      if ( !is.null(reps) ) {
        result <- fbed.lmm.reps(y = target, x = dataset, id = id, reps = reps, univ = ini, alpha = threshold, wei = wei, K = K)
      } else  result <- fbed.lmm(y = target, x = dataset, id = id, univ = ini, alpha = threshold, wei = wei, K = K) 
    } else {
      if ( !is.null(reps) ) {
        result <- fbed.glmm.reps(y = target, x = dataset, id = id, reps = reps, univ = ini, alpha = threshold, wei = wei, K = K, test = test)
      } else  result <- fbed.glmm(y = target, x = dataset, id = id, univ = ini, alpha = threshold, wei = wei, K = K, test = test)
    }  

    result$back.rem <- 0
    result$back.n.tests <- 0
    
    if ( backward ) {
      
      if (result$info[1, 1] > 0) {
        if ( test == "testIndGLMMOrdinal" ) {
          if ( !is.null(reps) ) {
            a <- glmm.ordinal.reps.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, reps = reps, threshold = threshold, wei = wei)
          } else  a <- glmm.ordinal.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, threshold = threshold, wei = wei)  
        } else if ( test == "testIndGLMMCR" ) {
          a <- glmm.cr.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, threshold = threshold, wei = wei)  
        } else if ( test == "testIndGLMMReg" ) {
          if ( !is.null(reps) ) {
            a <- lmm.reps.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, reps = reps, threshold = threshold, wei = wei)
          } else  a <- lmm.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, threshold = threshold, wei = wei)            
        } else {
          if ( !is.null(reps) ) {
            a <- glmm.reps.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, reps = reps, threshold = threshold, wei = wei, test = test)
          } else  a <- glmm.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, threshold = threshold, wei = wei, test = test)
        }
        
        if ( typeof(a) == "list" ) {
          result$back.rem <- result$res[a$info[, 1], 1]
          back.n.tests <- sum( dim(result$res)[1] : dim(a$mat)[1] )
          sel <- result$res[a$mat[, 1], 1] 
          stat <- a$mat[, 3]
          pval <- a$mat[, 2]
          result$res <- cbind(sel, stat, pval)
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime + a$runtime
        } else {
          back.rem <- 0
          back.n.tests <- 0
          result$back.rem <- back.rem
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime 
        }  ## end if ( typeof(a) == "list" ) 
      }   ## end if (result$info[1, 1] > 0)
    }  ## end if ( backward )
    
  } else {

    if (test == "testIndGLMMOrdinal") {
      if ( !is.null(reps) ) {
        result <- ebic.fbed.glmm.ordinal.reps(target, dataset, id, reps = reps, univ = ini, gam = gam, wei = wei, K = K)
      } else  result <- ebic.fbed.glmm.ordinal(target, dataset, id, univ = ini, gam = gam, wei = wei, K = K)
        
    } else if ( test == "testIndGLMMCR" ) {  
      result <- ebic.fbed.glmm.cr(target, dataset, id, univ = ini, gam = gam, wei = wei, K = K) 
    } else if (test == "testIndGLMMReg") {
      if ( !is.null(reps) ) {
        result <- ebic.fbed.lmm.reps(target, dataset, id, reps = reps, univ = ini, gam = gam, wei = wei, K = K)
      } else  result <- ebic.fbed.lmm(target, dataset, id, univ = ini, gam = gam, wei = wei, K = K)
    } else {
      if ( !is.null(reps) ) {
        result <- ebic.fbed.glmm.reps(target, dataset, id, reps = reps, univ = ini, gam = gam, wei = wei, K = K, test = test) 
      } else result <- ebic.fbed.glmm(target, dataset, id, univ = ini, gam = gam, wei = wei, K = K, test = test) 
    }  
    result$back.rem <- 0
    result$back.n.tests <- 0
    
    if ( backward ) {
      
      if (result$info[1, 1] > 0) {
        if ( test == "testIndGLMMOrdinal" ) {
          if ( !is.null(reps) ) {
            a <- ebic.glmm.ordinal.reps.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, reps = reps, wei = wei, gam = gam)
          } else  a <- ebic.glmm.ordinal.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, wei = wei, gam = gam)          
        } else if ( test == "testIndGLMMCR" ) {
          a <- ebic.glmm.cr.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, wei = wei, gam = gam) 
        } else { 
          if ( !is.null(reps) ) {
            a <- ebic.glmm.reps.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, reps = reps, wei = wei, gam = gam, test = test)
          } else  a <- ebic.glmm.bsreg(target, dataset[, result$res[, 1], drop = FALSE], id, wei = wei, gam = gam, test = test)
        }  
        if ( typeof(a) == "list" ) {
          back.n.tests <- sum( dim(result$res)[1] : length(a$mat[, 1]) )
          
          result$back.rem <- result$res[a$info[, 1], 1]
          sel <- result$res[ a$mat[, 1], 1]
          val <- a$mat[, 2]
          result$res <- cbind(sel, val)
          colnames(result$res) <- c("Vars", "eBIC")
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime + a$runtime
        } else {
          back.rem <- 0
          back.n.tests <- 0
          result$back.rem <- back.rem
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime 
        }  
      }   ## end if (result$info[1, 1] > 0) 
      
    }  ## end if ( backward )
    
  }  ## end if ( method == "LR" ) 
  
  }  ## end ( if length(K) > 1 )
  
  result
}

