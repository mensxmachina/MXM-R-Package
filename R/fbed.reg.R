fbed.reg <- function(target, dataset, ini = NULL, test = NULL, threshold = 0.05, wei = NULL, K = 0, method = "LR", gam = NULL, backward = TRUE) {
  
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any( is.na(dataset) ) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(za){ za[which(is.na(za))] = median(za, na.rm = TRUE) ; return(za) } ) 
    } else {
       poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
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
  
  dataset <- as.data.frame(dataset)
  
  if ( length(K) > 1 ) {
    
    result <- kfbed.reg(y = target, x = dataset, univ = ini, test = test, alpha = threshold, wei = NULL, K = K, method = method, gam = gam, backward = backward)
  
  } else {
   
  if (test == "gSquare") {
    
    result <- fbed.g2(y = target, x = dataset, alpha = threshold, univ = ini, K = K, backward = backward)
    
  } else {
    
  if (method == "LR") {
    
    if (test == "testIndFisher") {
      result <- Rfast::cor.fbed(target, as.matrix(dataset), alpha = threshold, K = K)
      result$univ <- NULL
      
    } else  result <- fbed.lr(y = target, x = dataset, alpha = threshold, univ = ini, test = test, wei = wei, K = K)  
    
    result$back.rem <- 0
    result$back.n.tests <- 0
    
    if ( backward ) {
      
      if (result$info[1, 1] > 0) {
        a <- MXM::bs.reg(target, dataset[, result$res[, 1], drop = FALSE], threshold = threshold, wei = wei, test = test)
        
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
        }  
      }   ## end if (result$info[1, 1] > 0)
    }  ## end if ( backward )
    
  } else {  ## end of method =="LR"
    
    if ( test == "testIndFisher" )    test <- "testIndReg"
    result <- fbed.ebic(y = target, x = dataset, test = test, univ = ini, gam = gam, wei = wei, K = K)  
    
    result$back.rem <- 0
    result$back.n.tests <- 0
    
    if ( backward ) {
      
      if ( result$info[1, 1] > 0 ) {
        a <- MXM::ebic.bsreg(target, dataset[, result$res[, 1], drop = FALSE], test = test, wei = wei, gam = gam) 
       
        if ( typeof(a) == "list" ) {
          back.n.tests <- sum( dim(result$res)[1] : length(a$mat[, 1]) )
          result$back.rem <- result$res[a$info[, 1], 1]
          sel <- result$res[ a$mat[, 1], 1]
          val <- a$mat[, 2]
          result$res <- cbind(sel, val)
          colnames(result$res) <- c("Vars", "eBIC difference")
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
    
  }  ## end of method == "eBIC"
  
  }  ## end if (test == "gSquare")
  
  }  ## end if ( length(K) > 1 )
  result
}