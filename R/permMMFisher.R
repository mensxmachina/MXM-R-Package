permMMFisher = function(target, dataset, xIndex, csIndex, wei = NULL, statistic = FALSE, univariateModels=NULL, hash = FALSE, stat_hash = NULL,
                      pvalue_hash = NULL, threshold = 0.05, R = 999) {
  # TESTINDFISHER Fisher Conditional Independence Test for continous class variables
  # PVALUE = TESTINDFISHER(Y, DATA, XINDEX, CSINDEX)
  # This test provides a p-value PVALUE for the NULL hypothesis H0 which is
  # X is independent by TARGET given CS. The pvalue is calculated following
  # Fisher's method (see reference below)
  # This method requires the following inputs
  #   TARGET: a numeric vector containing the values of the target (continuous) variable. 
  #   Its support can be R or any number betweeen 0 and 1, i.e. it contains proportions.
  #   DATASET: a numeric data matrix containing the variables for performing the test. They can be only be continuous variables. 
  #   XINDEX: the index of the variable whose association with the target we want to test. 
  #   CSINDEX: the indices if the variable to condition on. 
  # this method returns: the pvalue PVALUE, the statistic STAT.
  # Copyright 2012 Vincenzo Lagani and Ioannis Tsamardinos
  # R Implementation by Giorgos Athineou (10/2013)
  #if the test cannot performed succesfully these are the returned values
  pvalue = log(1)
  stat = 0
  
  if ( !is.list(target) ) {
    
    n = length( target )
    csIndex[which(is.na(csIndex))] = 0
    
    if ( hash )  {
      csIndex2 = csIndex[which(csIndex!=0)]
      csIndex2 = sort(csIndex2)
      xcs = c(xIndex,csIndex2)
      key = paste(as.character(xcs) , collapse=" ")
      
      if (is.null(stat_hash[key]) == FALSE) {
        stat = stat_hash[key]
        pvalue = pvalue_hash[key]
        results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
        return(results)
      }
    }
    #if the xIndex is contained in csIndex, x does not bring any new information with respect to cs
    if ( !is.na( match(xIndex, csIndex) ) ) {
      if ( hash ) {        #update hash objects
        stat_hash[key] <- 0       #.set(stat_hash , key , 0)
        pvalue_hash[key] <- log(1)     #.set(pvalue_hash , key , 1)
      }
      results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
      return(results)
    }
    #check input validity
    if ( any(xIndex < 0) || any(csIndex < 0) ) {
      message(paste("error in testIndFisher : wrong input of xIndex or csIndex"))
      results <- list(pvalue = pvalue, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
      return(results)
    }
    
    xIndex = unique(xIndex)
    csIndex = unique(csIndex)
    #extract the data
    x = dataset[ , xIndex]
    cs = dataset[ , csIndex, drop = FALSE]
    #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
    if ( length(cs) != 0 ) {
      if ( is.null(dim(cs)[2]) ) {     #cs is a vector
        if ( identical(x, cs) ) {    #if(!any(x == cs) == FALSE)
          if ( hash ) {        #update hash objects
            stat_hash[key] <- 0         #.set(stat_hash , key , 0)
            pvalue_hash[key] <- log(1)        #.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
          return(results)
        }
      } else { #more than one var
        for (col in 1:dim(cs)[2]) {
          if ( identical(x, cs[, col]) )  {    #if(!any(x == cs) == FALSE)
            if ( hash )  {       #update hash objects
              stat_hash[key] <- 0      #.set(stat_hash , key , 0)
              pvalue_hash[key] <- log(1)         #.set(pvalue_hash , key , 1)
            }
            results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
            return(results)
          }
        }
      }
    }
    res <- tryCatch(
      {
        #if the conditioning set (cs) is empty, we use a simplified formula
        if (length(cs) == 0) {
          if ( !is.null(univariateModels) ) {
            pvalue = univariateModels$pvalue[[xIndex]]
            stat = univariateModels$stat[[xIndex]]
            results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
            return(results)
          }
          #compute the correlation coefficient between x,target directly
          mod1 <- MASS::rlm( target ~ x, maxit = 2000, method = "MM" )
          mod2 <- MASS::rlm( x ~ target, maxit = 2000, method = "MM" )       
          e1 <- mod1$residuals
          e2 <- mod2$residuals
          res <- permcor(e1, e2, R )
          stat <- abs( res[1] )
          pvalue <- log(res[2])
        } else {
          #perform the test with the cs
          e1 <- resid( MASS::rlm( target ~ cs, maxit = 2000, method = "MM" ) ) 
          e2 <- resid( MASS::rlm( x ~ cs, maxit = 2000, method = "MM" ) )
          res <- permcor( e1, e2, R)
          stat <- abs( res[1] )
          pvalue <- log(res[2])
        }
        #last error check
        if ( is.na(pvalue) || is.na(stat) ) {
          pvalue = log(1)
          stat = 0
        } else {
          #update hash objects
          if( hash ) {
            stat_hash[key] <- stat        #.set(stat_hash , key , stat)
            pvalue_hash[key] <- pvalue        #.set(pvalue_hash , key , pvalue)
          }
        }
        results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
        return(results)
        
      },
      error=function(cond) {
        pvalue = log(1)
        stat = 0
        results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
        return(results)
      },
      finally={}
    )  
    
    
    ## meta-analytic approach
    ##################################
    
    
  } else {  
    
    D = length(target)
    pva = numeric(D)
    
    aa <- list()
    
    for ( i in 1:D ) {
      
      targ = target[[ i ]]
      data = dataset[[ i ]]  
      
      n = length( targ )
      csIndex[which(is.na(csIndex))] = 0
      
      if ( hash ) {
        csIndex2 = csIndex[which(csIndex!=0)]
        csIndex2 = sort(csIndex2)
        xcs = c(xIndex, csIndex2)
        key = paste(as.character(xcs) , collapse=" ")
        if (is.null(stat_hash[key]) == FALSE) {
          stat = stat_hash[key]
          pvalue = pvalue_hash[key]
          aa[[ i ]] <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
        }
      }
      
      #if the xIndex is contained in csIndex, x does not bring any new
      #information with respect to cs
      if ( !is.na( match(xIndex, csIndex) ) ) { 
        if ( hash ) {    #update hash objects
          stat_hash[key] <- 0    #.set(stat_hash , key , 0)
          pvalue_hash[key] <- log(1)   #.set(pvalue_hash , key , 1)
        }
        aa[[ i ]] <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
      }
      
      #check input validity
      if ( any(xIndex < 0) || any(csIndex < 0) ) {
        message(paste("error in testIndFisher : wrong input of xIndex or csIndex"))
        aa[[ i ]] <- list(pvalue = pvalue, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
      }
      
      xIndex = unique(xIndex)
      csIndex = unique(csIndex)
      #extract the data
      x = data[ , xIndex]
      cs = data[ , csIndex, drop = FALSE]
      #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
      if ( length(cs) != 0 ) {
        if ( is.null( dim(cs )[2]) ) {  #cs is a vector
          if (any(x != cs) == FALSE)  {  #if(!any(x == cs) == FALSE)
            if ( hash )  {  #update hash objects
              stat_hash[key] <- 0   #.set(stat_hash , key , 0)
              pvalue_hash[key] <- log(1)   #.set(pvalue_hash , key , 1)
            }
            aa[[ i ]] <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
          }
        } else { #more than one var
          for ( col in 1:ncol(cs) ) {
            if (any(x != cs[, col]) == FALSE) {   #if(!any(x == cs) == FALSE)
              if ( hash ) {    #update hash objects
                stat_hash[key] <- 0       #.set(stat_hash , key , 0)
                pvalue_hash[key] <- log(1)        #.set(pvalue_hash , key , 1)
              }
              aa[[ i ]] <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
            }
          }
        }
      }
      
      #if x or target is constant then there is no point to perform the test
      if ( Rfast::Var(x)  == 0 ) {
        if ( hash )  {  #update hash objects
          stat_hash[key] <- 0    #.set(stat_hash , key , 0)
          pvalue_hash[key] <- log(1)    #.set(pvalue_hash , key , 1)
        }
        aa[[ i ]] <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
      }
      
      #remove constant columns of cs
      cs = cs[, apply(cs, 2, var, na.rm=TRUE) != 0]
      
      aa[[ i ]] <- tryCatch(
        {
          #if the conditioning set (cs) is empty, we use a simplified formula
          if (length(cs) == 0)  {
            if ( !is.null(univariateModels) ) {
              pvalue = univariateModels$pvalue[[xIndex]]
              stat = univariateModels$stat[[xIndex]]
              results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
              return(results)
            }
            #compute the correlation coefficient between x, target directly
            mod1 <- MASS::rlm( targ ~ x, maxit = 2000, method = "MM" )
            mod2 <- MASS::rlm( x ~ targ, maxit = 2000, method = "MM" )       
            e1 <- mod1$residuals
            e2 <- mod2$residuals
            res <- permcor( e1, e2, R )
            stat <- abs( res[1] )
            pvalue <- log(res[2])
          } else{
            #perform the test with the cs
            e1 = resid( MASS::rlm( targ ~., data = data.frame( data[, csIndex] ), maxit = 2000, method = "MM" ) ) 
            e2 = resid( MASS::rlm( data[, xIndex] ~.,  data = data.frame( data[, csIndex] ), maxit = 2000, method = "MM" ) )
            res <- permcor( e1, e2, R )
            stat <- abs( res[1] )
            pvalue <- log(res[2])
          }
          #last error check
          if (is.na(pvalue) || is.na(stat) ) {
            pvalue = log(1)
            stat = 0
          } else {
            
            if( hash ) {
              stat_hash[key] <- stat    #.set(stat_hash , key , stat)
              pvalue_hash[key] <- pvalue    #.set(pvalue_hash , key , pvalue)
            }
          }
          #testerrorcaseintrycatch(4)
          list(pvalue = pvalue, nu = n, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
        },
        
        error=function(cond) {
          #    message(paste("warning in try catch of the testIndFisher test"))
          #    message("Here's the original message:")
          #    message(cond)
          #    for debug
          #    print("\nxIndex = \n")
          #    print(xIndex)
          #    print("\ncsindex = \n")
          #    print(csIndex)
          #   stop()
          #error case (we are pretty sure that the only error case is when x,cs are highly correlated and the inversion of the matrix is not possible)
          pvalue = log(1)
          stat = 0
          list(pvalue = pvalue, nu = n, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
        },
        finally={}
      ) 
      
    }
    
    pva <- numeric(D) 
    for ( j in 1:D )   pva[j] = -2 * aa[[ j ]]$pvalue
    stat = sum(pva)
    pvalue = pchisq( stat, 2 * D, lower.tail = FALSE, log.p = TRUE ) 
    
    if ( hash ) {
      stat_hash[key] <- stat
      pvalue_hash[key] <- pvalue  
    }
    
    res = list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash)
  }
  
  return(res)
}