testIndMMFisher = function(target, dataset, xIndex, csIndex, wei = NULL, statistic = FALSE, univariateModels=NULL , hash = FALSE, stat_hash = NULL,
 pvalue_hash = NULL) {
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
  # References
  # [1] Peter Spirtes, Clark Glymour, and Richard Scheines. Causation,
  # Prediction, and Search. The MIT Press, Cambridge, MA, USA, second
  # edition, January 2001.
  # Copyright 2012 Vincenzo Lagani and Ioannis Tsamardinos
  # R Implementation by Giorgos Athineou (10/2013)
  #if the test cannot performed succesfully these are the returned values
  pvalue = log(1);
  stat = 0;
  
  if ( !is.list(target) ) {
    n = length( target )
    csIndex[which(is.na(csIndex))] = 0;
    
    if ( hash )  {
      csIndex2 = csIndex[which(csIndex!=0)]
      csIndex2 = sort(csIndex2)
      xcs = c(xIndex,csIndex2)
      key = paste(as.character(xcs) , collapse=" ");
      
      if (is.null(stat_hash[key]) == FALSE) {
        stat = stat_hash[key];
        pvalue = pvalue_hash[key];
        results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    }
    #if the xIndex is contained in csIndex, x does not bring any new
    #information with respect to cs
    if ( !is.na( match(xIndex, csIndex) ) ) {
      if( hash )  {   #update hash objects
        stat_hash[key] <- 0;  #.set(stat_hash , key , 0)
        pvalue_hash[key] <- log(1);  #.set(pvalue_hash , key , 1)
      }
      results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
    
    #check input validity
    if( any(xIndex < 0) || any(csIndex < 0) ) {
      message(paste("error in testIndFisher : wrong input of xIndex or csIndex"))
      results <- list(pvalue = pvalue, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
    
    xIndex = unique(xIndex);
    csIndex = unique(csIndex);
    x = dataset[ , xIndex];
    cs = dataset[ , csIndex, drop = FALSE];
    #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
    if ( length(cs) != 0 ) {
      if ( is.null(dim(cs)[2]) )  {    #cs is a vector
        if (identical(x, cs) )  {    #if(!any(x == cs) == FALSE)
          if ( hash )  {   #update hash objects
            stat_hash[key] <- 0;#.set(stat_hash , key , 0)
            pvalue_hash[key] <- log(1);#.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      } else { #more than one var
        for (col in 1:dim(cs)[2])  {
          if ( identical(x, cs[, col]) )  { #if(!any(x == cs) == FALSE)
            if ( hash )  {     #update hash objects
              stat_hash[key] <- 0;    #.set(stat_hash , key , 0)
              pvalue_hash[key] <- log(1);    #.set(pvalue_hash , key , 1)
            }
            results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
            return(results);
          }
        }
      }
    }
    res <- tryCatch(
      {
        #if the conditioning set (cs) is empty, we use a simplified formula
        if ( length(cs) == 0 ) {
          if ( !is.null(univariateModels) ) {
            pvalue = univariateModels$pvalue[[xIndex]];
            stat = univariateModels$stat[[xIndex]];
            results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
            return(results);
          }
          #compute the correlation coefficient between x,target directly
           b1 = coef( MASS::rlm(target ~ x, maxit = 2000, method = "MM" ) )[2]
           b2 = coef( MASS::rlm(x ~ target, maxit = 2000, method = "MM" ) )[2]
           stat = sqrt( abs (b1 * b2) ) 
        } else {
          e1 = resid( MASS::rlm( target ~., data = data.frame( cs ), maxit = 2000, method = "MM" ) ) 
          e2 = resid( MASS::rlm( x ~.,  data = data.frame( cs ), maxit = 2000, method = "MM" ) )
          stat = cor(e1, e2) 
        }
        #lets calculate the p-value
        z = 0.5*log( (1+stat)/(1-stat) );
        dof = n - ncol( as.matrix(cs) ) - 3; #degrees of freedom
        stat = sqrt(dof) * abs(z);
        pvalue = log(2) + pt(stat, dof, lower.tail = FALSE, log.p = TRUE) ;  # ?dt for documentation
        #last error check
        if ( is.na(pvalue) || is.na(stat) ) {
          pvalue = log(1);
          stat = 0;
        } else {
          #update hash objects
          if( hash ) {
            stat_hash[key] <- stat;    #.set(stat_hash , key , stat)
            pvalue_hash[key] <- pvalue;   #.set(pvalue_hash , key , pvalue)
          }
        }
        #testerrorcaseintrycatch(4);
        results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      },
      error=function(cond) {
        #error case (we are pretty sure that the only error case is when x,cs are highly correlated and the inversion of the matrix is not possible)
        pvalue = log(1);
        stat = 0;
        results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      },
      finally={}
    )  
  
##########################
##  Meta-analytic approach
##################################
  
  
} else {  

  D = length(target)
  pva = numeric(D)
  aa <- list()
  for ( i in 1:D ) {
  
  targ = target[[ i ]]
  data = dataset[[ i ]]  
  #if the test cannot performed succesfully these are the returned values
  n = length( targ )
  csIndex[which(is.na(csIndex))] = 0;
  
  if( hash )  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csIndex2 = sort(csIndex2)
    xcs = c(xIndex, csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if(is.null(stat_hash[key]) == FALSE) {
      stat = stat_hash[key];
      pvalue = pvalue_hash[key];
      aa[[ i ]] <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    }
  }
  #if the xIndex is contained in csIndex, x does not bring any new
  #information with respect to cs
  if ( !is.na( match(xIndex, csIndex) ) ) {
    if( hash )  {    #update hash objects
      stat_hash[key] <- 0;#.set(stat_hash , key , 0)
      pvalue_hash[key] <- log(1);#.set(pvalue_hash , key , 1)
    }
    aa[[ i ]] <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  }
  #check input validity
  if( any(xIndex < 0) || any(csIndex < 0) ) {
    message(paste("error in testIndFisher : wrong input of xIndex or csIndex"))
    aa[[ i ]] <- list(pvalue = pvalue, z = 0, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  }
  
  xIndex = unique(xIndex);
  csIndex = unique(csIndex);
  x = data[ , xIndex];
  cs = data[ , csIndex];
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs) != 0 )  {
    if ( is.null( dim(cs )[2]) ) {    #cs is a vector
      if (any(x != cs) == FALSE) {   #if(!any(x == cs) == FALSE)
        if ( hash )  {    #update hash objects
          stat_hash[key] <- 0;   #.set(stat_hash , key , 0)
          pvalue_hash[key] <- log(1);   #.set(pvalue_hash , key , 1)
        }
        aa[[ i ]] <- list(pvalue = log(1), z = 0, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      }
    } else {   #more than one var
      for ( col in 1:ncol(cs) ) {
        if (any(x != cs[, col]) == FALSE)  {        #if(!any(x == cs) == FALSE)
          if ( hash )  {    #update hash objects
            stat_hash[key] <- 0;#.set(stat_hash , key , 0)
            pvalue_hash[key] <- log(1);#.set(pvalue_hash , key , 1)
          }
          aa[[ i ]] <- list(pvalue = log(1), z = 0, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        }
      }
    }
  }
  #if x or target is constant then there is no point to perform the test
  if ( Rfast::Var(x) == 0 ) {
    if( hash )  {   #update hash objects
      stat_hash[key] <- 0;       #.set(stat_hash , key , 0)
      pvalue_hash[key] <- log(1);     #.set(pvalue_hash , key , 1)
    }
    aa[[ i ]] <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  }
  #remove constant columns of cs
  cs = as.matrix(cs)
  cs = cs[, apply(cs, 2, var, na.rm=TRUE) != 0]
  
aa[[ i ]] <- tryCatch(
{
  #if the conditioning set (cs) is empty, we use a simplified formula
  if (length(cs) == 0)  {
    if ( !is.null(univariateModels) )  {
      pvalue = univariateModels$pvalue[[xIndex]];
      stat = univariateModels$stat[[xIndex]];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
    #compute the correlation coefficient between x, target directly
      b1 = coef( MASS::rlm(targ ~ x, method = "MM", maxit = 2000 ) )[2]
      b2 = coef( MASS::rlm(x ~ targ, method = "MM", maxit = 2000 ) )[2]
      stat = sqrt( abs (b1 * b2) ) 
    
  } else{
     #perform the test with the cs
      e1 = resid( MASS::rlm( targ ~., data = data.frame( cs ), method = "MM", maxit = 2000 ) ) 
      e2 = resid( MASS::rlm( x ~.,  data = data.frame( cs ), method = "MM", maxit = 2000 ) )
      stat = cor(e1, e2)
  }
  #lets calculate the p-value
  z = 0.5 * log( (1 + stat) / (1 - stat) );
  dof = n - NCOL(cs) - 3; #degrees of freedom
  stat = sqrt(dof) * abs(z) ; ## standard errot for Spearman
  pvalue = log(2) + pt(-stat, dof, log.p = TRUE) ;  # ?dt for documentation
  #last error check
  if ( is.na(pvalue) || is.na(stat) ) {
    pvalue = log(1);
    stat = 0;
  } else {
    #update hash objects
    if ( hash ) {
      stat_hash[key] <- stat; #.set(stat_hash , key , stat)
      pvalue_hash[key] <- pvalue; #.set(pvalue_hash , key , pvalue)
    }
  }
  #testerrorcaseintrycatch(4);
   list(pvalue = pvalue, z = z, nu = n, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
},

error=function(cond) {
  #error case (we are pretty sure that the only error case is when x,cs are highly correlated and the inversion of the matrix is not possible)
  pvalue = log(1);
  stat = 0;
  list(pvalue = pvalue, z = z, nu = n, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
 },
finally={}
) 
  
}

  if  ( !statistic ) {
    
	  pva <- numeric(D) 
    for ( j in 1:D )   pva[j] = -2 * aa[[ j ]]$pvalue
    stat = sum(pva)
    pvalue = pchisq( stat, 2 * D, lower.tail = FALSE, log.p = TRUE ) 
  
  } else {
    sta <- se <- numeric(D) 
	  cisa = ncol(cs)
    for ( j in 1:D )  {
	    sta[j] = aa[[ j ]]$z
	    se[j] = 1 / sqrt(aa[ j ]$nu -  cisa - 3 ) 
    }
  	sse <- sum(se)
	  stat <- (sta * se) / sqrt( sse )
	  pvalue <- log(2) + pnorm( -abs(stat), log.p = TRUE )
  }
  
  if ( hash ) {
    stat_hash[key] <- stat
    pvalue_hash[key] <- pvalue  
  }

 res = list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
}
  
 return(res)
}