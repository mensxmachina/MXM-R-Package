testIndPois = function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels = NULL , hash = FALSE, stat_hash = NULL, 
 pvalue_hash = NULL) {
  # TESTINDPOIS Conditional Independence Test for discrete class variables 
  # PVALUE = TESTINDPOIS(Y, DATA, XINDEX, CSINDEX)
  # This test provides a p-value PVALUE for the NULL hypothesis H0 which is
  # X is independent by TARGET given CS. The pvalue is calculated following
  # nested models
  # This method requires the following inputs
  #   TARGET: a numeric vector containing the values of the target (discrete) variable. 
  #   Its support can be R or any number betweeen 0 and 1, i.e. it contains proportions.
  #   DATASET: a numeric data matrix containing the variables for performing the test. They can be mixed variables. 
  #   XINDEX: the index of the variable whose association with the target we want to test. 
  #   CSINDEX: the indices if the variable to condition on. 
  # this method returns: the pvalue PVALUE, the statistic STAT.
  # References
  # [1] McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989.
  #initialization
  #if the test cannot performed succesfully these are the returned values
  pvalue = log(1);
  stat = 0;
  csIndex[which(is.na(csIndex))] = 0;
  
  if( hash )  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csIndex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if( !is.null(stat_hash[key]) ){
      stat = stat_hash[key];
      pvalue = pvalue_hash[key];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  
  if ( !is.na( match(xIndex, csIndex) ) ) {
    if( hash )  {      #update hash objects
      stat_hash[key] <- 0;    #.set(stat_hash , key , 0)
      pvalue_hash[key] <- log(1);   #.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if( any(xIndex < 0) || any(csIndex < 0) ) {
    message(paste("error in testIndPois : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  xIndex = unique(xIndex);
  csIndex = unique(csIndex);
  x = dataset[ , xIndex];
  cs = dataset[ , csIndex];
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs)!=0 )   {
    if ( is.null(dim(cs)[2]) ) {    #cs is a vector
      if ( identical(x, cs) )  {   #if(!any(x == cs) == FALSE)
        if ( hash ) {    #update hash objects
          stat_hash[key] <- 0;    #.set(stat_hash , key , 0)
          pvalue_hash[key] <- log(1);     #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for (col in 1:dim(cs)[2]) {
        if ( identical(x, cs[, col]) ) {    #if(!any(x == cs) == FALSE
          if ( hash ) {      #update hash objects
            stat_hash[key] <- 0;   #.set(stat_hash , key , 0)
            pvalue_hash[key] <- log(1);   #.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
  if (length(cs) == 0)  {
    fit2 = glm(target ~ x, poisson, weights = wei)
    stat = fit2$null.deviance - fit2$deviance
    dof = length( coef(fit2) ) - 1   
  } else {
    fit1 = glm(target ~., data = as.data.frame( cs), poisson, weights = wei, model = FALSE)
    fit2 = glm(target ~., data = as.data.frame( dataset[, c(csIndex, xIndex)] ), poisson, weights = wei, model = FALSE)
    stat = fit1$deviance - fit2$deviance
    dof = length( coef(fit2) ) - length( coef(fit1) )
  } 
    pvalue = pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE )  
  #last error check
  if ( is.na(pvalue) || is.na(stat) )   {
    pvalue = log(1);
    stat = 0;
  } else {
    #update hash objects
    if ( hash )  {
      stat_hash[key] <- stat;      #.set(stat_hash , key , stat)
      pvalue_hash[key] <- pvalue;        #.set(pvalue_hash , key , pvalue)
    }
  }
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results)  
}