permRQ = function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, 
                     pvalue_hash=NULL, threshold = 0.05, R = 999) {
  # TESTINDRQ Conditional Independence Test for numerical class variables 
  # PVALUE = TESTINDRQ(Y, DATA, XINDEX, CSINDEX)
  # This test provides a p-value PVALUE for the NULL hypothesis H0 which is
  # X is independent by TARGET given CS. The pvalue is calculated following
  # nested models
  # This method requires the following inputs
  #   TARGET: a numeric vector containing the values of the target (numerical) variable. 
  #   Its support can be R or any number betweeen 0 and 1, i.e. it contains proportions.
  #   DATASET: a numeric data matrix containing the variables for performing the test. They can be mixed variables. 
  #   XINDEX: the index of the variable whose association with the target we want to test. 
  #   CSINDEX: the indices if the variable to condition on. 
  # this method returns: the pvalue PVALUE, the statistic STAT.
  #if the test cannot performed succesfully these are the returned values
  pvalue = log(1)
  stat = 0;
  csIndex[which(is.na(csIndex))] = 0;
  thres <- threshold * R + 1
  
  if( hash )  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csIndex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if(is.null(stat_hash[key]) == FALSE)  {
      stat = stat_hash[key];
      pvalue = pvalue_hash[key];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #if the xIndex is contained in csIndex, x does not bring any new
  #information with respect to cs
  if ( !is.na( match(xIndex, csIndex) ) )  {
    if( hash )  {       #update hash objects
      stat_hash[key] <- 0;#.set(stat_hash , key , 0)
      pvalue_hash[key] <- 1;#.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = 1, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if ( any(xIndex < 0) || any(csIndex < 0) )  {
    message(paste("error in testIndRQ : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  
  xIndex = unique(xIndex);
  csIndex = unique(csIndex);
  #extract the data
  x = dataset[ , xIndex];
  cs = dataset[ , csIndex];
  if ( length(cs) != 0 ) {
    if ( is.null(dim(cs)[2]) )  {   #cs is a vector
      if (identical(x, cs) )  {  #if(!any(x == cs) == FALSE)
        if ( hash )  {   #update hash objects
          stat_hash[key] <- 0;  #.set(stat_hash , key , 0)
          pvalue_hash[key] <- 1;  #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = 1, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for (col in 1:dim(cs)[2])  {
        if (identical(x, cs[, col]) )  {    #if(!any(x == cs) == FALSE)
          if ( hash ) {      #update hash objects
            stat_hash[key] <- 0;#.set(stat_hash , key , 0)
            pvalue_hash[key] <- 1;#.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = 1, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
  res <- tryCatch(
{
  #if the conditioning set (cs) is empty, we use a simplified formula
  if (length(cs) == 0) {
    #compute the relationship between x,target directly
    #fit1 = fita
    fit1 = quantreg::rq(target ~ 1, weights = wei)
    fit2 = quantreg::rq(target ~ x, weights = wei)
    mod = anova(fit1, fit2, test = "rank")
    stat = as.numeric( mod[[1]][3] )
    step <- 0
    j <- 1		
    n <- length(target)
    while (j <= R & step < thres ) {
      xb <- sample(x, n)  
      bit2 = quantreg::rq(target ~ xb, weights = wei )
      ww = anova(fit1, bit2, test = "rank")
      step <- step + ( as.numeric( ww[[1]][3] ) > stat )
      j <- j + 1
    }
    pvalue <- log( (step + 1) / (R + 1) )
  } else {
    fit1 = quantreg::rq( target ~ cs, weights = wei )
    fit2 = quantreg::rq( target ~ cs + x, weights = wei )
    mod = anova(fit1, fit2, test = "rank")
    stat = as.numeric( mod[[1]][3] )
    step <- 0
    j <- 1		
    n <- length(target)
    while (j <= R & step < thres ) {
      xb <- sample(x, n)  
      bit2 = quantreg::rq(target ~ cs + xb, weights = wei )
      ww = anova(fit1, bit2, test = "rank")
      step <- step + ( as.numeric( ww[[1]][3] ) > stat )
      j <- j + 1
    }
    pvalue <- log( (step + 1) / (R + 1) )
  }
  #last error check
  if ( is.na(pvalue) || is.na(stat) )  {
    pvalue = log(1)
    stat = 0;
  } else {
    #update hash objects
    if ( hash ) {
      stat_hash[key] <- stat;#.set(stat_hash , key , stat)
      pvalue_hash[key] <- pvalue;#.set(pvalue_hash , key , pvalue)
    }
  }
  #testerrorcaseintrycatch(4);
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
  
},
error=function(cond) {
  pvalue = log(1)
  stat = 0;
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
},
finally={}
  )    
  return(res);
}