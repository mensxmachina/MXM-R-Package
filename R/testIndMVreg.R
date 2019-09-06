testIndMVreg = function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels=NULL , hash = FALSE, stat_hash=NULL, 
 pvalue_hash=NULL) {
  # TESTINDMVREG Conditional Independence Test for multivariate continous class variables 
  # PVALUE = TESTINDMVREG(Y, DATA, XINDEX, CSINDEX)
  # This test provides a p-value PVALUE for the NULL hypothesis H0 which is
  # X is independent by TARGET given CS. The pvalue is calculated following
  # nested models
  # This method requires the following inputs
  #   TARGET: a numeric matrix containing the values of the target (continuous) variable. 
  #   Its support can be R or any number betweeen 0 and 1, i.e. it contains proportions.
  #   DATASET: a numeric data matrix containing the variables for performing the test. They can be mixed variables. 
  #   XINDEX: the index of the variable whose association with the target we want to test. 
  #   CSINDEX: the indices if the variable to condition on. 
  # this method returns: the pvalue PVALUE, the statistic STAT.
  # References
  # [1] Norman R. Draper and Harry Smith. Applied Regression 
  # Analysis, Wiley, New York, USA, third edition, May 1998.
  # [2] Kanti V. Mardia, J. T. Kent and J. M. Bibby. Multivariate 
  # Analysis, Academic Press, New York, USA, 1979.
  #if the test cannot performed succesfully these are the returned values
  pvalue = log(1);
  stat = 0;
  csIndex[which(is.na(csIndex))] = 0;
  if ( hash )  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csIndex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if ( !is.null(stat_hash[key]) )   {
      stat = stat_hash[key];
      pvalue = pvalue_hash[key];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #if the xIndex is contained in csIndex, x does not bring any new information with respect to cs
  if ( !is.na(match(xIndex,csIndex)) )   {
    if ( hash ) {     #update hash objects
      stat_hash$key <- 0;#.set(stat_hash , key , 0)
      pvalue_hash$key <- log(1);#.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  
  #check input validity
  if( any(xIndex < 0) || any(csIndex < 0) ) {
    message(paste("error in testIndMVreg : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  
  xIndex = unique(xIndex);
  csIndex = unique(csIndex);
  x = dataset[ , xIndex];
  cs = dataset[ , csIndex];
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs) != 0 ) {
    if ( is.null(dim(cs)[2]) )  {     #cs is a vector
      if ( identical(x, cs) )  {    #if(!any(x == cs) == FALSE)
        if ( hash )  {    #update hash objects
          stat_hash$key <- 0;#.set(stat_hash , key , 0)
          pvalue_hash$key <- log(1);#.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for (col in 1:dim(cs)[2])  {
        if ( identical(x, cs[, col]) )  {   #if(!any(x == cs) == FALSE)
          if ( hash )  {    #update hash objects
            stat_hash$key <- 0;#.set(stat_hash , key , 0)
            pvalue_hash$key <- log(1);#.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
  #if the conditioning set (cs) is empty, we use a simplified formula
  if ( length(cs) == 0 )  {
    if ( !is.null(univariateModels) )  {
      pvalue = univariateModels$pvalue[[xIndex]];
      stat = univariateModels$stat[[xIndex]];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
    #compute the relationship between x,target directly
    fit2 = lm(target ~ x, weights = wei)
  } else  fit2 = lm(target ~., data = data.frame(dataset[ , c(csIndex, xIndex)] ), weights = wei )  
    if ( any( is.na( fit2$coefficients ) ) ) {
      stat <- 0
	    pvalue <- log(1)
	  } else {
	    fit1 = lm(target ~., data = data.frame( cs ), weights = wei )  
      mod = anova( fit1, fit2 ) 
      stat = mod[2, 5]     ## aproximate F test
      df1 = abs(mod[2, 2])
      df2 = mod[2, 1]
      pvalue = pf(stat, df1, df2, lower.tail= FALSE, log.p = TRUE)
    }
  #last error check
  if ( is.na(pvalue) | is.na(stat) ) {
    pvalue = log(1);
    stat = 0;
  } else {
    #update hash objects
    if( hash )  {
      stat_hash$key <- stat;         #.set(stat_hash , key , stat)
      pvalue_hash$key <- pvalue;           #.set(pvalue_hash , key , pvalue)
    }
  }
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
}