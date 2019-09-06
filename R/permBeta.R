permBeta = function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL,
                       threshold = 0.05, R = 999) {
  #   TESTINDBETA Conditional Independence Test based on beta regression for proportions
  #   provides a p-value PVALUE for the null hypothesis: X independent by target
  #   given CS. The pvalue is calculated by comparing a beta regression model based 
  #   on the conditioning set CS against a model containing both X and CS. 
  #   The comparison is performed through a chi-square test with some degrees 
  #   of freedom on the difference between the log-likelihoodss of the two models. 
  #   TESTINDBETA requires the following inputs:
  #   target: a column vector containing the values of the target variable. 
  #   target must be an integer vector, with values between 0 and 1 
  #   dataset: a numeric data matrix containing the variables for performing
  #   the conditional independence test. There can be mixed variables, i.e. continous and or categorical
  #   xIndex: the index of the variable whose association with the target
  #   must be tested. Can be mixed variables, either continous or categorical
  #   csIndex: the indices of the variables to condition on. 
  csIndex[which(is.na(csIndex))] = 0
  thres <- threshold * R + 1
  
  if ( hash )  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csIndex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if ( !is.null(stat_hash[key]) )  {
      stat = stat_hash[key];
      pvalue = pvalue_hash[key];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  pvalue = log(1)
  stat = 0;
  #information with respect to cs
  if ( !is.na(match(xIndex, csIndex)) ) {
    if ( hash )  {   #update hash objects
      stat_hash[key] <- 0;     #.set(stat_hash , key , 0)
      pvalue_hash[key] <- 1;        #.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = 1, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if ( any(xIndex < 0) || any(csIndex < 0) )  {
    message(paste("error in testIndBeta : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  xIndex = unique(xIndex);
  csIndex = unique(csIndex);
  x = dataset[ , xIndex];
  cs = dataset[ , csIndex];
  if ( length(cs) == 0 || any( is.na(cs) ) )  cs = NULL;
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs) != 0 )  {
    if ( is.null(dim(cs)[2]) )  {    #cs is a vector
      if ( identical(x, cs) )  {      #if(!any(x == cs) == FALSE)
        if ( hash )  {    #update hash objects
          stat_hash[key] <- 0;     #.set(stat_hash , key , 0)
          pvalue_hash[key] <- 1;     #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = 1, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for (col in 1:dim(cs)[2]) {
        if ( identical(x, cs[, col]) )  {    #if(!any(x == cs) == FALSE)
          if ( hash )  {    #update hash objects
            stat_hash[key] <- 0;    #.set(stat_hash , key , 0)
            pvalue_hash[key] <- 1;    #.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = 1, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
  #trycatch for dealing with errors
  res <- tryCatch(
{ 
  n <- length(target)
  #if the conditioning set (cs) is empty, we use the t-test on the coefficient of x.
  if ( length(cs) == 0 )  {
    #if the univariate models have been already compute
    if ( !is.null(univariateModels) )   {
      pvalue = univariateModels$pvalue[[xIndex]];
      stat = univariateModels$stat[[xIndex]];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
    #Fitting beta regressions
    if ( is.null(wei) )  fit1 = Rfast::beta.mle(target)  else  fit1 <- betamle.wei(target, wei)
    fit2 <- beta.reg(target, x, wei = wei)
    lik2 <- fit2$loglik
    stat <- 2 * lik2 - 2 * fit1$loglik
    step <- 0
    j <- 1		
    while (j <= R & step < thres ) {
      xb <- x[sample(n, n), ]  
      bit2 <- beta.reg(target, xb, wei = wei )  ;
      step <- step + ( bit2$loglik > lik2 )
      j <- j + 1
    }
    pvalue <-  log( (step + 1) / (R + 1) )
  } else {
    #Fitting beta regressions
    fit1 <- beta.reg(target, cs, wei = wei)
    fit2 <- beta.reg(target, cbind(cs, x) , wei = wei )  ;
  }
  lik2 <- fit2$loglik
  stat <- 2 * lik2 - 2 * fit1$loglik
  step <- 0
  j <- 1		
  while (j <= R & step < thres ) {
    xb <- x[sample(n, n), ]  
    bit2 <- beta.reg(target, cbind(cs, xb), wei = wei )  ;
    step <- step + ( bit2$loglik > lik2 )
    j <- j + 1
  }
  pvalue <-  log( (step + 1) / (R + 1) )

  #update hash objects
  if ( hash )  {
    stat_hash[key] <- stat;    #.set(stat_hash , key , stat)
    pvalue_hash[key] <- pvalue;       #.set(pvalue_hash , key , pvalue)
  }
  #last error check
  if (is.na(pvalue) || is.na(stat) ) {
    pvalue = log(1)
    stat = 0;
  } else {
    #update hash objects
    if ( hash )  {
      stat_hash[key] <- stat;      #.set(stat_hash , key , stat)
      pvalue_hash[key] <- pvalue;        #.set(pvalue_hash , key , pvalue)
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