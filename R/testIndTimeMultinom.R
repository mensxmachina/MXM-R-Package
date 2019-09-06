testIndTimeMultinom = function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels = NULL, hash = FALSE, stat_hash = NULL, 
                               pvalue_hash = NULL) {
  
  csIndex[ which( is.na(csIndex) ) ] = 0
  if (hash) {
    csIndex2 = csIndex[which(csIndex!=0)]
    csIndex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if ( !is.null(stat_hash[key]) ) {
      stat = stat_hash[key];
      pvalue = pvalue_hash[key];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #if the test cannot performed succesfully these are the returned values
  pvalue = log(1);
  stat = 0;
  #information with respect to cs
  if ( !is.na( match(xIndex, csIndex) ) ) {
    if ( hash) {  #update hash objects
      stat_hash[key] <- 0;     #.set(stat_hash , key , 0)
      pvalue_hash[key] <- log(1);     #.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if ( any(xIndex < 0) || any(csIndex < 0) ) {
    message(paste("error in testIndLogistic : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  
  n <- 0.5 * dim(dataset)[1]
  x <- cbind(dataset[1:n, xIndex], dataset[(n + 1):(2 * n), xIndex]) 
  cs <- cbind(dataset[1:n, csIndex], dataset[(n + 1):(2 * n), csIndex]) 

  if ( length(cs) == 0 | any( is.na(cs) ) | sum( csIndex == 0 ) > 0 )  cs = NULL;
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs) != 0 )  {
    if ( is.null(dim(cs)[2]) ) {  #cs is a vector
      if ( identical(x, cs) ) {  # if(!any(x == cs) == FALSE)
        if (hash) {     #update hash objects
          stat_hash[key] <- 0;    #.set(stat_hash , key , 0)
          pvalue_hash[key] <- log(1);    #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for ( col in 1:dim(cs)[2] ) {
        if ( identical( x, cs[, col] ) ) {  #if(!any(x == cs) == FALSE)
          if (hash) {    #update hash objects
            stat_hash[key] <- 0;   #.set(stat_hash , key , 0)
            pvalue_hash[key] <- log(1);   #.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
  
  if (length(cs) == 0) {
    #if the univariate models have been already computed
    if ( !is.null(univariateModels) ) {
      pvalue = univariateModels$pvalue[[xIndex]];
      stat = univariateModels$stat[[xIndex]];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
    fit1 <- nnet::multinom(target ~ 1, weights = wei, trace = FALSE)
    fit2 <- nnet::multinom(target ~ x, weights = wei, trace = FALSE)
  } else {
    fit1 <- nnet::multinom(target ~ cs, weights = wei, trace = FALSE)
    fit2 <- nnet::multinom(target ~ cs + x, weights = wei, trace = FALSE)
  }
  stat <- fit1$deviance - fit2$deviance
  dof <- length( coef(fit2) ) - length( coef(fit1) )
  pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE); 
  #update hash objects
  if (hash) {
    stat_hash[key] <- stat;      #.set(stat_hash , key , stat)
    pvalue_hash[key] <- pvalue;     #.set(pvalue_hash , key , pvalue)
  }
  #last error check
  if ( is.na(pvalue) || is.na(stat) ) {
    pvalue = log(1);
    stat = 0;
  } else {
    #update hash objects
    if (hash) {
      stat_hash[key] <- stat;      #.set(stat_hash , key , stat)
      pvalue_hash[key] <- pvalue;      #.set(pvalue_hash , key , pvalue)
    }
  }
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
}