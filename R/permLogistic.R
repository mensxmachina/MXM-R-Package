permLogistic = function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels = NULL, hash = FALSE, stat_hash = NULL, 
                           pvalue_hash = NULL, threshold = 0.05, R = 999) {

  csIndex[ which( is.na(csIndex) ) ] = 0
  thres <- threshold * R + 1
  
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
  pvalue = log(1)
  stat = 0;
  #information with respect to cs
  if ( !is.na( match(xIndex, csIndex) ) ) {
    if ( hash) {  #update hash objects
      stat_hash[key] <- 0;            #.set(stat_hash , key , 0)
      pvalue_hash[key] <- log(1);       #.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = 1, stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if ( any(xIndex < 0) || any(csIndex < 0) ) {
    message(paste("error in testIndLogistic : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  x = dataset[ , xIndex];
  cs = dataset[ , csIndex];
  if (length(cs) == 0 || any( is.na(cs) ) )  cs = NULL;
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs) != 0 )  {
    if ( is.null(dim(cs)[2]) ) {  #cs is a vector
      if ( identical(x, cs) ) {  #if(!any(x == cs) == FALSE)
        if (hash) {  #update hash objects
          stat_hash[key] <- 0;    #.set(stat_hash , key , 0)
          pvalue_hash[key] <- log(1);    #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for ( col in 1:dim(cs)[2] ) {
        if ( identical(x, cs[, col]) ) {  #if(!any(x == cs) == FALSE)
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
    fit2 = glm(target ~ x, binomial, weights = wei, model = FALSE)
    dev2 <- fit2$deviance
    stat = fit2$null.deviance - dev2
    if ( stat > 0 ) {
      step <- 0
      j <- 1		
      n <- length(target)
      while (j <= R & step < thres ) {
        xb <- sample(x, n)  
        bit2 <- glm(target ~ xb, binomial, weights = wei)  
        step <- step + ( bit2$deviance < dev2 )
        j <- j + 1
      }
      pvalue <- log( (step + 1) / (R + 1) )
     }  else pvalue <- log(1)	
  } else {
    fit1 <- glm(target ~ cs, binomial, weights = wei, model = FALSE)
    fit2 <- glm(target ~ cs + x, binomial, weights = wei, model = FALSE)
    dev2 <- fit2$deviance  
    stat <- fit1$deviance - dev2
    if (stat > 0) {
      step <- 0
      j <- 1		
      n <- length(target)
      while (j <= R & step < thres ) {
        xb <- sample(x, n)  
        bit2 <- glm(target ~ cs + xb, binomial, weights = wei, model = FALSE)  
        step <- step + ( bit2$deviance < dev2 )
        j <- j + 1
      }
      pvalue <- log( (step + 1) / (R + 1) )
    }  else  pvalue <- log(1)  
  }
  #update hash objects
  if (hash) {
    stat_hash[key] <- stat;      #.set(stat_hash , key , stat)
    pvalue_hash[key] <- pvalue;     #.set(pvalue_hash , key , pvalue)
  }
  #last error check
  if ( is.na(pvalue) || is.na(stat) ) {
    pvalue = log(1)
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