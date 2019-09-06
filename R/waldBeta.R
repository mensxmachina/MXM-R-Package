waldBeta <- function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels=NULL, 
                     hash = FALSE, stat_hash=NULL, pvalue_hash=NULL){
  #initialization
  csIndex[ which( is.na(csIndex) ) ] = 0
  
  if ( hash ) {
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
    if ( hash ) {  #update hash objects
      stat_hash[key] <- 0;    #.set(stat_hash , key , 0)
      pvalue_hash[key] <- log(1);     #.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if ( any(xIndex < 0) || any(csIndex < 0) ) {
    message(paste("error in testIndBeta : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #extract the data
  x <- dataset[ , xIndex];
  cs <- dataset[ , csIndex];
  if ( length(cs) == 0 || any( is.na(cs) ) )  cs <- NULL;
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs) != 0 )   {
    if ( is.null(dim(cs)[2]) ) {  #cs is a vector
      if (any(x != cs) == FALSE) { #if(!any(x == cs) == FALSE)
        if ( hash ) {     #update hash objects
          stat_hash[key] <- 0;          #.set(stat_hash , key , 0)
          pvalue_hash[key] <- log(1);          #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash = pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for (col in 1:dim(cs)[2]) {
        if (any(x != cs[,col]) == FALSE) { #if(!any(x == cs) == FALSE)
          if ( hash ) {      #update hash objects
            stat_hash[key] <- 0;          #.set(stat_hash , key , 0)
            pvalue_hash[key] <- log(1);         #.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash = pvalue_hash);
          return(results);
        }
      }
    }
  }
  #trycatch for dealing with errors
  res <- tryCatch(
    {
      #if the conditioning set (cs) is empty, we use the t-test on the coefficient of x
      if (length(cs) == 0) {
        #Fitting beta regression
        fit <- beta.mod(target, x, wei = wei)
      } else  fit <- beta.mod(target, dataset[, c(csIndex, xIndex)], wei = wei )
      res <- fit$be
      pr <- dim(res)[1]
      stat <- res[pr, 3]
      pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE) 
      #update hash objects
      if ( hash ) {
        stat_hash[key] <- stat;#.set(stat_hash , key , stat)
        pvalue_hash[key] <- pvalue;#.set(pvalue_hash , key , pvalue)
      }
      #last error check
      if ( is.na(pvalue) || is.na(stat) ) {
        pvalue <- log(1);
        stat <- 0;
      } else {
        #update hash objects
        if( hash ) {
          stat_hash[key] <- stat   # .set(stat_hash , key , stat)
          pvalue_hash[key] <- pvalue   # .set(pvalue_hash , key , pvalue)
        }
      }
      #testerrorcaseintrycatch(4);
      results <- list(pvalue = pvalue, stat = stat,  stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    },
    error = function(cond) {
      pvalue <- log(1);
      stat <- 0;
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    },
    #   warning=function(cond) {
    #     #do nothing, or
    #     message(paste("Warning in the testIndBeta testL"))
    #     message("Here's the original warning message:")
    #     message(cond)
    #   },
    finally = {}
  )
  
  return(res);
}