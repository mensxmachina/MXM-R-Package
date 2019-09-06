permLLR <- function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL,
                  threshold = 0.05, R = 999){
  # Conditional independence test based on the Log Likelihood ratio test
  if (!survival::is.Surv(target) )   stop('The survival test can not be performed without a Surv object target');
  csIndex[which(is.na(csIndex))] = 0;
  thres <- threshold * R + 1
  if ( hash )  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csIndex2 = sort(csIndex2)
    xcs = c(xIndex, csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if (is.null(stat_hash[key]) == FALSE) {
      stat = stat_hash[key];
      pvalue = pvalue_hash[key];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #initialization: these values will be returned whether the test cannot be carried out
  pvalue <- log(1);
  stat <- 0;
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  llr_results = NULL;
  llr_results_full = NULL;   
  event <- target[,2]            
  #retrieving the data
  x <- dataset[ , xIndex];
  #if the censored indicator is empty, a dummy variable is created
  numCases <- dim(dataset)[1];
  if (length(event) == 0)  event = vector('numeric', numCases) + 1;
  if ( length(csIndex) == 0 || sum(csIndex == 0, na.rm = TRUE) > 0 ) {
    llr_results <- survival::survreg( target ~ x, weights = wei, control = list(iter.max = 5000) )
    stat <- 2 * abs( diff(llr_results$loglik) )
    if (stat > 0) {
      step <- 0
      j <- 1	
      n <- length(x)
      while (j <= R & step < thres ) {
        xb <- sample(x, n)  
        bit2 <-  survival::survreg( target ~ xb, weights = wei, control = list(iter.max = 5000), dist = "loglogistic" )
        stat2 <- 2 * abs( diff(bit2$loglik) )
        step <- step + ( stat2 > stat )
        j <- j + 1
      }
      pvalue <- log( (step + 1) / (R + 1) )        
    } else pvalue <- log(1)
    
  } else {
    llr_results <- survival::survreg( target ~ ., data = as.data.frame( dataset[ , csIndex] ), weights = wei, control = list(iter.max = 5000), dist = "loglogistic" ) 
    llr_results_full <- survival::survreg( target ~ ., data = as.data.frame(  dataset[ , c(csIndex, xIndex)] ), weights = wei, control = list(iter.max = 5000), dist = "loglogistic" )
    res <- anova(llr_results, llr_results_full)
    stat <- abs( res[2, 6] );
    if (stat > 0) {
      j <- 1	
      step <- 0
      n <- length(x)
      while (j <= R & step < thres ) {
        xb <- sample(x, n)  
        bit2 <- survival::survreg(target ~., data = as.data.frame( cbind(dataset[ ,csIndex], xb ) ), weights= wei, control = list(iter.max = 5000), dist = "loglogistic" )
        step <- step + ( anova(llr_results, bit2)[2, 6] > stat )
        j <- j + 1
      }
      pvalue <- (step + 1) / (R + 1)
    }  
  }  
  if ( is.na(pvalue) || is.na(stat) ) {
    pvalue <- log(1);
    stat <- 0;
  } else {
    #update hash objects
    if( hash )  {
      stat_hash[key] <- stat;      #.set(stat_hash , key , stat)
      pvalue_hash[key] <- pvalue;     #.set(pvalue_hash , key , pvalue)
    }
  }
  #testerrorcaseintrycatch(4);
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
}
