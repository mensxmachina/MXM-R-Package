waldCR <- function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL){
  # Conditional independence test based on the Log Likelihood ratio test
  if ( !survival::is.Surv(target) )   stop('The survival test can not be performed without a Surv object target');
  csIndex[ which( is.na(csIndex) ) ] <- 0;
  
  if ( hash ) {
    csIndex2 <- csIndex[which(csIndex!=0)]
    csIndex2 <- sort(csIndex2)
    xcs <- c(xIndex,csIndex2)
    key <- paste(as.character(xcs) , collapse=" ");
    if ( is.null(stat_hash[key]) == FALSE ) {
      stat <- stat_hash[key];
      pvalue <- pvalue_hash[key];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #initialization: these values will be returned whether the test cannot be carried out
  pvalue <- log(1);
  stat <- 0;
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  cox_results <- NULL;
  cox_results_full <- NULL;
  event <- target[,2]
  x <- dataset[ , xIndex];
  numCases <- dim(dataset)[1];
  oop <- options(warn = -1) 
  on.exit( options(oop) )
  if ( length(event) == 0 )  event = vector('numeric',numCases) + 1;
      if ( length(csIndex) == 0 || sum(csIndex == 0, na.rm = TRUE) > 0 ) {
        cox_results <- survival::coxph(target ~ x, weights = wei )
        stat <- cox_results$wald.test
        pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE);
      } else {
        cox_results_full <- survival::coxph(target ~ ., data = as.data.frame(  dataset[ , c(csIndex, xIndex)] ), weights = wei) 
        res <- summary(cox_results_full)[[ 7 ]]
        pr <- dim(res)[1]
        stat <- res[pr, 4]^2
        pvalue < pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
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
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
}