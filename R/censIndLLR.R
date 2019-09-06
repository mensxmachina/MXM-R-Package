censIndLLR <- function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL){
  # Conditional independence test based on the Log Likelihood ratio test
  if (!survival::is.Surv(target) )   stop('The survival test can not be performed without a Surv object target');
  csIndex[which(is.na(csIndex))] = 0;
  if ( hash )  {
    csIndex2 <- csIndex[which(csIndex!=0)]
    csindex2 <- sort(csIndex2)
    xcs <- c(xIndex,csIndex2)
    key <- paste(as.character(xcs) , collapse=" ");
    if (is.null(stat_hash[[key]]) == FALSE) {
      stat <- stat_hash[[key]];
      pvalue <- pvalue_hash[[key]];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #initialization: these values will be returned whether the test cannot be carried out
  pvalue <- log(1);
  stat <- 0;
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  weibull_results <- NULL;
  weibull_results_full <- NULL;
  #timeIndex = dim(dataset)[2];
  event <- target[,2]
  x <- dataset[ , xIndex];
  #if the censored indicator is empty, a dummy variable is created
  numCases <- dim(dataset)[1];
  if (length(event) == 0)  event = vector('numeric',numCases) + 1;
      if ( length(csIndex) == 0 || sum(csIndex == 0, na.rm = TRUE) > 0 ) {
        llr_results <- try( survival::survreg( target ~ x, weights = wei, control = list(iter.max = 5000), dist = "loglogistic" ), silent = TRUE )
        if ( identical( class(llr_results), "try-error" ) ) {
          stat <- NA  
        } else {  
          dof <- length( coef(llr_results) ) - 1
          stat <- 2 * abs( diff(llr_results$loglik) )
          pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE);
        }  
      } else {
        llr_results <- try( survival::survreg(target ~ ., data = as.data.frame( dataset[ , c(csIndex)] ), weights = wei, control = list(iter.max = 5000), dist = "loglogistic" ), silent = TRUE ) 
        llr_results_full <- try( survival::survreg(target ~ ., data = as.data.frame(  dataset[ , c(csIndex, xIndex)] ), weights = wei, control = list(iter.max = 5000), dist = "loglogistic" ), silent = TRUE)
        if ( identical( class(llr_results), "try-error" )  |  identical( class(llr_results_full), "try-error" ) ) {
          stat <- NA
        } else {  
          res <- anova(llr_results, llr_results_full)
          stat <- abs( res[2, 6] );
          dof <- abs( res[2, 5] );
          pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        }  
      }  
      if ( is.na(pvalue) || is.na(stat) ) {
        pvalue <- log(1);
        stat <- 0;
      } else {
        #update hash objects
        if( hash )  {
          stat_hash[[key]] <- stat;      #.set(stat_hash , key , stat)
          pvalue_hash[[key]] <- pvalue;     #.set(pvalue_hash , key , pvalue)
        }
      }
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
}
