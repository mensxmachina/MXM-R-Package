testIndClogit = function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL){
  # Conditional independence test based on the Log Likelihood ratio test
  if ( !is.matrix(target) || ( is.matrix(target)  &  ncol(target) != 2 ) )   stop('The testIndClogit test can not be performed without a 2 column matrix target');
  csIndex[which(is.na(csIndex))] = 0;
  
  if( hash )  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csIndex2 = sort(csIndex2)
    xcs = c(xIndex, csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if( !is.null(stat_hash[key]) )  {
      stat = stat_hash[key];
      pvalue = pvalue_hash[key];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #initialization: these values will be returned whether the test cannot be carried out
  pvalue = log(1);
  stat = 0;
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  clogit_results = NULL;
  clogit_results_full = NULL;
  id = target[, 2] #the patient id
  x = dataset[ , xIndex];
  case = target[, 1]  ## case control, 1 is the case, 0 is the control

  res <- tryCatch(
    {
      if ( length(csIndex) == 0 || sum(csIndex == 0, na.rm = TRUE) > 0 ) {
        clogit_results <- survival::clogit(case ~ x + strata(id) )
        dof = length( coef(clogit_results) ) 
        stat = 2 * abs( diff(clogit_results$loglik) )
        pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE);
      } else {
        clogit_results <- survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , csIndex] ) )
        clogit_results_full <- survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , c(csIndex, xIndex)] ) )
        res = anova(clogit_results_full, clogit_results)
        stat = res[2, 2]
        dF = res[2, 3]
        pvalue = pchisq(stat, dF, lower.tail = FALSE, log.p = TRUE)
      }  
      if ( is.na(pvalue) || is.na(stat) ) {
        pvalue = log(1);
        stat = 0;
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
      
    },
    error=function(cond) {
      #   message(paste("error in try catch of the testIndZIP test"))
      #   message("Here's the original error message:")
      #   message(cond)
      #   #        #for debug
      #   #        print("\nxIndex = \n");
      #   #        print(xIndex);
      #   #        print("\ncsindex = \n");
      #   #        print(csIndex);
      #   stop();
      #error case
      pvalue = log(1);
      stat = 0;
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    },
    finally={}
  )
  return(res);
}
  