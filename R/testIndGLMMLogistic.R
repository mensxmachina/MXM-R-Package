testIndGLMMLogistic = function(target, reps = NULL, group, dataset, xIndex, csIndex,  wei =  NULL, univariateModels = NULL,
 hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, slopes = FALSE) {
  #   TESTINDGLMM Conditional Independence Test based on generalised linear mixed models for normal, binary ro discrete variables
  #   target: a vector containing the values of the target variable. 
  #   target must be a vector with percentages, binay data, numerical values or integers
  #   reps: a vector with the time points (if available)
  #   group: a vector indicating the groupings of the subjects.       
  #   dataset: a numeric data matrix containing the variables for performing
  #   the conditional independence test. They can be mixed variables, either continous or categorical
  #   xIndex: the index of the variable whose association with the target
  #   must be tested. Can be any type of variable, either continous or categorical.
  #   csIndex: the indices of the variables to condition on. They can be mixed variables, either continous or categorical
  #   this method returns: the pvalue PVALUE, the statistic STAT.
  #cast factor into numeric vector
  csIndex[which(is.na(csIndex))] = 0
  
  if( hash )  {
    csIndex2 <- csIndex[which(csIndex!=0)]
    csIndex2 <- sort(csIndex2)
    xcs <- c(xIndex,csIndex2)
    key <- paste(as.character(xcs) , collapse=" ");
    if( !is.null(stat_hash[key]) )  {
      stat <- stat_hash[key];
      pvalue <- pvalue_hash[key];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #if the test cannot performed succesfully these are the returned values
  pvalue <- log(1);
  stat <- 0; 
  #information with respect to cs
  if ( !is.na(match(xIndex, csIndex)) )  {
    if ( hash )  {      #update hash objects
      stat_hash$key <- 0;#.set(stat_hash , key , 0)
      pvalue_hash$key <- log(1);#.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if ( any(xIndex < 0) || any(csIndex < 0) ) {
    message(paste("error in testIndGLMM : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #extract the data
  x <- dataset[ , xIndex];
  cs <- dataset[ , csIndex];
  if ( length(cs) == 0 || any( is.na(cs) ) )  cs <- NULL;
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  oop <- options(warn = -1) 
  on.exit( options(oop) )
  if ( length(cs) != 0 )  {
    if ( is.null(dim(cs)[2]) )  {     #cs is a vector
      if ( identical(x, cs) )  {    #if(!any(x == cs) == FALSE)
        if ( hash )  {    #update hash objects
          stat_hash$key <- 0;           #.set(stat_hash , key , 0)
          pvalue_hash$key <- log(1);           #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for (col in 1:dim(cs)[2])  {
        if ( identical(x, cs[, col]) )  {   #if(!any(x == cs) == FALSE)
          if ( hash )  {    #update hash objects
            stat_hash$key <- 0;         #.set(stat_hash , key , 0)
            pvalue_hash$key <- log(1);         #.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
 #if the conditioning set (cs) is empty, we use the t-test on the coefficient of x.
  if (length(cs) == 0)  {
    #if the univariate models have been already compute
    if ( !is.null(univariateModels) )  {
      pvalue <- univariateModels$pvalue[[xIndex]];
      stat <- univariateModels$stat[[xIndex]];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
   if ( is.null(reps) ) {
     fit1 <- lme4::glmer( target ~ (1|group), weights = wei, family = binomial ) 
     fit2 <- lme4::glmer( target ~ (1|group) + x, weights = wei, family = binomial ) 
   } else {
      reps <- reps 
      if ( slopes ) {
        fit1 <- lme4::glmer( target ~ reps + (reps|group), weights = wei, family = binomial ) 
        fit2 <- lme4::glmer( target ~ reps + (reps|group) + x, weights = wei, family = binomial ) 
      } else {
        reps <- reps 
        fit1 <- lme4::glmer( target ~ reps + (1|group), weights = wei, family = binomial )  
        fit2 <- lme4::glmer( target ~ reps + (1|group) + x, weights = wei, family = binomial )  
      }
   }
   
  } else {
    if ( is.null(reps) ) {
      fit1 <- lme4::glmer( target ~ (1|group) + cs, weights = wei, family = binomial )    
      fit2 <- lme4::glmer( target ~ (1|group) + cs + x, weights = wei, family = binomial )    
    } else {
      reps <- reps 
      if (slopes ) {
        fit1 <- lme4::glmer( target ~ reps + (reps|group) + x, weights = wei, family = binomial )
        fit2 <- lme4::glmer( target ~ reps + (reps|group) + cs + x, weights = wei, family = binomial )
      } else {
        fit1 <- lme4::glmer( target ~ reps + (1|group) + cs, weights = wei, family = binomial) 
        fit2 <- lme4::glmer( target ~ reps + (1|group) + cs + x, weights = wei, family = binomial ) 
      }
    }  
  }
  #calculate the p value and stat.
  mod <- anova(fit1, fit2)
  stat <- mod[2, 6]
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
  #update hash objects
  if ( hash )  {
    stat_hash$key <- stat;   #.set(stat_hash, key, stat)
    pvalue_hash$key <- pvalue;    #.set(pvalue_hash, key, pvalue)
  }
  results <- list(pvalue = pvalue, stat = stat, stat_hash = stat_hash, pvalue_hash = pvalue_hash);
  return(results);
}