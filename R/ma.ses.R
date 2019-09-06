# R Implementation of the Statistical equivalence signatures algorithm (SES)
# as described in the "Discovering multiple, equivalent biomarker signatures" paper
# by Ioannis Tsamardinos, Vincenzo Lagani and Dimitris Pappas
# R Implementation by Giorgos Athineou (2013-2014)
# VERSION: 17/3/2014
# INPUTS
# target : the class variable , provide either a vector, an 1D matrix, an 1D array (same length with the data rows), a factor 
# or a formula. The data can be either continuous data in R, values within (0,1), binary [0,1], nominal or ordinal data.
# dataset : tha dataset , provide a data frame (columns = variables , rows = samples) or a matrix or an ExpressionSet.
# max_k : the maximum conditioning set which is used in the current conditional indepedence test used.
# threshold : the significance threshold ( must be in (0,1) ) for testing the null hypothesis of the generated pvalues.
# test : the conditional independence test we are going to use. the available conditional independence tests so far in this 
# implementation are:
#   "maFisher" : Fisher conditional independence test for continous targets (or proportions) and continuous predictors only
#   "testIndSpearman" : Fisher conditional independence test for continous targets (or proportions) and continuous predictors only (Spearman correlation is calculated first)
#   "testIndReg" : Conditional independence test based on regression for continous targets (or proportions) and mixed predictors using the F test
#   "testIndRQ" : Conditional Independence Test based on quantile (median) regression for numerical class variables and mixed predictors (F test)
#   "testIndLogistic" : Conditional Independence Test based on logistic regression for binary,categorical or ordinal class variables and mixed predictors
#   "testIndPois" : Conditional Independence Test based on Poisson regression for discrete class variables and mixed predictors (log-likelihood ratio test)
#   "testIndZIP" : Conditional Independence Test based on zero inflated poisson regression for discrete class variables and mixed predictors (log-likelihood ratio test)
#   "testIndNB" : Conditional Independence Test based on negative binomial regression for discrete class variables and mixed predictors (log-likelihood ratio test)
#   "testIndBeta" : Conditional Independence Test based on beta regression for proportions and mixed predictors (log likelihood ratio test)
#   "testIndMVreg" : Conditional Independence Test based on mu;ltivariate linear regression for Euclidean data and mixed predictors (log likelihood ratio test)
#   "gSquare" : Conditional Independence test based on the G test of independence (log likelihood ratio  test)
#   "censIndCR" : Conditional independence test for survival data based on the Log likelihood ratio test with mixed predictors (Cox regression)
# user_test : the user defined conditional independence test ( provide a closure type object )
# hash : a boolean variable whuch indicates whether (TRUE) or not (FALSE) to use the hash-based implementation of the statistics of SES.
# hashObject : a List with the hash objects (hash package) on which we have cached the generated statistics. 
#              SES requires this Object for the hash-based implementation of the statistics. This hashObject is produced or 
# updated by each run of SES (if hash == TRUE) and it can be reused in next runs of SES.
# there are default values for all of the parameters of the algorithm.
# OUTPUT <LIST>
# The output of the algorithm is a LIST with the following quantities (14) :
# selectedVars : the selected variables i.e. the dependent of the target variables.
# selectedVarsOrder : the increasing order of the selected variables due to their pvalues
# queues : the selected variable queues with the multiple statistically equivalent variables 
# (if you want to make multiple statistically equivalent signatures you 
# have to take one variable from each queue).
# signatures : all the possible combinations of the variables in the queues. One variable per queue in each signature. 
# (signature ~ the minimum subset with the most relevant features).
# hashObject : the hashObject with the cached statistic results of the current run.
# pvalues : the pvalues of all of the variables.
# stats : the stats of all of the variables.
# all_queues : the queues of all of the variables.
# data : the dataset used in the current run.
# target : the class variable used in the current run.
# test : the conditional independence test used in the current run.
# max_k : the max_k option used in the current run.
# threshold : the threshold option used in the current run.
# runtime : the run time of the algorithm.
# Conditional independence test arguments have to be in this exact fixed order : 
# target(target variable), data(dataset), xIndex(x index), csIndex(cs index),
# univariateModels(cached statistics for the univariate indepence test), hash(hash booleab), stat_hash(hash object), 
# output of each test: LIST of the generated pvalue, stat, and the updated hash objects.
#   if equal_case = 1 then, if we have more than one equivalent vars in z , we select the one with the most closer pvalue to the pvalue of cvar
#   if equal_case = 2 then, if we have more than one equivalent vars in z , we select the one with the most minimum pvalue (>a)
#   else in any other case, if we have more than one equivalent vars in z , we select the first one
# In this version we support the equal_case = 3.
# #hashObject
ma.ses <- function(target, dataset, ina, statistic = FALSE, max_k = 3, threshold = 0.05, test = NULL , ini = NULL, user_test = NULL, 
                  hash=FALSE, hashObject=NULL, ncores = 1) {
  ##############################
  # initialization part of SES #
  ##############################
  stat_hash <- NULL
  pvalue_hash <- NULL
  
  if ( hash )  {
    if (is.null(hashObject) )  {
      stat_hash <- Rfast::Hash()
      pvalue_hash <- Rfast::Hash()
    } else if ( class(hashObject) == "list" ) {
      stat_hash <- hashObject$stat_hash;
      pvalue_hash <- hashObject$pvalue_hash;
    } else   stop('hashObject must be a list of two hash objects (stat_hash, pvalue_hash)')
  }
  ##################################
  # target checking and initialize #
  ##################################
  targetID = -1;
  #check if the target is a string
  if (is.character(target) & length(target) == 1) {
    findingTarget <- target == colnames(dataset);      #findingTarget <- target %in% colnames(dataset);
    if ( !sum(findingTarget) == 1 ) {
      warning('Target name not in colnames or it appears multiple times');
      return(NULL);
    }
    targetID <- which(findingTarget);
    target <- dataset[ , targetID];
  }
  
  #checking if target is a single number
  if ( is.numeric(target) & length(target) == 1 ) {
    if (target > dim(dataset)[2]) {
      warning('Target index larger than the number of variables');
      return(NULL);
    }
    targetID <- target;
    target <- dataset[ , targetID];
  }
  ################################
  # test checking and initialize #
  ################################
  #cat("\nConditional independence test used: ");cat(test);cat("\n");
    #available conditional independence tests
    av_tests = c("testIndFisher", "testIndSpearman", "testIndReg", "testIndRQ", "testIndBeta", "censIndCR","censIndWR", 
	             "testIndClogit", "testIndLogistic", "testIndPois", "testIndNB", "testIndBinom", "gSquare", "auto" , 
				 "testIndZIP" , "testIndMVreg", "testIndMMReg", NULL);
    ci_test = test
    #cat(test)
    if(length(test) == 1)  {#avoid vectors, matrices etc

      test = match.arg(test, av_tests, TRUE);
      #convert to closure type
      if (test == "testIndFisher") {
        if ( is.data.frame(dataset) ) {
          if ( length( Rfast::which.is(dataset) ) > 0 ) {
            warning("Dataset contains categorical variables (factors). A regression model is advised to be used instead.")
          } else dataset <- as.matrix(dataset)
        }
        test = testIndFisher;
      }
      else if (test == "testIndSpearman") {
        if ( is.data.frame(dataset) ) {
          if ( length( Rfast::which.is(dataset) ) > 0 ) {
            warning("Dataset contains categorical variables (factors). A regression model is advised to be used instead.")
          } else dataset <- as.matrix(dataset)
        }
        target = rank(target)
        dataset = apply(dataset, 2, rank)  
        test = testIndSpearman;  ## Spearman is Pearson on the ranks of the data
      }
	    
    } else   stop('invalid test option');
  ###################################
  # options checking and initialize #
  ###################################
  #extracting the parameters
  max_k = floor(max_k);
  varsize = ncol(dataset);
  #option checking
  if ( (typeof(max_k) != "double" ) || max_k < 1 )   stop('invalid max_k option');
  if ( max_k > varsize )   max_k = varsize;
  if ( (typeof(threshold) != "double" ) || threshold <= 0 || threshold > 1 )   stop('invalid threshold option');
  #######################################################################################
  if ( !is.null(user_test) )  ci_test = "user_test";
  ## end of checking
  D = max(ina)
  targ = list()
  data = list()
  for (l in 1:D) {  
   targ[[ l ]] = target[ ina == l ]
   data[[ l ]] = dataset[ ina == l, ]
  }
  ####### 
  for (i in 1:D) {
    if ( is.null(data[[ l ]] ) || is.null(targ[[ l ]]) ) {  #|| (dim(as.matrix(target))[2] != 1 & class(target) != "Surv" ))
      stop('invalid dataset or target (class feature) arguments.');
    } else {
      targ[[ l ]] = targ[[ l ]];
    }
    #check for NA values in the dataset and replace them with the variable median or the mode
   da <- data[[ l ]]	 
   if ( any( is.na( da ) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix( da ) )  {
      da <- apply( da, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- unique( which( is.na(da), arr.ind = TRUE )[, 2] )
      for ( i in poia )  {
        xi <- da[, i]
        if ( is.numeric(xi) )  {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) ) {
          xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        }
        da[, i] <- xi
      }
    }
   }
   data[[ l ]] <- da
  }
  #call the main SES function after the checks and the initializations
  results = Internalmases(targ, data, statistic, max_k, log(threshold), test, ini, user_test, hash, varsize, stat_hash, pvalue_hash, targetID, ncores = ncores);
  masesoutput <-new("mases.output", selectedVars = results$selectedVars, selectedVarsOrder=results$selectedVarsOrder, queues=results$queues, signatures=results$signatures, hashObject=results$hashObject, pvalues=results$pvalues, stats=results$stats, univ = results$univ, max_k=results$max_k, threshold = results$threshold, runtime=results$runtime, test=ci_test);
  
  return(masesoutput);
}






