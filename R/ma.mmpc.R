ma.mmpc <- function(target, dataset, ina, statistic = FALSE, max_k = 3, threshold = 0.05, test = NULL, ini = NULL, user_test = NULL, 
                   hash=FALSE, hashObject=NULL, ncores = 1, backward = FALSE) {
  ##############################
  # initialization part of MMPC #
  ##############################
  stat_hash <- NULL
  pvalue_hash <- NULL
  
  if ( hash )  {
    if (is.null(hashObject) )  {
      stat_hash <- Rfast::Hash();
      pvalue_hash <- Rfast::Hash();
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
    findingTarget <- target == colnames(dataset);#findingTarget <- target %in% colnames(dataset);
    if (!sum(findingTarget)==1 ) {
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
      if(test == "testIndFisher") {

        if ( is.data.frame(dataset) ) {
          if ( length( Rfast::which.is(dataset) ) > 0 ) {
            warning("Dataset contains categorical variables (factors). A regression model is advised to be used instead.")
          } else dataset <- as.matrix(dataset)
        }
            
        test = testIndFisher;
      }
      else if(test == "testIndSpearman") {

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
  if ( (typeof(max_k) != "double" ) || max_k < 1 )   stop('invalid max_k option');
  if ( max_k > varsize )   max_k = varsize;
  if ( (typeof(threshold) != "double" ) || threshold <= 0 || threshold > 1 )   stop('invalid threshold option');
  #######################################################################################
  if(!is.null(user_test))  ci_test = "user_test";

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
  #call the main MMPC function after the checks and the initializations
  results = Internalmammpc(targ, data, statistic, max_k, log(threshold), test, ini, user_test, hash, varsize, stat_hash, pvalue_hash, targetID, ncores = ncores);
  mammpcoutput <-new("mammpc.output", selectedVars = results$selectedVars, selectedVarsOrder=results$selectedVarsOrder, hashObject=results$hashObject, pvalues=results$pvalues, stats=results$stats, univ=results$univ, max_k=results$max_k, threshold = results$threshold, runtime=results$runtime, test=ci_test); 
  return(mammpcoutput);
}






#########################################################################################################














