MMPC.gee <- function(target, reps = NULL, group, dataset, max_k = 3, threshold = 0.05 , test = NULL, ini = NULL, wei = NULL, user_test = NULL, 
                         hash=FALSE, hashObject=NULL, correl = "exchangeable", se = "jack", ncores = 1) {
  ##############################
  # initialization part of MMPC #
  ##############################
  runtime <- proc.time()
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
  ###################################
  # dataset checking and initialize #
  ###################################
  if ( is.null(dataset) || is.null(target) ) {
    stop('invalid dataset or target (class feature) arguments.');
  } else  target <- target;
  #check for NA values in the dataset and replace them with the variable mean
  if ( any( is.na(dataset) ) ) {
    warning("The dataset contains missing values and they were replaced automatically by the variable (column) median.")
    dataset <- apply(dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) });
  }
  ##################################
  # target checking and initialize #
  ##################################
  targetID <-  -1;
  #check if the target is a string
  if (is.character(target) & length(target) == 1) {
    findingTarget <- target == colnames(dataset);       #findingTarget <- target %in% colnames(dataset);
    if ( !sum(findingTarget)== 1 ) {
      warning('Target name not in colnames or it appears multiple times');
      return(NULL);
    }
    targetID <- which(findingTarget);
    target <- dataset[ , targetID];
  }
  
  #checking if target is a single number
  if (is.numeric(target) & length(target) == 1){
    if (target > dim(dataset)[2]){
      warning('Target index larger than the number of variables');
      return(NULL);
    }
    targetID <- target;
    target <- dataset[ , targetID];
  }
  ################################
  # test checking and initialize #
  ################################
  la <- length( unique( as.numeric(target) ) )
  
  if (typeof(user_test) == "closure") {
    test <- user_test;
  } else {
    #auto detect independence test in case of not defined by the user and the test is null or auto
    if ( is.null(test) || test == "auto" )  {
      if ( la > 2  &  sum(target - round(target) ) != 0 ) {
        test <- "testIndGEEReg"
      } else if (la == 2) {
        test <- "testIndGEELogistic"
      } else if ( la > 2  &  sum(target - round(target)) == 0 ) {
        test <- "testIndGEEPois"
      } 
    }  
    #available conditional independence tests
    av_tests <- c("testIndGEEReg", "testIndGEELogistic", "testIndGEEPois", "testIndGEEGamma", "auto",  NULL);
    ci_test <- test
    if (length(test) == 1) {   #avoid vectors, matrices etc
      test = match.arg(test, av_tests, TRUE);
      #convert to closure type
      if(test == "testIndGEEReg") {
        test <- testIndGEEReg;
      } else if(test == "testIndGEELogistic") {
        test <- testIndGEELogistic;
      } else if (test == "testIndGEEPois") {
        test <- testIndGEEPois;
      } else if (test == "testIndGEEGamma") {
        test <- testIndGEEGamma
      } else if (test == "testIndGEENormLog") {
        test <- testIndGEENormLog
      }  
      
    } else   stop('invalid test option');
  }
  ###################################
  # options checking and initialize #
  ###################################
  #extracting the parameters
  max_k <- floor(max_k);
  varsize <- dim(dataset)[2];
  #option checking
  if ( (typeof(max_k)!="double") || max_k < 1 )   stop('invalid max_k option');
  if ( max_k > varsize )   max_k = varsize;
  if ( (typeof(threshold) != "double") || threshold <= 0 || threshold > 1 )   stop('invalid threshold option');
  #######################################################################################
  if( !is.null(user_test) )  ci_test = "user_test";
  #call the main MMPC.gee function after the checks and the initializations
  results <- InternalMMPC.gee(target, reps, group, dataset, max_k, log(threshold), test, ini, wei, user_test, hash, varsize, stat_hash, 
                             pvalue_hash, targetID, correl = correl, se = se, ncores = ncores);
  
  runtime <- proc.time() - runtime
  MMPC.gee.output <- new("MMPC.gee.output", selectedVars = results$selectedVars, selectedVarsOrder = results$selectedVarsOrder, 
                             hashObject = results$hashObject, pvalues = results$pvalues, stats = results$stats, univ = results$univ, 
                             max_k = results$max_k, threshold = results$threshold, n.tests = results$n.tests, runtime = runtime, 
                             test = ci_test, correl = results$correl, se = results$se);
  return(MMPC.gee.output);
}


