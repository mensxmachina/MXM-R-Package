MMPC.glmm <- function(target, reps = NULL, group, dataset, max_k = 3, threshold = 0.05 , test = NULL, ini = NULL, wei = NULL, user_test = NULL, 
                         hash=FALSE, hashObject=NULL, slopes = FALSE, ncores = 1) {
  ##############################
  # initialization part of MMPC #
  ##############################
  runtime <- proc.time()
  stat_hash <- NULL;
  pvalue_hash <- NULL;
  
  if ( hash )  {
    if (is.null(hashObject) )  {
      stat_hash <- Rfast::Hash();
      pvalue_hash <- Rfast::Hash();
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
      if ( survival::is.Surv(target) ) {
        test <- "testIndGLMMCR"
      } else if ( la > 2 & sum(target - round(target) ) != 0  &  is.null(wei)  &  !slopes  &  is.null(reps) ) {
        test <- "testIndLMM"
      } else if (la == 2) {
        test = "testIndGLMMLogistic"
      } else if ( la > 2  &  sum(target - round(target)) == 0 ) {
        test <- "testIndGLMMPois"
      } else if ( is.ordered(target) ) {
        test <- "testIndGLMMOrdinal"
      } else test <- "testIndGLMMReg"
    }  
    #available conditional independence tests
    av_tests <- c("testIndGLMMReg", "testIndGLMMLogistic", "testIndGLMMPois", "testIndGLMMGamma", 
                 "testIndGLMMNormLog", "testIndGLMMOrdinal", "testIndGLMMCR", "testIndLMM", "auto",  NULL);
    ci_test <- test
    if (length(test) == 1) {   #avoid vectors, matrices etc
      test <- match.arg(test, av_tests, TRUE);
      #convert to closure type
      if ( test == "testIndGLMMReg" ) {
        test <- testIndGLMMReg;
      } else if ( test == "testIndGLMMReg"  &  is.null(wei)  &  (!slopes) ) {
        test <- testIndLMM;
      } else if ( test == "testIndGLMMLogistic" ) {
          test <- testIndGLMMLogistic;
      } else if ( test == "testIndGLMMPois" ) {
        test <- testIndGLMMPois;
      } else if ( test == "testIndGLMMGamma" ) {
        test <- testIndGLMMGamma;
      } else if ( test == "testIndGLMMNormLog" ) {
        test <- testIndGLMMNormLog;
      } else if ( test == "testIndGLMMOrdinal" ) {
        test <- testIndGLMMOrdinal;
      } else if ( test == "testIndGLMMCR" ) {
        test <- testIndGLMMCR;
     } else if ( test == "testIndLMM" ) {
        test <- testIndLMM
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
  #call the main MMPC.glmm function after the checks and the initializations
  results <- InternalMMPC.glmm(target, reps, group, dataset, max_k, log(threshold), test, ini, wei, user_test, hash, varsize, 
                                   stat_hash, pvalue_hash, targetID, slopes = slopes, ncores = ncores);
  
  runtime <- proc.time() - runtime
  MMPC.glmm.output <-new("MMPC.glmm.output", selectedVars = results$selectedVars, selectedVarsOrder = results$selectedVarsOrder, 
                             hashObject=results$hashObject, pvalues=results$pvalues, stats=results$stats, univ = results$univ, 
                             max_k = results$max_k, threshold = results$threshold, n.tests = results$n.tests, runtime = runtime, 
                             test = ci_test, slope = slopes);
  return(MMPC.glmm.output);
}


