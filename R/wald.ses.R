wald.ses <- function(target, dataset, max_k = 3, threshold = 0.05, test = NULL, ini = NULL, wei = NULL, user_test = NULL, 
                    hash = FALSE, hashObject = NULL, ncores = 1, backward = FALSE) {
  ##############################
  # initialization part of SES #
  ##############################
  runtime <- proc.time()
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
  ###################################
  # dataset checking and initialize #
  ###################################
  if ( !is.null(dataset) ) {
    if( sum( class(target) == "matrix") == 1 ) {
      if( sum( class(target) == "Surv") == 1 )  stop('Invalid dataset class. For survival analysis provide a dataframe-class dataset');
    }
  }
  if( is.null(dataset) || is.null(target) ) {    #|| (dim(as.matrix(target))[2] != 1 & class(target) != "Surv" ))
    stop('invalid dataset or target (class feature) arguments.');
  } else  target <- target;
    if ( any( is.na(dataset) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
  }
  ##################################
  # target checking and initialize #
  ##################################
  targetID <-  -1;
  #check if the target is a string
  if (is.character(target) & length(target) == 1) {
    findingTarget <- target == colnames(dataset);#findingTarget <- target %in% colnames(dataset);
    if (!sum(findingTarget)==1){
      warning('Target name not in colnames or it appears multiple times');
      return(NULL);
    }
    targetID <- which(findingTarget);
    target <- dataset[ , targetID];
  }
  #checking if target is a single number
  if (is.numeric(target) & length(target) == 1) {
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
  la <- length( unique( as.numeric(target) ) )
  
  if(typeof(user_test) == "closure") {
    test <- user_test;
  } else {
    #auto detect independence test in case of not defined by the user and the test is null or auto
    if (is.null(test) || test == "auto") {
      
      if ( la == 2 )   target <- as.factor(target)
      #if target is a factor then use the Logistic test
      if ( is.factor(target) & la == 2  | la == 2 )  {
        test = "waldLogistic";
      } else if (is.factor(target) & la > 2) {
        test = "waldOrdinal";

      } else if ( ( is.numeric(target) || is.integer(target) ) & survival::is.Surv(target) == FALSE ) {
        
        if ( sum( floor(target) - target ) == 0  &  la > 2 )  test = "waldPois";
        
      } else if (survival::is.Surv(target) == TRUE) {
        test = "waldCR";
      } else   stop('Target must be a factor, vector or a Surv object');
    }
    #available conditional independence tests
    av_tests = c("waldBeta", "waldCR", "waldWR", "waldER", "waldLLR", "waldClogit", "waldLogistic", "waldPois", "waldNB", 
                 "waldBinom", "auto", "waldZIP", "waldMMReg", "waldIGreg", "waldOrdinal", "waldGamma", 
				 "waldNormLog", "waldTobit", "waldQPois", "waldQBinom", NULL);
    ci_test = test
    #cat(test)
    if (length(test) == 1) {      #avoid vectors, matrices etc
      test = match.arg(test, av_tests, TRUE);
      #convert to closure type
      if (test == "waldBeta") {
        test = waldBeta;
        
      } else if (test == "waldMMReg") {
        test = waldMMReg;
        
      } else if (test == "waldIGreg") {
        test = waldIGreg;
        
      } else if (test == "waldPois") { 
        test = waldPois;

      } else if (test == "waldNB") {
        test = waldNB;
        
      } else if (test == "waldGamma") {  
        test = waldGamma;
        
      } else if (test == "waldNormLog") {  
        test = waldNormLog;
        
      } else if (test == "waldZIP") {
        test = waldZIP;
        
      } else if (test == "waldTobit") { ## Tobit regression
        test = waldTobit;
        
      }  else if (test == "waldCR") {
        test = waldCR;
        
      } else if (test == "waldWR") {
        test = waldWR;
        
      } else if (test == "waldER") {
        test = waldER;
        
      } else if (test == "waldLLR") {
        test = waldLLR;
        
      } else if (test == "waldBinom") {
        test = waldBinom;
        
      } else if (test == "waldLogistic") {
        test = waldLogistic;
        
      } else if (test == "waldOrdinal") {
        test = waldOrdinal;
        
      } else if (test == "waldQPois") {
        test = waldQPois;
        
      } else if (test == "waldQBinom") {
        test = waldQBinom;
      }
      #more tests here
    } else {
      stop('invalid test option');
    }
  }
  ###################################
  # options checking and initialize #
  ###################################
  #extracting the parameters
  max_k = floor(max_k);
  varsize = ncol(dataset);
  #option checking
  if ( (typeof(max_k)!="double") || max_k < 1 )   stop('invalid max_k option');
  if ( max_k > varsize )    max_k = varsize;
  if ( (typeof(threshold) != "double" ) || threshold <= 0  || threshold > 1 )    stop('invalid threshold option');
  
  if ( !is.null(user_test) )   ci_test = "user_test";
  #call the main SES function after the checks and the initializations
  results <- wald.Internalses(target, dataset, max_k, log(threshold), test, ini, wei, user_test, hash, varsize, stat_hash, pvalue_hash, 
                             targetID, ncores = ncores);
  
  varsToIterate <- results$selectedVarsOrder
  
  if ( backward  & length( varsToIterate ) > 0  ) {
    varsOrder <- results$selectedVarsOrder
    bc <- mmpcbackphase(target, dataset[, varsToIterate, drop = FALSE], test = test, wei = wei, max_k = max_k, threshold = threshold)
    met <- bc$met
    results$selectedVars <- varsOrder[met]
    results$selectedVarsOrder <- varsOrder[met]
    results$signatures <- results$signatures[, met, drop = FALSE]
    results$pvalues[varsToIterate] <- bc$pvalues
    results$n.tests <- results$n.tests + bc$counter
  }
  
  runtime <- proc.time() - runtime
  SESoutput <-new("SESoutput", selectedVars = results$selectedVars, selectedVarsOrder = results$selectedVarsOrder, queues = results$queues, 
                  signatures = results$signatures, hashObject = results$hashObject, pvalues = results$pvalues, stats = results$stats, 
                  univ = results$univ, max_k = results$max_k, threshold = results$threshold, n.tests = results$n.tests, 
                  runtime = runtime, test=ci_test);
  return(SESoutput);
}










