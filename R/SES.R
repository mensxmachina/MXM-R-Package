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
#   "testIndFisher" : Fisher conditional independence test for continous targets (or proportions) and continuous predictors only
#   "testIndReg" : Conditional independence test based on regression for continous targets (or proportions) and mixed predictors using the F test
#   "testIndRQ" : Conditional Independence Test based on quantile (median) regression for numerical class variables and mixed predictors (F test)
#   "testIndLogistic" : Conditional Independence Test based on logistic regression for binary,categorical or ordinal class variables and mixed predictors
#   "testIndPois" : Conditional Independence Test based on Poisson regression for discrete class variables and mixed predictors (log-likelihood ratio test)
#   "testIndNB" : Conditional Independence Test based on Negative Binomial regression for discrete class variables and mixed predictors (log-likelihood ratio test)
#   "testIndBeta" : Conditional Independence Test based on beta regression for proportions and mixed predictors (log likelihood ratio test)
#   "gSquare" : Conditional Independence test based on the G test of independence (log likelihood ratio  test)
#   "censIndLR" : Conditional independence test for survival data based on the Log likelihood ratio test with mixed predictors
# user_test : the user defined conditional independence test ( provide a closure type object )
# hash : a boolean variable whuch indicates whether (TRUE) or not (FALSE) to use the hash-based implementation of the statistics of SES.
# hashObject : a List with the hash objects (hash package) on which we have cached the generated statistics. 
#              SES requires this Object for the hash-based implementation of the statistics. This hashObject is produced or 
# updated by each run
#              of SES (if hash == TRUE) and it can be reused in next runs of SES.
# robust : Should the Fisher or the normal regression be robustly estimated? TRUE or FALSE 

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
# target(target variable), data(dataset), xIndex(x index), csIndex(cs index), dataInfo(list), 
# univariateModels(cached statistics for the univariate indepence test), hash(hash booleab), stat_hash(hash object), 
# pvalue_hash(hash object), robust=robust
# example: test(target, data, xIndex, csIndex, dataInfo=NULL, univariateModels=NULL, hash=FALSE, stat_hash=NULL, pvalue_hash=NULL, robust=robust = FALSE)
# output of each test: LIST of the generated pvalue, stat, flag and the updated hash objects.

# equal_case variable inside the code : it determines the method of the equivalent estimation
#   if equal_case = 1 then, if we have more than one equivalent vars in z , we select the one with the most closer pvalue to the pvalue of cvar
#   if equal_case = 2 then, if we have more than one equivalent vars in z , we select the one with the most minimum pvalue (>a)
#   else in any other case, if we have more than one equivalent vars in z , we select the first one
# In this version we support the equal_case = 3.
# 
# library(gRbase) #faster
# 
# #hashObject
# library(hash)

SES = function(target , dataset , max_k = 3 , threshold = 0.05 , test = NULL , user_test = NULL, hash=FALSE, hashObject=NULL, robust = FALSE)
{
  ##############################
  # initialization part of SES #
  ##############################
  faster = 0;
  #assign("gRbaseON",0,envir = .GlobalEnv)
  options(warn=-1)
  if(requireNamespace("gRbase", quietly = TRUE, warn.conflicts = FALSE) == TRUE)
  {
    #assign("gRbaseON",1,envir = .GlobalEnv)
    faster = 1;
  }
  options(warn=0);

  equal_case = 3;
  stat_hash = NULL;
  pvalue_hash = NULL;
  
  if(hash == TRUE)
  {
    if(requireNamespace("hash"))
    {
      if(is.null(hashObject))
      {
        stat_hash = hash();
        pvalue_hash = hash();
      }else if(class(hashObject) == "list"){
        stat_hash = hashObject$stat_hash;
        pvalue_hash = hashObject$pvalue_hash;
      }else{
        stop('hashObject must be a list of two hash objects (stat_hash, pvalue_hash)')
      }
    }else{
      cat('The hash version of SES requires the hash package');
      return(NULL);
    }
  }
  
  dataInfo = NULL;
  
  ###################################
  # dataset checking and initialize #
  ###################################
  
  if(!is.null(dataset))
  {
    if(class(dataset) == "matrix")
    {
      if(class(target) == "Surv")
      {
        stop('Invalid dataset class. For survival analysis provide a dataframe-class dataset');
      }
    }
    
    #check if dataset is an ExpressionSet object of Biobase package
    if(class(dataset) == "ExpressionSet")
    {
      #get the elements (numeric matrix) of the current ExpressionSet object.
      dataset = Biobase::exprs(dataset);
      dataset = t(dataset);#take the features as columns and the samples as rows
#     }else if(is.data.frame(dataset)){
#       if(class(target) != "Surv")
#       {
#         dataset = as.matrix(dataset);
#       }
    }else if((class(dataset) != "matrix") && (is.data.frame(dataset) == FALSE) ){
      stop('Invalid dataset class. It must be either a matrix, a dataframe or an ExpressionSet');
    }
  }
    if(is.null(dataset) || is.null(target) || (dim(as.matrix(target))[2] != 1 && class(target) != "Surv" ))
    {
      stop('invalid dataset or target (class feature) arguments.');
    }else{
      target = target;
    }
  
  #check for NA values in the dataset and replace them with the variable mean
  if(any(is.na(dataset)) == TRUE)
  {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values and they were replaced automatically by the variable (column) mean.")
    dataset = apply(dataset, 2, function(x){ x[which(is.na(x))] = mean(x,na.rm = TRUE) ; return(x)});
  }
  
  ##################################
  # target checking and initialize #
  ##################################
  
  targetID = -1;
  
  #check if the target is a string
  if (is.character(target) && length(target) == 1){
    findingTarget <- target == colnames(dataset);#findingTarget <- target %in% colnames(dataset);
    if(!sum(findingTarget)==1){
      warning('Target name not in colnames or it appears multiple times');
      return(NULL);
    }
    targetID <- which(findingTarget);
    target <- dataset[ , targetID];
  }
  
  #checking if target is a single number
  if (is.numeric(target) && length(target) == 1){
    if(target > dim(dataset)[2]){
      warning('Target index larger than the number of variables');
      return(NULL);
    }
    targetID <- target;
    target <- dataset[ , targetID];
  }
  
  ################################
  # test checking and initialize #
  ################################
  
  if(typeof(user_test) == "closure")
  {
    test = user_test;
  }else{
    #auto detect independence test in case of not defined by the user and the test is null or auto
    if(is.null(test) || test == "auto")
    {
      #if target is a factor then use the Logistic test
      if("factor" %in% class(target))
      {
        test = "testIndLogistic";
        if(is.ordered(target) == TRUE)
        {
          dataInfo$target_type = "ordinal";
          cat('\nTarget variable type: Ordinal')
        }else{
          if(length(unique(target)) == 2)
          {
            dataInfo$target_type = "binary"
            cat('\nTarget variable type: Binomial')
          }else{
            dataInfo$target_type = "nominal"
            cat('\nTarget variable type: Nominal')
          }
        }
      }else if(class(target) == "numeric" || class(target) == "matrix"){
        if(class(target) == "matrix")
        {
          if(dim(target)[2]!=1)
          {
            stop('Target can not be a matrix')
          }
        }
        
        if(identical(floor(target),target) == TRUE)
        {
          test = "testIndPois";
        }else{
          test = "testIndFisher";  
        }
      }else if(survival::is.Surv(target) == TRUE){
        test = "censIndLR";
      }else{
        stop('Target must be a factor, vector, matrix with one column or a Surv object');
      }
    }
    
    if(test == "testIndLogistic")
    {
      if(is.ordered(target) == TRUE)
      {
        dataInfo$target_type = "ordinal";
        cat('\nTarget variable type: Ordinal')
        
        if(requireNamespace("ordinal", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The testIndLogistic test requires the ordinal package for the ordered logistic regression method. Please install it.");
          return(NULL);
        }
        
      }else{
        if(length(unique(target)) == 2)
        {
          dataInfo$target_type = "binary"
          cat('\nTarget variable type: Binomial')
        }else{
          dataInfo$target_type = "nominal"
          cat('\nTarget variable type: Nominal')
          
          if(requireNamespace("nnet", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
          {
            cat("The testIndLogistic test requires the nnet package for the multinomial logistic regression method. Please install it.");
            return(NULL);
          }
        }
      }
    }
    
    cat("\nConditional independence test used: ");cat(test);cat("\n");
    
    #available conditional independence tests
    av_tests = c("testIndFisher", "testIndReg", "testIndRQ", "testIndBeta", "censIndLR", "testIndLogistic", "testIndPois", "testIndNB", "gSquare", "auto" ,  NULL);
    
    if(length(test) == 1) #avoid vectors, matrices etc
    {
      test = match.arg(test , av_tests ,TRUE);
      #convert to closure type
      if(test == "testIndFisher")
      {
        #an einai posostiaio target
        if ( all(target>0 & target<1) ){
          target = log(target/(1-target)) ## logistic normal 
        }
        
        test = testIndFisher;
      }
      else if (test == "testIndReg") ## It uses the F test
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        
        #an einai posostiaio target
        if ( all(target>0 & target<1) ){
          target = log(target/(1-target)) ## logistic normal 
        }
        
        test = testIndReg;
      }
      else if(test == "testIndBeta") ## beta regression for proportions
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        test = testIndBeta;
        if(requireNamespace("betareg", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The testIndBeta requires the betareg package. Please install it.");
          return(NULL);
        }
      }
      else if(test == "testIndRQ") ## beta regression for proportions
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        
        #an einai posostiaio target
        if ( all(target>0 & target<1) ){
          target = log(target/(1-target)) ## logistic normal 
        }
        
        test = testIndRQ;
        if(requireNamespace("quantreg", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The testIndRQ requires the quantreg package. Please install it.");
          return(NULL);
        }
      }
      else if (test == "testIndPois") ## Poisson regression
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        test = testIndPois;
      }
      else if (test == "testIndNB") ## Negative binomial regression
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        
        test = testIndNB;
      }
      else if(test == "censIndLR")
      {
        test = censIndLR;
        if(requireNamespace("survival", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The censIndLR requires the survival package. Please install it.");
          return(NULL);
        }
      }

      else if(test == "testIndLogistic")
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        test = testIndLogistic;
      }
      else if(test == "gSquare")
      {
        test = gSquare;
        if(requireNamespace("pcalg", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The gSquare test requires the pcalg package. Please install it.");
          return(NULL);
        }
      }
      #more tests here
    }else{
      stop('invalid test option');
    }
  }
  
  ###################################
  # options checking and initialize #
  ###################################
  
  #extracting the parameters
  max_k = floor(max_k);
  varsize = dim(dataset)[[2]];
  
  #option checking
  if((typeof(max_k)!="double") || max_k < 1)
  {
    stop('invalid max_k option');
  }
  if(max_k > varsize)
  {
    max_k = varsize;
  }
  if((typeof(threshold)!="double") || threshold <= 0 || threshold > 1)
  {
    stop('invalid threshold option');
  }
  if(typeof(equal_case)!="double")
  {
    stop('invalid equal_case option');
  }
  
  #######################################################################################
  
  #call the main SES function after the checks and the initializations
  results = InternalSES(target, dataset, max_k, threshold , test, equal_case, user_test, dataInfo, hash, varsize, stat_hash, pvalue_hash, targetID, faster, robust=robust);
  
  SESoutput <-new("SESoutput", selectedVars = results$selectedVars, selectedVarsOrder=results$selectedVarsOrder, queues=results$queues, signatures=results$signatures, hashObject=results$hashObject, pvalues=results$pvalues, stats=results$stats, max_k=results$max_k, threshold = results$threshold, runtime=results$runtime);
  
  return(SESoutput);
  
}

#########################################################################################################

InternalSES = function(target , dataset , max_k = 3 , threshold = 0.05 , test = NULL , equal_case = 3 , user_test = NULL , dataInfo = NULL , hash=FALSE, varsize, stat_hash, pvalue_hash, targetID, faster, robust=robust)
{
  #get the current time
  runtime = proc.time();
  
  #######################################################################################
  
  #univariate feature selection test
  
  if(is.loaded("fisher_uv") == TRUE && identical(test, testIndFisher) == TRUE && robust == FALSE)
  {
    a = .Fortran("fisher_uv", R = as.integer(dim(dataset)[1]), C = as.integer(dim(dataset)[2]), y = target, dataset = dataset,cs_cols = as.integer(0), pvalues = as.double(rep(0,dim(dataset)[2])), stats = as.double(rep(0,dim(dataset)[2])), targetID = as.integer(targetID))
    univariateModels = NULL;
    univariateModels$pvalue = a$pvalues;
    univariateModels$stat = a$stats;
    univariateModels$flag = rep(1,dim(dataset)[2]);
  	univariateModels$stat_hash = stat_hash;
  	univariateModels$pvalue_hash = pvalue_hash;
  }else{  
    univariateModels = univariateScore(target , dataset , test, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, targetID, robust=robust);
  }
  
  #univariateModels = univariateScore(target , dataset , test, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, targetID, robust=robust);
  pvalues = univariateModels$pvalue;
  stats = univariateModels$stat;
  flags = univariateModels$flag;
  stat_hash = univariateModels$stat_hash;
  pvalue_hash = univariateModels$pvalue_hash;
  #if we dont have any associations , return
  if(min(pvalues , na.rm=TRUE) > threshold) #or min(pvalues, na.rm=TRUE)
  {
    cat('No associations!');
    
    results = NULL;
    results$selectedVars = c();
    class(results$selectedVars) = "numeric";
    results$selectedVarsOrder = c();
    class(results$selectedVarsOrder) = "numeric";
    results$queues = c();
    class(results$queues) = 'list';
    results$signatures = matrix(nrow=1,ncol=1);
    class(results$signatures) = 'matrix';
    results$hashObject = NULL;
    class(results$hashObject) = 'list';
    
    results$pvalues = pvalues;
    results$stats = stats;
    results$max_k = max_k;
    results$threshold = threshold;
    runtime = proc.time() - runtime;
    results$runtime = runtime;
    
    return(results);
  }
  
  
  #Initialize the data structs
  selectedVars = rep(0,varsize);
  selectedVarsOrder = rep(0,varsize);
  queues = vector('list' , varsize);
  
  queues <- lapply(1:varsize , function(i){queues[[i]] = i;})
  
  #select the variable with the highest association
  selectedVar = which(flags == 1 & stats == stats[[which.max(stats)]]);
  selectedVars[selectedVar] = 1;
  selectedVarsOrder[selectedVar] = 1; #CHANGE
  
  #lets check the first selected var
  #cat('First selected var: %d, p-value: %.6f\n', selectedVar, pvalues[selectedVar]);
  
  #remaining variables to be considered
  remainingVars = rep(1,varsize);
  remainingVars[selectedVar] = 0;
  remainingVars[pvalues > threshold] = 0;
  if (targetID > 0){
    remainingVars[targetID] = 0;
  }
  
  ################ main ses loop ################
  
  #main SES loop
  #loop until there are not remaining vars
  loop = any(as.logical(remainingVars));
  
  while(loop)
  {
    #lets find the equivalences
    IdEq_results <- IdentifyEquivalence(equal_case , queues , target , dataset , test , threshold , max_k , selectedVars , pvalues , stats , remainingVars , univariateModels, selectedVarsOrder, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, faster, robust=robust);
    queues = IdEq_results$queues;
    selectedVars = IdEq_results$selectedVars;
    remainingVars = IdEq_results$remainingVars;
    pvalues = IdEq_results$pvalues;
    stats = IdEq_results$stats;
    stat_hash=IdEq_results$stat_hash;
    pvalue_hash=IdEq_results$pvalue_hash;
    
    #lets find the variable with the max min association
    max_min_results = max_min_assoc(target, dataset , test , threshold , max_k , selectedVars , pvalues , stats , remainingVars , univariateModels, selectedVarsOrder, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, faster, robust=robust);
    selectedVar = max_min_results$selected_var;
    selectedPvalue = max_min_results$selected_pvalue;
    remainingVars = max_min_results$remainingVars;
    pvalues = max_min_results$pvalues;
    stats = max_min_results$stats;
    stat_hash=max_min_results$stat_hash;
    pvalue_hash=max_min_results$pvalue_hash;
    
    #if the selected variable is associated with target , add it to the selected variables
    if(selectedPvalue <= threshold)
    {
      selectedVars[selectedVar] = 1;
      selectedVarsOrder[selectedVar] = max(selectedVarsOrder) + 1;
      remainingVars[selectedVar] = 0;
    }
    
    loop = any(as.logical(remainingVars));
  }
  
  #lets find the variables to be discarded
  IdEq_results <- IdentifyEquivalence(equal_case , queues , target , dataset , test , threshold , max_k , selectedVars , pvalues , stats , remainingVars , univariateModels, selectedVarsOrder, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, faster, robust=robust);
  queues = IdEq_results$queues;
  selectedVars = IdEq_results$selectedVars;
  pvalues = IdEq_results$pvalues;
  stats = IdEq_results$stats;
  remainingVars = IdEq_results$remainingVars;
  stat_hash=IdEq_results$stat_hash;
  pvalue_hash=IdEq_results$pvalue_hash;
  
  selectedVarsOrder[which(!selectedVars)] = varsize;#
  numberofSelectedVars = sum(selectedVars);#
  selectedVarsOrder = sort(selectedVarsOrder);#
  #   selectedVars = selectedVarsOrder[1:numberofSelectedVars];
  
  #queues correctness
  all_queues = queues
  queues = queues[which(selectedVars==1)];
  
  queues <- lapply(1:length(queues) , function(i){queues[[i]] = unique(queues[[i]]);});
  
  #adjusting the results
  if(targetID > 0)
  {
    toAdjust <- which(selectedVars > targetID);
    selectedVars[toAdjust] = selectedVars[toAdjust] + 1;
  }
  
  
  results = NULL;
  results$selectedVars = which(selectedVars == 1);
  
  svorder = sort(pvalues[results$selectedVars] , index.return = TRUE);
  svorder = results$selectedVars[svorder$ix];
  results$selectedVarsOrder = svorder;
  
  results$queues = queues;
  results$signatures = as.matrix(do.call(expand.grid, results$queues))
  hashObject = NULL;
  hashObject$stat_hash = stat_hash;
  hashObject$pvalue_hash = pvalue_hash;
  results$hashObject = hashObject;
  class(results$hashObject) = 'list';
  
  results$pvalues = pvalues;
  results$stats = stats;
#   results$all_queues = all_queues;
#   already known
#   results$data = dataset;
#   results$target = target;
#   results$test = test;
  results$max_k = max_k;
  results$threshold = threshold;
  
  runtime = proc.time() - runtime;
  results$runtime = runtime;
  
  
  
  return(results);
}

#univariate feature selection ( uncoditional independence )

univariateScore = function(target , dataset , test, hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, targetID, robust=robust)
{
  #how many tests
  nTests = ncol(dataset);
  
  #data structure to be returned
  univariateModels = NULL;
  univariateModels$pvalue = rep(1,nTests)
  univariateModels$stat = rep(0,nTests)
  univariateModels$flag = rep(0,nTests);
  #univariateModels$uniModelFit = rep(NA,nTests);
  
  test_results = NULL;
  #for way to initialize the univariateModel
  #FOR LOOP IS FASTER THAN VAPPLY IN THIS CASE (try apply withm margin 2)
  for(i in 1:nTests)
  {
    #arguments order for any CI test are fixed
    if (i != targetID){
      test_results = test(target , dataset , i, 0 , dataInfo=dataInfo, hash=hash, stat_hash=stat_hash, pvalue_hash=pvalue_hash, robust=robust)
      univariateModels$pvalue[[i]] = test_results$pvalue;
      univariateModels$stat[[i]] = test_results$stat;
      univariateModels$flag[[i]] = test_results$flag;
      univariateModels$stat_hash = test_results$stat_hash
      univariateModels$pvalue_hash = test_results$pvalue_hash      
    }else{
      univariateModels$pvalue[[i]] = 1;
      univariateModels$stat[[i]] = 0;
      univariateModels$flag[[i]] = 1;
    }
  }
  
  return(univariateModels);
}

#########################################################################################################

#just like matlab's nchoosek but with transposed result
#Its a slightly different from combn() 
#(ex. (nchoosekm(4,2) != nchoosekm(1:4,2) like nchoosek in matlab , combn(4,2) == combn(1:4,2)))
nchoosekm = function(cs , k, faster) #i can also pass the compFun arg for selecting
{ 
  if(length(cs) == 1) #if not vector
  {
    res = choose(cs , k); #or nchoosek
  }else{
    if(faster == 1)
    {
      res = gRbase::combnPrim(cs,k); #combs(as.vector(cs),k); #combnPrim
    }else
    {
      res = combn(cs,k);
    }

  }
  return(res);
}

IdentifyEquivalence = function(equal_case , queues , target , dataset , test , threshold , max_k , selectedVars , pvalues , stats , remainingVars , univariateModels, selectedVarsOrder, hash, dataInfo, stat_hash, pvalue_hash, faster, robust=robust)
{ 
  varsToBeConsidered = which(selectedVars==1 | remainingVars==1); #CHANGE
  lastvar = which(selectedVarsOrder == max(selectedVarsOrder))[1]; #CHANGE
  
  #for every variable to be considered
  for(cvar in varsToBeConsidered)
  {
    #if var is the last one added, no sense to perform the check
    if(cvar == lastvar) #CHANGE
    {
      next;
    }
    
    #the maximum conditioning set
    selectedVarsIDs = which(selectedVars == 1);
    cs = setdiff(selectedVarsIDs , cvar);
    k = min(c(max_k , length(cs)));
    
    #for all possible subsets at most k size
    #problem in 1:k if K=0 - solve with while temporary
    klimit = 1;
    while(klimit <= k)
    {
      #set to investigate
      tempCS = setdiff(cs, lastvar)#CHANGE
      if(klimit == 1) #CHANGE
      {
        subsetcsk = as.matrix(lastvar); #CHANGE
      }else{
        subsetcsk = as.matrix(nchoosekm(tempCS,klimit-1,faster)); #CHANGE
        numSubsets = dim(subsetcsk)[2]; #CHANGE
        subsetcsk = rbind(subsetcsk, lastvar*rep(1,numSubsets));#CHANGE
      }
      
      #flag to get out from the outermost loop
      breakFlag = FALSE;
      
      #or combs or nchoosekm
      #subsetcsk = as.matrix(nchoosekm(cs,klimit));
      
      for(i in 1:ncol(subsetcsk))
      {
        z = subsetcsk[,i];
        z = t(t(z));
        
        cur_results = test(target , dataset , cvar, z , dataInfo=dataInfo, univariateModels, hash = hash, stat_hash, pvalue_hash, robust=robust);
        stat_hash = cur_results$stat_hash;
        pvalue_hash = cur_results$pvalue_hash;
        
        #check if the pvalues and stats should be updated
        if(compare_p_values(pvalues[[cvar]] , cur_results$pvalue , stats[[cvar]] , cur_results$stat))
        {
          pvalues[[cvar]] = cur_results$pvalue;
          stats[[cvar]] = cur_results$stat;
        }
        
        #if there is any subset that make var independent from y,
        #then let's throw away var; moreover, we also look for
        #equivalent variables. Note that we stop after the first 
        #z such that pvalue_{var, y | z} > threshold
        if(cur_results$flag && cur_results$pvalue > threshold)
        {
          remainingVars[[cvar]] = 0;
          selectedVars[[cvar]] = 0;
          queues = identifyTheEquivalent(equal_case , queues , target , dataset , cvar , z , test , threshold , univariateModels , pvalues, hash, dataInfo, stat_hash, pvalue_hash, robust=robust);
          breakFlag = TRUE;
          break;
        }
      }
      if(breakFlag == TRUE)
      {
        break;
      }else
      {
        klimit = klimit + 1;
      }
    }
  }
  results <- list(pvalues = pvalues , stats = stats , queues = queues , selectedVars = selectedVars , remainingVars = remainingVars, stat_hash = stat_hash, pvalue_hash = pvalue_hash);
  return(results);
}

#########################################################################################################

identifyTheEquivalent = function(equal_case , queues , target , dataset , cvar , z , test , threshold , univariateModels , pvalues, hash, dataInfo, stat_hash, pvalue_hash, robust=robust)
{
  z = t(z);
  
  #case 3
  #if we have more than one equivalent vars in z , we select the first one
  #for loop if we dont use the below lapply function
  #equalsY = 0;
  for(i in 1:ncol(z))
  {
    w = z[,i];
    w = t(t(w));
    zPrime = c(setdiff(z , w) , cvar);
    
    cur_results = test(target , dataset , w, zPrime , dataInfo=dataInfo, univariateModels, hash = hash, stat_hash, pvalue_hash, robust=robust);
    
    if(cur_results$flag & (cur_results$pvalue > threshold))
    {  
      queues[[w]] = as.matrix(c(queues[[w]] , queues[[cvar]]));
      break;
      #equalsY = equalsY+1;
    }
  }
  #cat("\nEquals Ys in the current z subset = %d",equalsY);
  return(queues);
}

apply_ideq = function(i , queues , target , dataset , cvar , z , test , threshold , univariateModels, hash, dataInfo, stat_hash, pvalue_hash, robust=robust)
{
  w = z[,i];
  w = t(t(w));
  zPrime = c(setdiff(z , w) , cvar);
  
  cur_results = test(target , dataset , w, zPrime , dataInfo=dataInfo, univariateModels, hash = hash, stat_hash, pvalue_hash, robust=robust);
  
  if(cur_results$flag & (cur_results$pvalue > threshold))
  {
    queues[[w]] = as.matrix(c(queues[[w]] , queues[[cvar]]));
    return(queues[[w]]);
  }else{
    return(NA);
  }
}

#########################################################################################################

compare_p_values = function(pval, pval2, stat, stat2)
{
  if(length(pval) == 0 | length(pval2) == 0 | length(stat) == 0 | length(stat2) ==0)
  {
    return(FALSE);
  }else{
    if(is.na(pval2)==TRUE | is.na(stat2)==TRUE | is.na(pval)==TRUE | is.na(stat)==TRUE)
    {
      pval2 = 0.0;
      return(FALSE);#(pval < pval2);
    }else{
      if (pval <= 2e-16 | pval2 <= 2e-16){
        return(stat > stat2);
      }else{
        return(pval < pval2);
      }
    }
  }
}

#########################################################################################################

max_min_assoc = function(target, dataset , test , threshold , max_k , selectedVars , pvalues , stats , remainingVars , univariateModels, selectedVarsOrder, hash, dataInfo, stat_hash, pvalue_hash, faster, robust=robust)
{
  #Initialize
  selected_var = -1;
  selected_pvalue = 2;
  selected_stat = 0;
  
  varsToIterate = which(remainingVars==1);
  for(cvar in varsToIterate)
  {
    mma_res = min_assoc(target, dataset , test , max_k , cvar , selectedVars , pvalues , stats , univariateModels , selectedVarsOrder, hash, dataInfo, stat_hash, pvalue_hash, faster, robust=robust);
    pvalues = mma_res$pvalues;
    stats = mma_res$stats;
    stat_hash = mma_res$stat_hash;
    pvalue_hash = mma_res$pvalue_hash;
    
    
    if(mma_res$pvalue > threshold)
    {
      remainingVars[[cvar]] = 0;
    }
    
    if(compare_p_values(mma_res$pvalue , selected_pvalue , mma_res$stat , selected_stat))
    {
      selected_var = cvar;
      selected_pvalue = mma_res$pvalue;
      selected_stat = mma_res$stat;
    }
  }
  results <- list(selected_var = selected_var , selected_pvalue = selected_pvalue , remainingVars = remainingVars , pvalues = pvalues , stats = stats, stat_hash=stat_hash, pvalue_hash = pvalue_hash);
  return(results); 
}

#########################################################################################################

min_assoc = function(target , dataset , test ,  max_k , cvar , selectedVars , pvalues , stats , univariateModels , selectedVarsOrder, hash, dataInfo, stat_hash, pvalue_hash, faster, robust=robust)
{
  #initialization
  #baseline values
  #   ma_pvalue = univariateModels$pvalue[[cvar]];
  #   ma_stat = univariateModels$stat[[cvar]];
  ma_pvalue = pvalues[[cvar]]; #CHANGE
  ma_stat = stats[[cvar]]; #CHANGE
  
  selectedVars = which(selectedVars==1);
  #max size of the condiotioning test
  k = min(c(max_k , length(selectedVars)));
  
  ck = 1;
  while(ck<=k)
  {
    #lastvar = unique(which(selectedVarsOrder == max(selectedVarsOrder)));
    lastvar = which(selectedVarsOrder == max(selectedVarsOrder))[1]; #CHANGE
    
    tempCS = setdiff(selectedVars, lastvar) #CHANGE
    if(ck == 1) #CHANGE
    {
      subsetcsk = as.matrix(lastvar); #CHANGE
    }else{
      subsetcsk = as.matrix(nchoosekm(tempCS,ck-1,faster)); #CHANGE
      numSubsets = dim(subsetcsk)[2]; #CHANGE
      subsetcsk = rbind(subsetcsk, lastvar*rep(1,numSubsets)); #CHANGE
    }
    
    #or combs or nchoosekm
    #subsetcsk = as.matrix(nchoosekm(1:length(selectedVars),ck));
    
    #subsetcsk = t(subsetcsk);
    for(i in 1:ncol(subsetcsk))
    {
      s = subsetcsk[,i];
      s = t(t(s));
      
      cur_results = test(target , dataset , cvar, selectedVars[s] , dataInfo=dataInfo, univariateModels, hash = hash, stat_hash, pvalue_hash, robust=robust);
      stat_hash = cur_results$stat_hash;
      pvalue_hash = cur_results$pvalue_hash;
      
      #check if the pvalues and stats should be updated
      if(cur_results$flag == 1 & !compare_p_values(cur_results$pvalue, ma_pvalue, cur_results$stat , ma_stat))
      {
        ma_pvalue = cur_results$pvalue;
        pvalues[[cvar]] = cur_results$pvalue;
        
        ma_stat = cur_results$stat;
        stats[[cvar]] = cur_results$stat;
      }
    }
    ck = ck+1;
  }
  results <- list(pvalue = ma_pvalue , stat = ma_stat , pvalues = pvalues , stats = stats, stat_hash=stat_hash, pvalue_hash  = pvalue_hash);
  return(results);
}

# .onAttach <- function(libname, pkgname){
#   # do whatever needs to be done when the package is loaded
#   packageStartupMessage( "Loading MXM package version 0.2, thanks for downloading." )
#   #load the dll files from the fortran code for the package
#   #library.dynam("MXM", pkgname, libname)
# }
