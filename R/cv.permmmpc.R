#Cross Validation for mmpc
#INPUT
#target as in mmpc
#dataset as in mmpc
#kfolds: number of folds (integer)
#folds: already defined folds of the data to use (a list). If NULL the folds created internally with the same function
#alphas: vector of mmpc alphas hyper parameters used in CV. Default is c(0.1, 0.05, 0.01)
#maxk_s: vector of mmpc max_ks parameters used in CV. Default is c(3, 2)
#task: character, it can be "C" for classification (logistic regression classifier), "R" for regression (linear regression classifier), "S" for cox survival analysis (cox regression classifier)
#metric: a metric function provided by the user or auto defined due to the task. It may be NULL or a function in the form of other metric functions (e.g., mse.mxm). For example the default for the classification task is auc.mxm but the user can also define acc.mxm (based on the accuracy metric) that is supported on the package. Or the user can make his own metric function that follows the signature and the inputs, outputs of ours.
#modeler: a modeling function provided by the user or auto defined due to the task if it is NULL (e.g., lm.mxm)
#mmpc_test: A function object that defines the test used in the mmpc function (see mmpc help page for more). If it is NULL, its is auto defined due to the task.
#OUTPUT
#a list called best_model with the below slots
#cv_results_all: a list with the predictions, performances and the signatures for each fold of each configuration (i.e cv_results_all[[3]]$performances[1] indicates the performance of the 1st fold with the 3d configuration of mmpc)
#best_performance: numeric, the best average performance
#best_configuration: the best configuration of mmpc (a list with the slots id, a, max_k)
cv.permmmpc <- function(target, dataset, wei = NULL, kfolds = 10, folds = NULL, alphas = c(0.1, 0.05, 0.01), max_ks = c(3, 2), task = NULL, 
                        metric = NULL, metricbbc = NULL, modeler = NULL, mmpc_test = NULL, R = 999, ncores = 1, B = 1) 
{
  if ( ncores > 1 ) {  ## multi-threaded task
    result = cvpermmmpc.par(target, dataset, wei = wei, kfolds = kfolds, folds = folds, alphas = alphas, max_ks = max_ks, task = task, metric = metric, modeler = modeler, mmpc_test = mmpc_test, R = R, ncores = ncores)
    
  } else { ## one core task
    
    if (is.null(alphas) )  alphas <- c(0.1, 0.05, 0.01)
    if (is.null(max_ks) )  max_ks <- c(3, 2)  
    
    alphas = sort(alphas, decreasing = TRUE)
    max_ks = sort(max_ks, decreasing = TRUE)
    nAlpha <- length(alphas);
    nMax_ks <- length(max_ks);
    #defining the mmpc configurations
    nmmpcConfs <- nAlpha * nMax_ks;
    mmpc_configurations <- vector("list" , nmmpcConfs);
    i <- 0;
    for (a in alphas) {
      for (k in max_ks) {
        configuration <- NULL;
        i <- i + 1;
        configuration$id <- i;
        configuration$a <- a;
        configuration$max_k <- k;
        mmpc_configurations[[i]] <- configuration;
      }
    } 
    
    if ( is.null(folds) ) {
      if (task == "R" ) {
        folds <- generatefolds(target, nfolds = kfolds, stratified = FALSE, seed = FALSE)
      } else if (task == "S") {
        folds <- generatefolds(target[, 1], nfolds = kfolds, stratified = FALSE, seed = FALSE)
      } else   folds <- generatefolds(target, nfolds = kfolds, stratified = TRUE, seed = FALSE)
    } else  kfolds <- length( folds );
    
    if ( is.null(task) ) {
      stop("Please provide a valid task argument 'C'-classification, 'R'-regression, 'S'-survival.")
      #to do: select automatically the appropriate task due to the data, target
    } else if (task == 'C') {
      
      #Classification task (logistic regression)
      if ( is.null(metric) ) {
        metricFunction <- auc.mxm;
      } else   metricFunction <- metric;
      
      if ( is.null(modeler) ) {
        modelerFunction <- glm.mxm;
      } else   modelerFunction <- modeler;
      
      if ( is.null(mmpc_test) ) {
        test <- 'permLogistic';
      } else  test <- mmpc_test;
      
    } else if (task == 'R') {
      
      #Regression task 
      if ( is.null(metric) ) {
        metricFunction <- mse.mxm;
      } else  metricFunction <- metric;
      
      if ( is.null(modeler) ) {
        modelerFunction <- lm.mxm;
      } else  modelerFunction <- modeler;
      
      if (is.null(mmpc_test) ) {
        test = 'permFisher';
      } else  test <- mmpc_test;
      
    } else if(task == 'S') {
      
      #cox survival analysis (cox regression)
      if ( is.null(metric) ) {
        metricFunction <- ci.mxm;
      } else  metricFunction <- metric;
      
      if ( is.null(modeler) ) {
        modelerFunction <- coxph.mxm;
      } else  modelerFunction <- modeler;
      
      if ( is.null(mmpc_test) ) {
        test = "permCR";
      } else  test <- mmpc_test;
      
    } else  stop("Please provide a valid task argument 'C'-classification, 'R'-regression, 'S'-survival.")
    
    nmmpcConfs = length(mmpc_configurations)
    #merging mmpc configuration lists and create the general cv results list
    conf_mmpc <- vector("list" , nmmpcConfs)
    for (i in 1:nmmpcConfs) {
      conf_mmpc[[i]]$configuration <- mmpc_configurations[[i]]
      conf_mmpc[[i]]$preds <- vector('list', kfolds)
      conf_mmpc[[i]]$performances <- vector('numeric', kfolds)
      conf_mmpc[[i]]$variables <- vector('list', kfolds)
    }
    
    tic <- proc.time()
    
    for (k in 1:kfolds) {
      #print(paste('CV: Fold', k, 'of', kfolds));
      train_samples <- c();
      for(i in which(c(1:kfolds) != k))  train_samples = c( train_samples, folds[[ i ]] ) 
      #leave one fold out each time as a test set and the rest as train set
      train_set <- dataset[train_samples, ] #Set the training set
      train_target <- target[train_samples]
      wtrain <- wei[train_samples]
      test_set <- dataset[ folds[[k]], ] #Set the validation set
      test_target <- target[ folds[[k]] ]
      #mmpc hashmap
      mmpcHashMap = NULL;
      mmpcini = NULL
      
      #for each conf of mmpc
      for (mmpc_conf_id in 1:nmmpcConfs) {
        
        #mmpc options
        threshold <- mmpc_configurations[[mmpc_conf_id]]$a;
        max_k <- mmpc_configurations[[mmpc_conf_id]]$max_k;
        #running mmpc
        results <- perm.mmpc(train_target, train_set, max_k, threshold, test = test, ini = mmpcini, wei = wtrain, hash = TRUE, hashObject = mmpcHashMap, R = R)
        mmpcini <- results@univ
        mmpcHashMap <- results@hashObject;
        variables <- results@selectedVars;
        conf_mmpc[[mmpc_conf_id]]$variables[[k]] <- variables
        #get the data of the reference signature (i.e the selected variables)
        curr_sign <- as.vector(variables)
        #curr_sign <- as.matrix(results@selectedVars) #in case that the signature slot is not returned due to lack of memory. See Internalmmpc final part.
        sign_data <- train_set[, curr_sign, drop = FALSE]
        sign_test <- test_set[, curr_sign, drop = FALSE]
        
        if( length(variables) > 0 ) {
          #generate a model due to the task and find the performance
          #logistic model for a classification task, linear model for the regression task and a cox model for the survival task
          moda <- modelerFunction(train_target, sign_data, sign_test, wei = wtrain)
          preds <- moda$preds
          theta <- moda$theta
        } else  {
          moda <- modelerFunction(train_target, rep(1, nrow(sign_data)), rep(1, nrow(sign_test)), wei = wtrain)
          preds <- moda$preds
          theta <- moda$theta
        }
        if ( is.null(preds) ) {
          conf_mmpc[[mmpc_conf_id]]$preds[[k]] <- NULL
          conf_mmpc[[mmpc_conf_id]]$performances[k] <- NA
        } else {
          performance = metricFunction(preds, test_target, theta)
          conf_mmpc[[mmpc_conf_id]]$preds[[k]] <- preds
          conf_mmpc[[mmpc_conf_id]]$performances[k] <- performance
        }
      }
      #clear the hashmap and garbages
      if ( !is.null(mmpcHashMap$pvalue_hash) )   mmpcHashMap$pvalue_hash <- NULL
      if ( !is.null(mmpcHashMap$stat_hash) )     mmpcHashMap$stat_hash <- NULL
    }
    
    #finding the best performance for the metric  
    index = 1;
    best_perf = mean(conf_mmpc[[1]]$performances, na.rm = TRUE);
    for ( i in 2:length(conf_mmpc) ) {
      averagePerf <- mean(conf_mmpc[[i]]$performances, na.rm = TRUE);
      if ( !is.na(averagePerf)  &  !is.na(best_perf) ) {
        if (averagePerf < best_perf) {
          best_perf <- averagePerf;
          index <- i;
        }
      }
    }
    
    #recording the best results
    best_model <- NULL
    best_model$cv_results_all <- conf_mmpc;
    best_model$best_performance <- best_perf
    #TT
    mat <- matrix(nrow = length(best_model[[ 1 ]]), ncol = kfolds)
    for ( i in 1:nrow(mat) )  mat[i, ] <- as.vector( best_model[[ 1 ]][[ i ]]$performances )  
    
    opti <- rowMeans(mat, na.rm = TRUE)
    bestpar <- which.max(opti)
    best_model$best_configuration = conf_mmpc[[bestpar]]$configuration
    best_model$best_performance <- max( opti )
    best_model$bbc_best_performance <- NULL
    
    if ( B > 1) {
      if (task == "S")  {
        n <- 0.5 * length(target) 
      } else  n <- length(target)
      predictions <- matrix(0, nrow = n, ncol = nmmpcConfs)
      for (i in 1:nmmpcConfs)  predictions[, i] <- unlist( best_model$cv_results_all[[ i ]]$preds )
      best_model$bbc_best_performance <- MXM::bbc(predictions, target[unlist(folds)], metric = metricbbc, B = B )$bbc.perf
    }

    best_model$runtime <- proc.time() - tic 
    result <- best_model
  }
  
  result 
}














