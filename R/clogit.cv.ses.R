#Cross Validation for SES
#INPUT
#target as in SES
#dataset as in SES
#kfolds: number of folds (integer)
#folds: already defined folds of the data to use (a list). If NULL the folds created internally with the same function
#alphas: vector of SES alphas hyper parameters used in CV. Default is c(0.1, 0.05, 0.01)
#maxk_s: vector of SES max_ks parameters used in CV. Default is c(3, 2)
#task: character, it can be "C" for classification (logistic regression classifier), "R" for regression (linear regression classifier), "S" for cox survival analysis (cox regression classifier)
#metric: a metric function provided by the user or auto defined due to the task. It may be NULL or a function in the form of other metric functions (e.g., mse.mxm). For example the default for the classification task is auc.mxm but the user can also define acc.mxm (based on the accuracy metric) that is supported on the package. Or the user can make his own metric function that follows the signature and the inputs, outputs of ours.
#modeler: a modeling function provided by the user or auto defined due to the task if it is NULL (e.g., lm.mxm)
#ses_test: A function object that defines the test used in the SES function (see SES help page for more). If it is NULL, its is auto defined due to the task.
#OUTPUT
#a list called best_model with the below slots
#cv_results_all: a list with the predictions, performances and the signatures for each fold of each configuration (i.e cv_results_all[[3]]$performances[1] indicates the performance of the 1st fold with the 3d configuration of SES)
#best_performance: numeric, the best average performance
#best_configuration: the best configuration of SES (a list with the slots id, a, max_k)
clogit.cv.ses <- function(target, dataset, kfolds = 10, folds = NULL, alphas = c(0.1, 0.05, 0.01), max_ks = c(3, 2), 
                   metricbbc = NULL, B = 1, ncores = 1)
{
  
  if ( ncores > 1 ) {  ## multi-threaded task
   ##  result <- clogit.cvses.par(target, dataset, kfolds = kfolds, folds = folds, alphas = alphas, max_ks = max_ks, 
   ##                             metricbbc = metricbbc, B = 1, ncores = ncores)
    
  } else { ## one core task     

    if ( is.null(alphas) )  alphas <- c(0.1, 0.05, 0.01)
    if ( is.null(max_ks) )  max_ks <- c(3, 2)  
    
    alphas <- sort(alphas, decreasing = TRUE)
    max_ks <- sort(max_ks, decreasing = TRUE)
    nAlpha <- length(alphas);
    nMax_ks <- length(max_ks);
    #defining the SES configurations
    nSESConfs <- nAlpha*nMax_ks;
    SES_configurations <- vector("list" , nSESConfs);
    i <- 0;
    for (a in alphas) {
      for (k in max_ks) {
        configuration <- NULL;
        i <- i + 1;
        configuration$id <- i;
        configuration$a <- a;
        configuration$max_k <- k;
        SES_configurations[[i]] <- configuration;
      }
    }
    
    if ( is.null(folds) ) {
      index <- target[, 2]
      folds <- generatefolds( unique(index), nfolds = kfolds, stratified = FALSE, seed = FALSE)
      for (i in 1:kfolds) {
        j <- 1
        a <- NULL
        while ( j < length(folds[[ i ]]) ) {
          a <- c( a, which( index == folds[[i]][j] ) )
          j <- j + 1
        }
        folds[[ i ]] <- a
      }  ##  end  for (i in 1:kfolds) { 
    } else  kfolds <- length( folds )
    
    metricFunction <- mci.mxm
    modelerFunction <- clogit.mxm
    test <- "testIndClogit"

    nSESConfs <- length(SES_configurations)
    conf_ses <- vector("list", nSESConfs)
    
    for(i in 1:nSESConfs){
      conf_ses[[i]]$configuration <- SES_configurations[[i]]
      conf_ses[[i]]$preds <- vector('list', kfolds)
      conf_ses[[i]]$performances <- vector('numeric', kfolds)
      conf_ses[[i]]$signatures <- vector('list', kfolds)
    }
    ####################
    ## Start the CV procedure
    ####################
    tic <- proc.time()
    
    for (k in 1:kfolds) {
      #print(paste('CV: Fold', k, 'of', kfolds));
      train_samples <- c();
      for ( i in which(c(1:kfolds) != k) )   train_samples = c( train_samples, folds[[ i ]] )
      #leave one fold out each time as a test set and the rest as train set
      train_set <- dataset[train_samples, ] #Set the training set
      train_target <- target[train_samples, ]
      test_set <- dataset[folds[[k]], ] #Set the validation set
      test_target <- target[ folds[[k]], ]
      #SES hashmap
      SESHashMap = NULL;
      sesini = NULL
      #for each conf of SES
      for(ses_conf_id in 1:nSESConfs){
        #SES options
        threshold <- SES_configurations[[ses_conf_id]]$a;
        max_k <- SES_configurations[[ses_conf_id]]$max_k;
        #running SES
        results <- SES(train_target, train_set, max_k, threshold, test = test, ini = sesini, hash = TRUE, hashObject = SESHashMap)
        sesini <- results@univ
        SESHashMap <- results@hashObject;
        signatures <- results@signatures;
        #recording the selected signatures
        conf_ses[[ses_conf_id]]$signatures[[k]] <- signatures;
        #get the data of the reference signature (i.e the selected variables)
        curr_sign <- as.matrix(signatures[1, ])
        #curr_sign <- as.matrix(results@selectedVars) #in case that the signature slot is not returned due to lack of memory. See InternalSES final part.
        sign_data <- train_set[, curr_sign, drop = FALSE]
        sign_test <- test_set[, curr_sign, drop = FALSE]
        
        if ( dim(signatures)[1] >= 1 & length(results@selectedVars ) > 0 ) {
          #generate a model due to the task and find the performance
          #logistic model for a classification task, linear model for the regression task and a cox model for the survival task
          preds <- modelerFunction(train_target, sign_data, sign_test)$preds
        } else  {
          preds <- NULL
        }  
        
        if ( is.null(preds) ) {
          conf_ses[[ses_conf_id]]$preds[[k]] <- NULL
          conf_ses[[ses_conf_id]]$performances[k] <- NA
        } else {
          performance <- metricFunction(preds, test_target)
          conf_ses[[ses_conf_id]]$preds[[k]] <- preds
          conf_ses[[ses_conf_id]]$performances[k] <- performance
        }
      }
      #clear the hashmap and garbages
      if ( !is.null(SESHashMap$pvalue_hash) )   SESHashMap$pvalue_hash <- NULL
      if ( !is.null(SESHashMap$stat_hash) )     SESHashMap$stat_hash <- NULL
    }
    #finding the best performance for the metric  
    index = 1;
    best_perf = mean(conf_ses[[1]]$performances, na.rm = TRUE);
    
    for ( i in 2:length(conf_ses) ) {
      averagePerf <- mean( conf_ses[[i]]$performances, na.rm = TRUE );
      if ( !is.na(averagePerf)   &  !is.na(best_perf) ) {
        if ( averagePerf < best_perf ) {
          best_perf <- averagePerf;
          index <- i;
        }
      }
    }
    #recording the best results
    best_model <- NULL
    best_model$cv_results_all <- conf_ses;
    best_model$best_performance <- best_perf
    #TT
    mat <- matrix(nrow = length(best_model[[ 1 ]]), ncol = kfolds)
    for ( i in 1:dim(mat)[1] )  mat[i, ] <- as.vector( best_model[[ 1 ]][[ i ]]$performances )  
    
    opti <- Rfast::rowmeans(mat)
    bestpar <- which.max(opti)
    best_model$best_configuration <- conf_ses[[bestpar]]$configuration
    best_model$best_performance <- max( opti )
    best_model$bbc_best_performance <- NULL
    
    #if ( B > 1) {
    #  if (task == "S")  {
    #    n <- 0.5 * length(target) 
    #  } else  n <- length(target)
    #  predictions <- matrix(0, nrow = n, ncol = nSESConfs)
    #  for (i in 1:nSESConfs)  predictions[, i] <- unlist( best_model$cv_results_all[[ i ]]$preds )
    #  best_model$bbc_best_performance <- MXM::bbc(predictions, target[unlist(folds)], metric = metricbbc, B = B )$bbc.perf
    #}
    
    best_model$runtime <- proc.time() - tic 
    
    result <- best_model
  }
  
  result
}














## matched C-Index for matched case controls (conditional logistic regression)
mci.mxm <- function(predictions, test_target) {
  case <- test_target[, 1]  ## case control, 0 is the control  
  id <- test_target[, 2] #the patient id
  ti <- numeric( length( unique(id) ) )
  for ( i in 1:length(ti) ) {
    est <- predictions[ id == id[i] ]
    index <- case[id == id[i]]  
    ti[i] <- mean(est[index == 1] >est[index == 0] ) + 0.5 * mean(est[index == 1] == est[index == 0] )
  }
  mean(ti)
}


## conditional logistic regression
clogit.mxm <- function(train_target, sign_data, sign_test) {
  x <- model.matrix(~., data = data.frame(sign_data) )[, -1]
  id <- train_target[, 2] #the patient id
  case <- train_target[, 1]  ## case control, 0 is the control
  sign_model <- survival::clogit(case ~ . + strata(id), data = data.frame(x) )
  x <- model.matrix(~., data = data.frame(sign_test) )[, -1]
  preds <- x %*% sign_model$coefficients
  list(preds = preds)
}

