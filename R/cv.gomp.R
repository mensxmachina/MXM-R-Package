cv.gomp <- function(target, dataset, kfolds = 10, folds = NULL, tol = seq(4, 9, by = 1), task = "C", metric = NULL,
                    metricbbc = NULL, modeler = NULL, test = NULL, method = "ar2", B = 1) {

    if ( is.null(tol) )   tol <- seq(4, 9, by = 1)
    tol <- sort(tol)
    ntol <- length(tol);
    sel.vars <- list()
    cv_results_all <- list()
    nama <- paste("tol=", tol, sep = "")
    
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
    } else if(task == 'C') {
      
      #Classification task (logistic regression)
      if ( is.null(metric) ) {
        metricFunction <- auc.mxm;
      } else   metricFunction <- metric;
      
      if ( is.null(modeler) ) {
        modelerFunction <- glm.mxm;
      } else   modelerFunction <- modeler;
      
      if ( is.null(test) ) {
        test <- 'testIndLogistic';
      } else  test <- test;
      
    } else if (task == 'R') {
      
      #Regression task (linear regression)
      if ( is.null(metric) ) {
        metricFunction <- mse.mxm;
      } else  metricFunction <- metric;
      
      if ( is.null(modeler) ) {
        modelerFunction <- lm.mxm;
      } else  modelerFunction <- modeler;
      
      if ( is.null(test) ) {
        test = 'testIndFisher';
      } else  test <- test;
      
    } else if(task == 'S') {
      
      #cox survival analysis (cox regression)
      if ( is.null(metric) ) {
        metricFunction <- ci.mxm;
      } else  metricFunction <- metric;
      
      if ( is.null(modeler) ) {
        modelerFunction <- coxph.mxm;
      } else  modelerFunction <- modeler;
      
      if ( is.null(test) ) {
        test <- "censIndCR";
      } else  test <- test;
      
    } else  stop("Please provide a valid task argument 'C'-classification, 'R'-regression, 'S'-survival.")
    
    tic <- proc.time()
    
    for (k in 1:kfolds) {
      #print(paste('CV: Fold', k, 'of', kfolds));
      train_samples <- c();
      for ( i in which(c(1:kfolds) != k) )  train_samples = c( train_samples, folds[[ i ]] ) 
      #leave one fold out each time as a test set and the rest as train set
      train_set <- dataset[train_samples, ] #Set the training set
      train_target <- target[train_samples]
      test_set <- dataset[ folds[[ k ]], ] #Set the validation set
      test_target <- target[ folds[[ k ]] ]
      dm <- dim(test_set)
      
      results <- gomp.path(target = train_target, dataset = train_set, tol = tol, test = test, method = method)
      sel.vars[[ k ]] <-  results$res[-1, -(ntol + 1)]
      cv_results_all[[ k ]] <- list()
      cv_results_all[[ k ]]$preds <- matrix(0, dm[1], ntol )
      colnames(cv_results_all[[ k ]]$preds) <- nama
      cv_results_all[[ k ]]$performances <- numeric(ntol)
      names(cv_results_all[[ k ]]$performances) <- nama
      cv_results_all[[ k ]]$selectedVars <- sel.vars[[ k ]]
      colnames(cv_results_all[[ k ]]$selectedVars) <- nama
      
      for ( j in 1:ntol ) {
        
        variables <- sel.vars[[ k ]][, j]
        sign_data <- train_set[, variables, drop = FALSE]
        sign_test <- test_set[, variables, drop = FALSE]
        
        if ( sum( variables > 0 ) > 0 ) {
          #generate a model due to the task and find the performance
          #logistic model for a classification task, linear model for the regression task and a cox model for the survival task
          moda <- modelerFunction(train_target, sign_data, sign_test, wei = NULL)
          preds <- moda$preds
          theta <- moda$theta
        } else  {
          moda <- modelerFunction(train_target, rep(1, nrow(sign_data)), rep(1, nrow(sign_test)), wei = NULL)
          preds <- moda$preds
          theta <- moda$theta
        }
        if ( is.null(preds) ) {
          cv_results_all[[ k  ]]$preds[, j] <- NULL
          cv_results_all[[ k ]]$performances[j] <- NA
        } else {
          performance <- metricFunction(preds, test_target, theta)
          cv_results_all[[ k ]]$preds[, j] <- preds
          cv_results_all[[ k ]]$performances[j] <- performance
        }
      }  ##  end for ( i in 1:ntol ) {  
      
    }  ## end for (k in 1:kfolds) {
    
    bbc_best_performance <- NULL
    
    if (B > 1) {
       n <- dim(dataset)[1] 
       predictions <- cv_results_all[[ 1 ]]$preds
       for ( i in 2:kfolds )  predictions <- rbind(predictions, cv_results_all[[ i ]]$preds )
       bbc_best_performance <- MXM::bbc(predictions, target[unlist(folds)], metric = metricbbc, B = B )$bbc.perf
    }
    
    runtime <- proc.time() - tic
    perf <- matrix(0, nrow = kfolds, ncol = ntol)
    for (i in 1:kfolds)  perf[i, ] <- cv_results_all[[ i ]]$performances
    perf <- colMeans(perf, na.rm = TRUE)
    best_performance <- max(perf)
    best_configuration <- tol[ which.max(perf) ]
    list(cv_results_all = cv_results_all, best_performance = best_performance, best_configuration = best_configuration, 
         bbc_best_performance = bbc_best_performance, runtime = runtime) 
}














