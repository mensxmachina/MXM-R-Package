wald.univregs <- function(target, dataset, targetID = - 1, test = NULL, user_test = NULL, wei = NULL, ncores = 1) {
  
  la <- length( unique(target) )
  univariateModels <- list();
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]
  if (targetID != -1 )   {
    target <- dataset[, targetID]
    dataset[, targetID] <- rnorm(length(target) )
  }  
  ina <- NULL
  id <- NULL
  id <- Rfast::check_data(dataset)
  if ( sum(id > 0) )  dataset[, id] <- rnorm(rows * length(id) )
  
  if ( !is.null(user_test) ) {
    univariateModels <- wald.univariateScore(target, dataset, test = user_test, wei, targetID)
    
  } else if ( identical(test, waldBeta) ) {  ## Beta regression
    mod <- MXM::wald.betaregs(target, dataset, wei, logged = TRUE, ncores = ncores)
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]

  } else if ( identical(test, waldMMReg) ) {  ## M (Robust) linear regression
    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- numeric(cols)
      for ( i in 1:cols ) {
        fit <- MASS::rlm(target ~ dataset[, i], maxit = 2000, method = "MM" )
        stat[i] <- summary(fit)[[ 4 ]][2, 3]^2
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)

    } else {
      
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = 1:cols, .combine = rbind, .packages = "MASS") %dopar% {
        fit = MASS::rlm(target ~ dataset[, i], maxit = 2000, method = "MM" )
        res <- summary(fit)[[ 4 ]]
        return( res[2, 3]^2 )
      } 
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }   

  } else if ( identical(test, waldLogistic) ) {  ## logistic regression
    
    mod <- MXM::wald.logisticregs(target, dataset, wei = wei, logged = TRUE) 
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]

  } else if ( identical(test, waldOrdinal) ) {  ## Zero-inflated Poisson regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- numeric(cols)
      for ( i in 1:cols ) {
        fit <- ordinal::clm(target ~ dataset[, i], weights = wei)
        stat[i] <- summary(fit)[[ 5 ]][2, 3]^2
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = 1:cols, .combine = rbind, .packages = "ordinal") %dopar% {
        fit = ordinal::clm(target ~ dataset[, i], weights = wei)
        return( summary(fit)[[ 5 ]][2, 3]^2 )
      }
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }

  } else if ( identical(test, waldBinom)  ) {  ## Logistic regression
    wei <- target[, 2] 
    y <- target[, 1] / wei

    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- numeric(cols)
      for ( i in 1:cols ) {
        fit <- glm( y ~ dataset[, i], binomial, weights = wei )
        stat[i] <- summary(fit)[[ 12 ]][2, 3]^2
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)

    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      wei <- target[, 2] 
      y <- target[, 1] / wei
      mod <- foreach::foreach(i = 1:cols, .combine = rbind) %dopar% {
          fit <- glm( y ~ dataset[, i], binomial, weights = wei )
          return( summary(fit)[[ 12 ]][2, 3]^2 )
      }
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, waldPois) ) {  ## Poisson regression
    mod <- wald.poissonregs( target, dataset, wei = wei, logged = TRUE ) 
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]

  } else if ( identical(test, waldGamma) ) {  ## Gamma regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      stat = numeric(cols)
      for ( i in 1:cols ) {
        fit = glm( target ~ dataset[, i], family = Gamma(link = log), weights = wei )
        stat[i] = summary(fit)[[ 12 ]][2, 3]^2 / summary(fit)[[ 14 ]]
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = 1:cols, .combine = rbind) %dopar% {
        fit <- glm( target ~ dataset[, i], family = Gamma(link = log), weights = wei )
        return( summary(fit)[[ 12 ]]  [2, 3]^2 / summary(fit)[[ 14 ]] )
      }
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, waldNormLog) ) {      ## Gaussian regression with a log link
    
    if ( ncores <= 1 | is.null(ncores) ) {
      stat = numeric(cols)
      for ( i in 1:cols ) {
        fit = glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
        stat[i] = summary(fit)[[ 12 ]][2, 3]^2 / summary(fit)[[ 14 ]]
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = 1:cols, .combine = rbind) %dopar% {
        fit <- glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
        return( summary(fit)[[ 12 ]]  [2, 3]^2 / summary(fit)[[ 14 ]] )
      }
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, waldNB) ) {  ## Zero-inflated Poisson regression

    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- numeric(cols)
      for ( i in 1:cols ) {
        fit = MASS::glm.nb( target ~ dataset[, i], weights = wei )
        stat[i] = summary(fit)[[ 11 ]][2, 3]^2
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)

    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = 1:cols, .combine = rbind, .packages = "MASS") %dopar% {
          fit = MASS::glm.nb( target ~ dataset[, i], weights = wei )
          return( summary(fit)[[ 11 ]][2, 3]^2 )
      }
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, waldZIP) ) {  ## Zero-inflated Poisson regression
    mod <- wald.zipregs(target, dataset, wei, logged = TRUE, ncores = ncores) 
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]

  } else if ( identical(test, waldIGreg) ) {  ## Poisson regression

    if ( ncores <= 1 | is.null(ncores) ) {
      stat = numeric(cols)
      for ( i in 1:cols ) {
        fit = glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
        stat[i] = summary(fit)[[ 12 ]][2, 3]^2
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)

    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = 1:cols, .combine = rbind) %dopar% {
          fit <- glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
          return( summary(fit)[[ 12 ]]  [2, 3]^2 )
      }
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, waldCR) ) {  ## Cox regression

    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- numeric(cols)
      for ( i in 1:cols) {
         fit = survival::coxph( target ~ dataset[, i], weights = wei)
         stat[i] = fit$wald.test
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)

    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
          fit <- survival::coxph( target ~ dataset[, i], weights = wei)
          return( fit$wald.test )
      }
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, waldWR) ) {  ## Weibull regression

    if ( ncores <= 1 | is.null(ncores) ) {
      stat = numeric(cols)
      for ( i in 1:cols ) {
        fit = survival::survreg( target ~ dataset[, i], weights = wei, control = list(iter.max = 5000) )
        stat[i] = summary(fit)[[ 9 ]][2, 3]^2
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)

    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
          fit = survival::survreg( target ~ dataset[, i], weights = wei, control = list(iter.max = 5000) )
          return( summary(fit)[[ 9 ]][2, 3]^2 )
      }
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, waldER) ) {  ## Exponential regression

    if ( ncores <= 1 | is.null(ncores) ) {
      stat = numeric(cols)
      for ( i in 1:cols ) {
        fit = survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
        stat[i] = summary(fit)[[ 9 ]][2, 3]^2
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)

    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
          fit = survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
          return( summary(fit)[[ 9 ]][2, 3]^2 )
      }
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, waldLLR) ) {  ## Log-logistic regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      stat = numeric(cols)
      for ( i in 1:cols ) {
        fit = survival::survreg( target ~ dataset[, i], weights = wei, control = list(iter.max = 5000), dist = "loglogistic" )
        stat[i] = summary(fit)[[ 9 ]][2, 3]^2
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit = survival::survreg( target ~ dataset[, i], weights = wei, control = list(iter.max = 5000), dist = "loglogistic" )
        return( summary(fit)[[ 9 ]][2, 3]^2 )
      }
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, waldTobit) ) {  ## Tobit regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- numeric(cols)
      for ( i in 1:cols ) {
        fit <- survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
        stat[i] <- summary(fit)[[ 9 ]][2, 3]^2
      }
      univariateModels$stat <- stat
      univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit <- survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
        return( summary(fit)[[ 9 ]][2, 3]^2 )
      }
      parallel::stopCluster(cl)
      univariateModels$stat <- mod
      univariateModels$pvalue <- pchisq(mod, 1, lower.tail = FALSE, log.p = TRUE)
    }   
    
  } else   univariateModels <- NULL
  
  if ( !is.null(univariateModels) )  {
    if (targetID != - 1) {
      univariateModels$stat[targetID] <- 0
      univariateModels$pvalue[targetID] <- log(1)
    }
    if ( sum(id>0) > 0 ) {
      univariateModels$stat[id] <- 0
      univariateModels$pvalue[id] <- log(1)
    }
  }
  
  univariateModels
}