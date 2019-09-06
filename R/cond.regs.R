cond.regs <- function(target, dataset, xIndex, csIndex, test = NULL, wei = NULL, ncores = 1) {
  
  if ( identical(csIndex, 0) ) {
    models <- MXM::univregs(target = target, dataset = dataset, test = test, wei = wei, ncores = ncores)  
  } else {
  
  if ( length(xIndex) > 0 ) {
      
  models <- list()
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]
  id <- NULL
  # id <- Rfast::check_data(dataset)
  # if ( sum(id > 0) )  dataset[, id] <- rnorm(rows * length(id) )
  
  oop <- options(warn = -1) 
  on.exit( options(oop) )
	
  if ( identical(test, testIndBeta) ) {  ## Beta regression
    lik2 <- dof <- numeric( cols )
    fit1 <- beta.reg(target, dataset[, csIndex, drop = FALSE])
    lik1 <- fit1$loglik
    d1 <- length(fit1$be)
    for (i in xIndex) {
      fit2 <- beta.reg(target, dataset[, c(csIndex, xIndex[i] ) ])
      lik2[i] <- fit2$loglik
      dof[i] <- length(fit2$be)
    }
    models$stat <- 2 * (lik2 - lik1)
    models$pvalue <- pchisq(stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
    
  } else if ( identical(test, testIndMMReg) ) {  ## M (Robust) linear regression
    fit1 <- MASS::rlm( target ~., data = dataset[, csIndex, drop = FALSE], maxit = 2000, method = "MM" )
    lik1 <- as.numeric( logLik(fit1) )
    d1 <- length(coef(fit1))
    
    if ( ncores <= 1 | is.null(ncores) ) {
      lik2 <- dof <- numeric( cols )
      for ( i in xIndex ) {
        fit2 <- MASS::rlm( target ~., data = dataset[, c(csIndex, i)], maxit = 2000, method = "MM" )
        lik2[i] <- as.numeric( logLik(fit2) )
        dof[i] <- length( coef(fit2) ) 
      } 
      models$stat <- 2 * (lik2 - lik1)
      models$pvalue <- pchisq(models$stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "MASS") %dopar% {
        fit2 <- MASS::rlm( target ~., data = dataset[, c(csIndex, i)], maxit = 2000, method = "MM" )
        lik2 <- as.numeric( logLik(fit2) )
        return( c(lik2, length( coef(fit2) ) ) )
      }  
      parallel::stopCluster(cl)
      models$stat <- as.vector( 2 * ( mod[, 1] - lik1) )
      models$pvalue <- pchisq(models$stat, mod[, 2] - d1, lower.tail = FALSE, log.p = TRUE)
    }   
    
  } else if ( identical(test, testIndReg) ) {  ## linear regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- pvalue <- numeric(cols)
      fit1 <- lm( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei, y = FALSE, model = FALSE )
      
      for ( i in xIndex) {
        fit2 <- lm( target ~., data = dataset[, c(csIndex, i)], weights = wei, y = FALSE, model = FALSE )
        tab <- anova(fit1, fit2)
        stat[i] <- tab[2, 5] 
        df1 <- tab[2, 3]    
        df2 <- tab[2, 1]
        pvalue[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
      }
      models$stat <- stat
      models$pvalue <- pvalue
            
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      fit1 <- lm( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei, y = FALSE, model = FALSE )
      mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
        ww <- lm( target ~., data = dataset[, c(csIndex, i)], weights = wei, y = FALSE, model = FALSE )
        tab <- anova( fit1, ww )
        stat <- tab[2, 5] 
        df1 <- tab[2, 3]  
        df2 = tab[2, 1]
        pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE )
        return( c(stat, pval) )
      }
      parallel::stopCluster(cl)
      models$stat <- mod[, 1]
      models$pvalue <- mod[, 2]
    }   

  } else if ( identical(test, testIndMVreg)  &  !is.null(wei) ) {  ## Weighted linear regression
  
    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- pval <- numeric(cols)
      fit1 <- lm( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei, y = FALSE, model = FALSE )
      
      for ( i in xIndex ) {
        fit2 <- lm( target ~., data = dataset[, c(csIndex, i)], weights = wei, y = FALSE, model = FALSE )
        tab <- anova(fit1, fit2)
        stat[i] <- tab[2, 5] 
        df1 <- tab[2, 6]   
        df2 <- tab[2, 7]
        pvalue[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
      }
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      fit1 <- lm( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei, y = FALSE, model = FALSE )
      
      mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
        ww <- lm( target ~., data = dataset[, c(csIndex, i)], weights = wei, y = FALSE, model = FALSE )
        tab <- anova( fit1, ww )
        stat <- tab[2, 5] 
        df1 <- tab[2, 6]   
        df2 <- tab[2, 7]
        pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE )
        return( c(stat, pval) )
      }
      parallel::stopCluster(cl)
      stat <- mod[, 1]
      pvalue <- mod[, 2]
    }   
    
  } else if ( identical(test, testIndOrdinal) ) {  ## ordinal regression
    fit1 <- ordinal::clm( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei )
    lik1 <- as.numeric( logLik(fit1) )
    d1 <- length( coef(fit1) )
    
    if ( ncores <= 1 | is.null(ncores) ) {
      lik2 <- dof <- numeric( cols )
      for ( i in xIndex ) {
        fit2 <- ordinal::clm( target ~., data = dataset[, c(csIndex, i)], weights = wei )
        lik2[i] <- as.numeric( logLik(fit2) )
        dof[i] <- length( coef(fit2) ) - d1
      }
      models$stat <- 2 * (lik2 - lik1)
      models$pvalue <- pchisq(models$stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "ordinal") %dopar% {
        fit2 <- ordinal::clm( target ~., data = dataset[, c(csIndex, i)], weights = wei )
        lik2[i] <- as.numeric( logLik(fit2) )
        return( c(lik2, length( coef(fit2) ) ) )
      } 
      parallel::stopCluster(cl)
      models$stat <- 2 * (mod[, 1] - lik1)
      pvalue <- pchisq(models$stat, mod[, 2] - d1, lower.tail = FALSE, log.p = TRUE)
      
    }
    
  } else if ( identical(test, testIndMultinom) ) {  ## multinomial regression
    
    target = as.factor( as.numeric( target ) )
    fit1 <- nnet::multinom( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei, trace = FALSE )
    lik1 = as.numeric( logLik(fit1) )
    d1 = length( coef(fit1) )
    
    if ( ncores <= 1 | is.null(ncores) ) {
      lik2 <- dof <- numeric( cols )
      for ( i in xIndex ) {
        fit2 = nnet::multinom( target ~., data = dataset[, c(csIndex, i)], weights = wei, trace = FALSE )
        lik2[i] <- as.numeric( logLik(fit2) )
        dof[i] <- length( coef(fit2) ) 
      }
      models$stat <- 2 * (lik2 - lik1)
      models$pvalue <- pchisq(models$stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "nnet") %dopar% {
        fit2 = nnet::multinom( target ~., data = dataset[, c(csIndex, i)], weights = wei, trace = FALSE )
        lik2 = as.numeric( logLik(fit2 ) )
        return( c(lik2, length( coef(fit2) ) ) )
        
      }
      parallel::stopCluster(cl)
      models$stat <- 2 * (mod[, 1] - lik1)
      models$pvalue <- pchisq(models$stat, mod[, 2] - d1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, testIndLogistic) ) {  ## Logistic regression
    fit1 <- glm( target ~., data = dataset[, csIndex, drop = FALSE], binomial, weights = wei )
    lik1 <- fit1$deviance
    d1 <- length(fit1$coefficients)
    
    if ( ncores <= 1 | is.null(ncores) ) {
      lik2 <- dof <- numeric( cols )
      for ( i in xIndex ) {
        fit2 = glm( target ~., data = dataset[, c(csIndex, i)], binomial, weights = wei )
        lik2[i] <- fit2$deviance
        dof[i] <- length( fit2$coefficients ) 
      }
      models$stat <- lik1 - lik2
      models$pvalue <- pchisq(models$stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
        fit2 = glm( target ~., data = dataset[, c(csIndex, i)], binomial, weights = wei )
        lik2 = fit2$deviance
        return( c(lik2, length( fit2$coefficients ) ) )
      }
      parallel::stopCluster(cl)
      models$stat <- lik1 - mod[, 1]
      models$pvalue <- pchisq(models$stat, mod[, 2] - d1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, testIndBinom) ) {  ## Binomial regression
    wei <- target[, 2] 
    y <- target[, 1] / wei
    fit1 <- glm( y ~., data = dataset[, csIndex, drop = FALSE], binomial, weights = wei )
    lik1 <- fit1$deviance
    d1 <- length(fit1$coefficients)
    
    if ( ncores <= 1 | is.null(ncores) ) {
      lik2 <- dof <- numeric( cols )
      for ( i in xIndex ) {
        fit2 <- glm( y ~., data = dataset[, c(csIndex, i)], binomial, weights = wei )
        lik2[i] <- fit2$deviance
        dof[i] <- length( coef(fit2) ) 
      }
      models$stat <- lik1 - lik2
      models$pvalue <- pchisq(models$stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      wei <- target[, 2] 
      y <- target[, 1] / wei
      mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
        fit2 = glm( y ~., data = dataset[, c(csIndex, i)], binomial, weights = wei )
        lik2 = as.numeric( logLik(fit2) )
        return( c(lik2, length( fit2$coefficients ) ) )
      }
      parallel::stopCluster(cl)
      models$stat <- as.vector( lik1 - mod[, 1] )
      models$pvalue <- pchisq(models$stat, mod[, 2] - d1, lower.tail = FALSE, log.p = TRUE)
    }

  } else if ( identical(test, testIndPois) ) {  ## Poisson regression
    fit1 <- glm( target ~., data = dataset[, csIndex, drop = FALSE], poisson, weights = wei )
    lik1 = fit1$deviance
    d1 <- length(fit1$coefficients)
    
    if ( ncores <= 1 | is.null(ncores) ) {
      lik2 <- dof <- numeric( cols )
      for ( i in xIndex ) {
        fit2 = glm( target ~., data = dataset[, c(csIndex, i)], poisson, weights = wei )
        lik2[i] <- fit2$deviance
        dof[i] <- length( fit2$coefficients ) 
      }
      models$stat <- lik1 - lik2
      models$pvalue <- pchisq(models$stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
        fit2 = glm( target ~., data = dataset[, c(csIndex, i)], poisson, weights = wei )
        return( c(fit2$deviance, length( fit2$coefficients ) ) )
      }
      parallel::stopCluster(cl)
      models$stat <- as.vector( lik1 - mod[, 1] )
      models$pvalue <- pchisq(models$stat, mod[, 2] - d1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, testIndNB) ) {  ## Negative binomial regression
    lik1 <- MASS::glm.nb( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei )$twologlik
    d1 <- length(fit1$coefficients)
    
    if ( ncores <= 1 | is.null(ncores) ) {
      lik2 <- dof <- numeric( cols )
      for ( i in xIndex ) {
        fit2 = MASS::glm.nb( target ~., data = dataset[, c(csIndex, i)], weights = wei )
        lik2[i] <- fit2$twologlik
        dof[i] <- length( fit2$coefficients ) 
      }
      models$stat <- lik2 - lik1
      models$pvalue <- pchisq(models$stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "MASS") %dopar% {
        fit2 = MASS::glm.nb( target ~., data = dataset[, c(csIndex, i)], weights = wei )
        return( c(fit2$twologlik, length( fit2$coefficients ) ) )
      }
      parallel::stopCluster(cl)
      models$stat <- as.vector(mod[, 1]) - lik1
      models$pvalue <- pchisq(models$stat, mod[, 2] - d1, lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, testIndQBinom)   ) {  ## Quasi Binomial regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      fit1 <- glm( target ~ ., data = dataset[, csIndex, drop = FALSE], family = quasibinomial(link = logit), weights = wei )
      for ( i in xIndex ) {
        fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = quasibinomial(link = logit), weights = wei )
        tab <- anova(fit1, fit2, test = "F")
        stat[i] <- tab[2, 5]
        df1 <- tab[2, 3]   
        df2 <- tab[2, 1]
        pvalue[i] <- pf(tab[2,5] ,df1,df2, lower.tail = FALSE, log.p = TRUE)
      }
      models$stat <- stat
      models$pvalue <- pvalue
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      fit1 <- glm( target ~ ., data = dataset[, csIndex, drop = FALSE], family = quasibinomial(link = logit), weights = wei )
      mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
        fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = quasibinomial(link = logit), weights = wei )
        tab <- anova(fit1, fit2, test = "F")
        stat <- tab[2, 5]
        df1 <- tab[2, 3]   
        df2 <- tab[2, 1]
        pvalue <- pf(tab[2, 5], df1, df2, lower.tail = FALSE, log.p = TRUE)
        return( c(stat, pvalue ) )
      }
      parallel::stopCluster(cl)
      models$stat <- mod[, 1]
      models$pvalue <- mod[, 2]
    }
    
  } else if ( identical(test, testIndQPois)   ) {  ## Quasi Poisson regression
    stat <- pvalue <- numeric(cols)
    fit1 <- glm( target ~ ., data = dataset[, csIndex, drop = FALSE], family = quasipoisson(link = log), weights = wei )
    
    if ( ncores <= 1 | is.null(ncores) ) {
      for ( i in xIndex ) {
        fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = quasipoisson(link = log), weights = wei )
        tab <- anova(fit1, fit2, test = "F")
        stat[i] <- tab[2, 5]
        df1 <- tab[2, 3]   
        df2 <- tab[2, 1]
        pvalue[i] <- pf(tab[2,5] ,df1,df2, lower.tail = FALSE, log.p = TRUE)
      }
      models$stat <- stat
      models$pvalue <- pvalue
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      fit1 <- glm( target ~ ., data = dataset[, csIndex, drop = FALSE], family = quasipoisson(link = log), weights = wei )
      mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
        fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = quasipoisson(link = log), weights = wei )
        tab <- anova(fit1, fit2, test = "F")
        stat <- tab[2, 5]
        df1 <- tab[2, 3]   
        df2 <- tab[2, 1]
        pvalue <- pf(tab[2, 5], df1, df2, lower.tail = FALSE, log.p = TRUE)
        return( c(stat, pvalue ) )
      }
      parallel::stopCluster(cl)
      models$stat <- mod[, 1]
      models$pvalue <- mod[, 2]
    }
    
  } else if ( identical(test, testIndNormLog) ) {  ## Normal log link regression
    fit1 <- glm(target ~., data = dataset[, csIndex, drop = FALSE], family = gaussian(link = log), weights = wei)

    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- pvalue <- numeric( cols )
      for ( i in xIndex ) {
        fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = gaussian(link = log), weights = wei )
        tab <- anova(fit1, fit2, test = "F")
        stat[i] <- tab[2, 5]
        df1 <- tab[2, 3]   
        df2 <- tab[2, 1]
        pvalue[i] <- pf(tab[2,5] ,df1,df2, lower.tail = FALSE, log.p = TRUE)
      }
      models$stat <- stat
      models$pvalue <- pvalue
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
        fit2 <- glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
        tab <- anova(fit1, fit2, test = "F")
        stat <- tab[2, 5]
        df1 <- tab[2, 3]   
        df2 <- tab[2, 1]
        pvalue <- pf(tab[2, 5], df1, df2, lower.tail = FALSE, log.p = TRUE)
        return( c(stat, pvalue ) )
      }
      parallel::stopCluster(cl)
      models$stat <- mod[, 1]
      models$pvalue <- mod[, 2]
    }
    
  } else if ( identical(test, testIndGamma)   ) {  ## Gamma regression
    fit1 <- glm(target ~., data = dataset[, csIndex, drop = FALSE], family = Gamma(link = log), weights = wei)

    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- pvalue <- numeric( cols )
      for ( i in xIndex ) {
        fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = Gamma(link = log), weights = wei )
        tab <- anova(fit1, fit2, test = "F")
        stat[i] <- tab[2, 5]
        df1 <- tab[2, 3]   
        df2 <- tab[2, 1]
        pvalue[i] <- pf(tab[2,5] ,df1,df2, lower.tail = FALSE, log.p = TRUE)
      }
      models$stat <- stat
      models$pvalue <- pvalue
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
        fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = Gamma(link = log), weights = wei )
        tab <- anova(fit1, fit2, test = "F")
        stat <- tab[2, 5]
        df1 <- tab[2, 3] 
        df2 <- tab[2, 1]
        pvalue <- pf(tab[2, 5], df1, df2, lower.tail = FALSE, log.p = TRUE)
        return( c(stat, pvalue ) )
      }
      parallel::stopCluster(cl)
      models$stat <- mod[, 1]
      models$pvalue <- mod[, 2]
    } 
    
  } else if ( identical(test, testIndIGreg) ) {  ## Inverse Gaussian regression
    fit1 <- glm( target ~., data = dataset[, csIndex, drop = FALSE], family = inverse.gaussian(link = log), weights = wei )

    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- pvalue <- numeric( cols )
      for ( i in xIndex ) {
        fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = inverse.gaussian(link = log), weights = wei )
        tab <- anova(fit1, fit2, test = "F")
        stat[i] <- tab[2, 5]
        df1 <- tab[2, 3] 
        df2 <- tab[2, 1]
        pvalue[i] <- pf(tab[2,5] ,df1,df2, lower.tail = FALSE, log.p = TRUE)
      }
      models$stat <- stat
      models$pvalue <- pvalue
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
        fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = inverse.gaussian(link = log), weights = wei )
        tab <- anova(fit1, fit2, test = "F")
        stat <- tab[2, 5]
        df1 <- tab[2, 3] 
        df2 <- tab[2, 1]
        pvalue <- pf(tab[2, 5], df1, df2, lower.tail = FALSE, log.p = TRUE)
        return( c(stat, pvalue ) )
      }
      parallel::stopCluster(cl)
      models$stat <- mod[, 1]
      models$pvalue <- mod[, 2]
    }
    
  } else if ( identical(test, testIndZIP) ) {  ## Zero-inflated Poisson regression
    lik2 <- dof <- numeric( cols )
    fit1 <- zip.reg(target, dataset[, csIndex, drop = FALSE])
    lik1 <- fit1$loglik
    d1 <- length(fit1$be)
    for (i in xIndex) {
      fit2 <- zip.reg(target, dataset[, c(csIndex, xIndex[i] ) ])
      lik2[i] <- fit2$loglik
      dof[i] <- length(fit2$be)
    }
    models$stat <- 2 * (lik2 - lik1)
    models$pvalue <- pchisq(stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
    
  } else if ( identical(test, testIndRQ) ) {  ## Median (quantile) regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      fit1 <- quantreg::rq(target ~., data = dataset[, csIndex, drop = FALSE], weights = wei)
      stat <- pval <- numeric(cols)
      for ( i in xIndex ) {
        fit2 = quantreg::rq(target ~., data = dataset[, c(csIndex, i)], weights = wei )
        ww = anova(fit1, fit2, test = "rank")
        df1 = as.numeric( ww[[1]][1] )
        df2 = as.numeric( ww[[1]][2] )
        stat[i] = as.numeric( ww[[1]][3] )
        pval[i] = pf(stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE)
      }
      models$stat <- stat
      models$pvalue <- pval
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      fit1 <- quantreg::rq(target ~., data = dataset[, csIndex, drop = FALSE], weights = wei)
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "quantreg") %dopar% {
        fit2 = quantreg::rq(target ~., data = dataset[, c(csIndex, i)], weights = wei )
        ww = anova(fit1, fit2, test = "rank")
        df1 = as.numeric( ww[[1]][1] )
        df2 = as.numeric( ww[[1]][2] )
        stat = as.numeric( ww[[1]][3] )
        pval = pf(stat, df1, df2, lower.tail = FALSE, log.p = TRUE)
        return( c(stat, pval ) )
      }
      parallel::stopCluster(cl)
      models$stat <- as.vector( mod[, 1] )
      models$pvalue <- as.vector( mod[, 2] )
    }
    
  } else if ( identical(test, censIndCR) ) {  ## Cox regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- dof <- numeric(cols)
      fit1 <- survival::coxph( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei)
      for ( i in xIndex ) {
        fit2 = survival::coxph( target ~., data = dataset[, c(csIndex, i)], weights = wei)
        res <- anova(fit1, fit2)
        stat[i] <- res[2, 2]
        dof[i] <- res[2, 3]
      }
      models$stat <- stat
      models$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      fit1 <- survival::coxph( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei)
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::coxph( target ~., data = dataset[, c(csIndex, i)], weights = wei )
        res <- anova(fit2)
        return( c(res[2, 2], res[2, 3] ) )
      }
      parallel::stopCluster(cl)
      models$stat <- mod[, 1]
      models$pvalue <- pchisq(mod[, 1], mod[, 2], lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, censIndWR) ) {  ## Weibull regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- dof <- numeric(cols)
      fit1 <- survival::survreg( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei )
      for ( i in xIndex ) {
        fit2 = survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei )
        res <- anova(fit1, fit2)
        stat[i] <- res[2, 6]
        dof[i] <- res[2, 5]
      }
      models$stat <- stat
      models$pvalue <- pchisq(models$stat, dof, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      fit1 <- survival::survreg( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei )
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei )
        res <- anova(fit1, fit2)
        ## stat <- res[2, 6]
        ## dof <- res[2, 5]
        return( c( res[2, 6],res[2, 5] ) )
      }
      parallel::stopCluster(cl)
      models$stat <- mod[, 1]
      models$pvalue <- pchisq(models$stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, censIndER) ) {  ## Exponential regression
    if ( ncores <= 1 | is.null(ncores) ) {
      fit1 <- survival::survreg( target ~., data = dataset[, csIndex, drop = FALSE], dist = "exponential", weights = wei )
      stat <- dof <- numeric(cols)
      
      for ( i in xIndex ) {
        fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], dist = "exponential", weights = wei )
        res <- anova(fit1, fit2)
        stat[i] <- res[2, 6]
        dof[i] <- res[2, 5]
      }
      models$stat <- stat
      models$pvalue <- pchisq(models$stat, dof, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      fit1 <- survival::survreg( target ~., data = dataset[, csIndex, drop = FALSE], dist = "exponential", weights = wei )
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], dist = "exponential", weights = wei )
        res <- anova(fit1, fit2)
        ## stat <- res[2, 6]
        ## dof <- res[2, 5]
        return( c( res[2, 6],res[2, 5] ) )
      }
      parallel::stopCluster(cl)
      models$stat <- mod[, 1]
      models$pvalue <- pchisq(models$stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, censIndLLR) ) {  ## Log-logistic regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      stat <- dof <- numeric(cols)
      fit1 <- survival::survreg( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei, dist = "loglogistic" )
      for ( i in xIndex ) {
        fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei, dist = "loglogistic" )
        res <- anova(fit1, fit2)
        stat[i] <- res[2, 6]
        dof[i] <- res[2, 5]
      }
      models$stat <- stat
      models$pvalue <- pchisq(models$stat, dof, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      fit1 <- survival::survreg( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei, dist = "loglogistic" )
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei, dist = "loglogistic" )
        res <- anova(fit1, fit2)
        ## stat <- res[2, 6]
        ## dof <- res[2, 5]
        return( c( res[2, 6],res[2, 5] ) )
      }
      parallel::stopCluster(cl)
      models$stat <- mod[, 1]
      models$pvalue <- pchisq(models$stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, testIndTobit) ) {  ## Tobit regression

    if ( ncores <= 1 | is.null(ncores) ) {
      fit1 <- survival::survreg( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei, dist = "gaussian" )
      stat <- dof <- numeric(cols)
      for ( i in xIndex ) {
        fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei, dist = "gaussian" )
        res <- anova(fit1, fit2)
        stat[i] <- res[2, 6]
        dof[i] <- res[2, 5]
      }
      models$stat <- stat
      models$pvalue <- pchisq(models$stat, dof, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      fit1 <- survival::survreg( target ~., data = dataset[, csIndex, drop = FALSE], weights = wei, dist = "gaussian" )
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
        fit2 = survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei, dist = "gaussian" )
        res <- anova(fit1, fit2)
        return( c(res[2, 6], res[2, 5] ) )
      }
      parallel::stopCluster(cl)
      models$stat <- mod[, 1]
      models$pvalue <- pchisq(models$stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, testIndClogit) ) {  ## Conditional logistic regression
    subject <- target[, 2] #the patient id
    case <- as.logical(target[, 1])  ## case 

    if ( ncores <= 1 | is.null(ncores) ) {
      fit1 <- survival::clogit( case ~., data = dataset[, csIndex, drop = FALSE] + strata(subject) ) 
      stat <- dof <- numeric(cols)
      for ( i in xIndex ) {
        fit2 <- survival::clogit( case ~., data = dataset[, c(csIndex, i)] + strata(subject) ) 
        res <- anova(fit1, fit2)
        stat[i] <- res[2, 2]
        dof[i] <- res[2, 3] 
      }
      pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      fit1 <- survival::clogit( case ~., data = dataset[, csIndex, drop = FALSE] + strata(subject) ) 
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::clogit( case ~., data = dataset[, c(csIndex, i)] + strata(subject) ) 
        res <- anova(fit1, fit2)
        return( c(res[2, 2], res[2, 3]) )
      }
      parallel::stopCluster(cl)
      models$stat <- as.vector( mod[, 1] )
      models$pvalue <- pchisq(mod[, 1], mod[, 2], lower.tail = FALSE, log.p = TRUE)
    }
    
  } else if ( identical(test, testIndSPML) ) {  ## Circular regression
    if ( !is.matrix(target) )   target <- cbind( cos(target), sin(target) )
    fit1 <- Rfast::spml.reg( target, dataset[, csIndex] )
    lik1 <- 2 * fit1$loglik

    if ( ncores <= 1 | is.null(ncores) ) {
      lik2 <- numeric( cols )
      for ( i in xIndex ) {
        fit2 <- try( Rfast::spml.reg( target, dataset[, c(csIndex, i)] ), silent = TRUE )
        if ( identical( class(fit2), "try-error" ) )   {
          lik2[i] <- lik1
        } else   lik2[i] <- fit2$loglik
      }
      models$stat <- 2 * lik2 - lik1
      models$pvalue <- pchisq(models$stat, 2, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      mod <- foreach::foreach(i = xIndex, .combine = rbind, .pcakages = "Rfast", .export = "spml.reg") %dopar% {
        fit2 <- Rfast::spml.reg( target, dataset[, c(csIndex, i)] )
        return( fit2$loglik )
      }
      parallel::stopCluster(cl)
      models$stat <- 2 * as.vector(mod) - lik1
      models$pvalue <- pchisq(models$stat, 2, lower.tail = FALSE, log.p = TRUE)
    }
    
  }  else   models <- NULL  ## end of all if (test == )
  
  models$stat[ - xIndex ] <- 0
  models$pvalue[ - xIndex ] <- log(1)
  
  } else {
    models <- list()
    models$stat <- NULL
    models$pvalue <- NULL
  }  ##  end  if  ( length(xIndex) > 0 )  
    
  }  ## end if ( identical(csIndex, 0) )
  
  models
}