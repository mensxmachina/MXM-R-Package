univregs <- function(target, dataset, targetID = -1, test = NULL, user_test = NULL, wei = NULL, ncores = 1) {
  
  univariateModels <- list();
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]
  if (targetID != -1 ) {
     target <- dataset[, targetID]
     dataset[, targetID] <- rbinom(rows, 1, 0.5)
  }   
  id <- NULL
  if ( !identical(test, testIndFisher) & !identical(test, testIndSpearman) ) {
    ina <- NULL
    id <- Rfast::check_data(dataset)
    if ( sum(id > 0) )  dataset[, id] <- rnorm(rows * length(id) )
  }  
  la <- length( unique(target) )
  
if ( !is.null(user_test) ) {
  univariateModels <- univariateScore(target, dataset, test = user_test, wei, targetID)

} else if ( identical(test, testIndFisher) )  { ## Pearson's correlation 
  a <- as.vector( cor(target, dataset) )
  dof <- rows - 3; #degrees of freedom
  wa <- 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof)
  id <- which( is.na(a) )
  if ( length(id) > 0)  wa[id] <- 0
  univariateModels$stat <- wa;
  univariateModels$pvalue <- log(2) + pt( abs(wa), dof, lower.tail = FALSE, log.p = TRUE) ;

} else if ( identical(test, testIndSpearman) ) {  ## Spearman's correlation
  a <- as.vector( cor(target, dataset) )
  dof <- rows - 3; #degrees of freedom
  wa <- 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof) / 1.029563
  id <- which( is.na(a) )
  if ( length(id) > 0)  wa[id] <- 0
  univariateModels$stat <- wa 
  univariateModels$pvalue <- log(2) + pt( abs(wa), dof, lower.tail = FALSE, log.p = TRUE);

} else if ( identical(test, gSquare) ) {  ## G^2 test
  z <- cbind(dataset, target)
  if ( !is.matrix(z) )   z <- as.matrix(z)
  dc <- Rfast::colrange(z, cont = FALSE)
  a <- Rfast::g2tests(data = z, x = 1:cols, y = cols + 1, dc = dc)
  stat <- a$statistic
  univariateModels$stat <- stat
  univariateModels$pvalue <- pchisq(stat, a$df, lower.tail = FALSE, log.p = TRUE)

} else if ( identical(test, testIndBeta) ) {  ## Beta regression
  mod <- beta.regs(target, dataset, wei, logged = TRUE, ncores = ncores)
  univariateModels$stat <- mod[, 1]
  univariateModels$pvalue <- mod[, 2]

} else if ( identical(test, testIndMMReg) ) {  ## M (Robust) linear regression
  fit1 <- MASS::rlm(target ~ 1, maxit = 2000, method = "MM")
  lik1 <- as.numeric( logLik(fit1) )
  lik2 <- numeric(cols)
  dof <- numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <- MASS::rlm(target ~ dataset[, i], maxit = 2000, method = "MM" )
      lik2[i] <- as.numeric( logLik(fit2) )
      dof[i] <- length( coef(fit2) ) - 1
    } 
    stat <- 2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "MASS") %dopar% {
       fit2 <- MASS::rlm(target ~ dataset[, i], maxit = 2000, method = "MM" )
       lik2 <- as.numeric( logLik(fit2) )
       return( c(lik2, length( coef(fit2) ) ) )
    }  
    parallel::stopCluster(cl)
    stat <- as.vector( 2 * abs(lik1 - mod[, 1]) )
    dof <- as.vector( mod[, 2] ) - 1 
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
  }   
  
} else if ( identical(test, testIndReg)  &  !is.null(wei) ) {  ## Weighted linear regression
  
  univariateModels <- list();
  stat <- pval <- numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <- lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
      tab <- anova(fit2)
      stat[i] <- tab[1, 4] 
      df1 <- tab[1, 1]    ;   df2 = tab[2, 1]
      pval[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
    }
    univariateModels$stat <- stat
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
        ww <- lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
        tab <- anova( ww )
        stat <- tab[1, 4] 
        df1 <- tab[1, 1]   ;  df2 = tab[2, 1]
        pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE )
        return( c(stat, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]
  }   
  
} else if ( identical(test, testIndReg)  &  is.null(wei) ) {  ## linear regression
  mod <- Rfast::regression(dataset, target, logged = TRUE)
  univariateModels$stat <- mod[, 1]
  univariateModels$pvalue <- mod[, 2]

} else if ( identical(test, testIndMVreg) ) {  ## Weighted linear regression
  
  univariateModels = list();
  stat = pval = numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
      tab = anova(fit2)
      stat[i] = tab[2, 3] 
      df1 = tab[2, 4]    ;   df2 = tab[2, 5]
      pval[i] = pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
    }
    univariateModels$stat <- stat
    univariateModels$pvalue <- pval
    
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
      ww <- lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
      tab <- anova( ww )
      stat <- tab[2, 3] 
      df1 <- tab[2, 4]   ;  df2 = tab[2, 5]
      pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE )
      return( c(stat, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]
  }   
  
} else if ( identical(test, testIndOrdinal) ) {  ## ordinal regression
  lik2 <- numeric(cols)
  dof <- numeric(cols)
  fit1 <- ordinal::clm(target ~ 1, weights = wei)
  lik1 <- as.numeric( logLik(fit1) )
  df1 <- length( coef(fit1) )
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      mat <- model.matrix(target ~ dataset[, i] )
      fit2 <- ordinal::clm.fit(target, mat, weights = wei)
      lik2[i] <- as.numeric( fit2$logLik )
      dof[i] <- length( coef(fit2) ) - df1
    }
    stat = 2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "ordinal") %dopar% {
        mat <- model.matrix(target ~ dataset[, i] )
        fit2 <- ordinal::clm.fit(target, mat, weights = wei)
        lik2 <- as.numeric( fit2$logLik )
        return( c(lik2, length( coef(fit2) ) ) )
    } 
    parallel::stopCluster(cl)
    stat = 2 * (mod[, 1] - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, mod[, 2] - df1, lower.tail = FALSE, log.p = TRUE)

  }
  
} else if ( identical(test, testIndMultinom) ) {  ## multinomial regression
  
  target = as.factor( as.numeric( target ) )
  lik2 = numeric(cols)
  dof = numeric(cols)
  fit1 = nnet::multinom(target ~ 1, trace = FALSE, weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  df1 = length( coef(fit1) )

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = nnet::multinom(target ~ dataset[, i], trace = FALSE, weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - df1
    }
    stat = 2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "nnet") %dopar% {
        fit2 = nnet::multinom(target ~ dataset[, i], weights = wei)
        lik2 = as.numeric( logLik(fit2 ) )
        return( c(lik2, length( coef(fit2) ) ) )

    }
    parallel::stopCluster(cl)
    stat = 2 * (mod[, 1] - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, mod[, 2] - df1, lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndLogistic)  &  is.matrix(dataset)  &  is.null(wei) ) {  ## logistic regression
  if ( is.factor(target) )   target <- as.numeric(target) - 1
  mod <- Rfast::univglms( target, dataset, oiko = "binomial", logged = TRUE )
  univariateModels$stat <- mod[, 1]
  univariateModels$pvalue <- mod[, 2]
  
} else if ( identical(test, testIndLogistic)  &  is.data.frame(dataset)  &  is.null(wei) ) {  ## logistic regression
  if ( is.factor(target) )   target <- as.numeric(target) - 1
  mod <- Rfast::univglms2( target, dataset, oiko = "binomial", logged = TRUE )
  univariateModels$stat <- mod[, 1]
  univariateModels$pvalue <- mod[, 2]

} else if ( identical(test, testIndLogistic)  &  !is.null(wei) ) {  ## Logistic regression
  fit1 = glm(target ~ 1, binomial, weights = wei)
  lik1 = fit1$deviance
  lik2 = numeric(cols)
  dof = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = glm( target ~ dataset[, i], binomial, weights = wei )
      lik2[i] = fit2$deviance
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = lik1 - lik2
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
        fit2 = glm( target ~ dataset[, i], binomial, weights = wei )
        lik2 = fit2$deviance
        return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    parallel::stopCluster(cl)
    stat = lik1 - mod[, 1]
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndBinom) ) {  ## Binomial regression
  wei <- target[, 2] 
  y <- target[, 1] / wei
  fit1 = glm(y ~ 1, binomial, weights = wei)
  lik1 = fit1$deviance
  lik2 = numeric(cols)
  dof = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = glm( y ~ dataset[, i], binomial, weights = wei )
      lik2[i] = fit2$deviance
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = lik1 - lik2
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    wei <- target[, 2] 
    y <- target[, 1] / wei
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
        fit2 = glm( y ~ dataset[, i], binomial, weights = wei )
        lik2 = as.numeric( logLik(fit2) )
        return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    parallel::stopCluster(cl)
    stat = as.vector( lik1 - mod[, 1] )
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndPois)  &  is.matrix(dataset)  &  is.null(wei) ) {  ## Poisson regression
  mod <- Rfast::univglms( target, dataset, oiko = "poisson", logged = TRUE ) 
  univariateModels$stat <- mod[, 1]
  univariateModels$pvalue <- mod[, 2]
  
} else if ( identical(test, testIndPois)  &  is.data.frame(dataset)  &  is.null(wei) ) {  ## Poisson regression
  mod <- Rfast::univglms2( target, dataset, oiko = "poisson", logged = TRUE ) 
  univariateModels$stat <- mod[, 1]
  univariateModels$pvalue <- mod[, 2]

} else if ( identical(test, testIndPois)  &  !is.null(wei) ) {  ## Poisson regression
  fit1 = glm(target ~ 1, poisson, weights = wei)
  lik1 = fit1$deviance
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in ina ) {
      fit2 = glm( target ~ dataset[, i], poisson, weights = wei )
      lik2[i] = fit2$deviance
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = lik1 - lik2
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
        fit2 = glm( target ~ dataset[, i], poisson, weights = wei )
        return( c(fit2$deviance, length( coef(fit2) ) - 1 ) )
    }
    parallel::stopCluster(cl)
    stat = as.vector( lik1 - mod[, 1] )
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndNB) ) {  ## Negative binomial regression
  lik1 <- MASS::glm.nb( target ~ 1, weights = wei )$twologlik

  if ( ncores <= 1 | is.null(ncores) ) {
    lik2 <- dof <- numeric(cols)
    for ( i in 1:cols ) {
      fit2 = MASS::glm.nb( target ~ dataset[, i], weights = wei )
      lik2[i] = fit2$twologlik
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = lik2 - lik1
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "MASS") %dopar% {
        fit2 = MASS::glm.nb( target ~ dataset[, i], weights = wei )
        return( c(fit2$twologlik, length( coef(fit2) ) - 1 ) )
    }
    parallel::stopCluster(cl)
    stat <- as.vector(mod[, 1]) - lik1
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }

} else if ( identical(test, testIndNormLog) ) {  ## Normal log link regression
  fit1 = glm(target ~ 1, family = gaussian(link = log), weights = wei)
  lik1 = fit1$deviance
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols
  phi <- numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in ina ) {
      fit2 = glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
      lik2[i] = fit2$deviance
      phi[i] <- summary(fit2)[[14]] 
      dof[i] = length( fit2$coefficients )
    }
    stat = (lik1 - lik2 ) / (dof - 1) / phi 
    univariateModels$stat <- stat
    univariateModels$pvalue <- pf(stat, dof - 1, rows - dof, lower.tail = FALSE, log.p = TRUE)
    
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
      fit2 = glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
      return( c(fit2$deviance, length( fit2$coefficients ), summary(fit2)[[14]] ) )
    }
    parallel::stopCluster(cl)
    stat = as.vector( lik1 - mod[, 1] ) / (mod[, 2] - 1) / mod[, 3] 
    univariateModels$stat <- stat
    univariateModels$pvalue <- pf(stat, mod[, 2] - 1, rows - mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
   
} else if ( identical(test, testIndGamma)   ) {  ## Gamma regression
  fit1 = glm(target ~ 1, family = Gamma(link = log), weights = wei)
  lik1 = fit1$deviance
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols
  phi <- numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in ina ) {
      fit2 = glm( target ~ dataset[, i], family = Gamma(link = log), weights = wei )
      lik2[i] = fit2$deviance
      phi[i] = summary(fit2)[[ 14 ]]
      dof[i] = length( fit2$coefficients)
    }
    stat = (lik1 - lik2) / (dof - 1) / phi
    univariateModels$stat <- stat
    univariateModels$pvalue <- pf(stat, dof - 1, rows - dof, lower.tail = FALSE, log.p = TRUE)
    
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
      fit2 = glm( target ~ dataset[, i], family = Gamma(link = log), weights = wei )
      return( c(fit2$deviance, length( fit2$coefficients ), summary(fit2)[[14]] ) )
    }
    parallel::stopCluster(cl)
    stat = as.vector( lik1 - mod[, 1] ) / (mod[, 2] - 1) / mod[, 3] 
    univariateModels$stat <- stat
    univariateModels$pvalue <- pf(stat, mod[, 2] - 1, rows - mod[, 2], lower.tail = FALSE, log.p = TRUE)
  } 
  
} else if ( identical(test, testIndZIP) ) {  ## Zero-inflated Poisson regression
  moda <- zip.regs(target, dataset, wei, logged = TRUE, ncores = ncores) 
  univariateModels$stat <- moda[, 1]
  univariateModels$pvalue <- moda[, 2]

} else if ( identical(test, testIndRQ) ) {  ## Median (quantile) regression
  
  fit1 = quantreg::rq(target ~ 1, weights = wei)
  stat = pval = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = quantreg::rq(target ~ dataset[, i], weights = wei )
      ww = anova(fit1, fit2, test = "rank")
      df1 = as.numeric( ww[[1]][1] )
      df2 = as.numeric( ww[[1]][2] )
      stat[i] = as.numeric( ww[[1]][3] )
      pval[i] = pf(stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE)
    }
    univariateModels$stat <- stat
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "quantreg") %dopar% {
        fit2 = quantreg::rq(target ~ dataset[, i], weights = wei )
        ww = anova(fit1, fit2, test = "rank")
        df1 = as.numeric( ww[[1]][1] )
        df2 = as.numeric( ww[[1]][2] )
        stat = as.numeric( ww[[1]][3] )
        pval = pf(stat, df1, df2, lower.tail = FALSE, log.p = TRUE)
        return( c(stat, pval ) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- as.vector( mod[, 1] )
    univariateModels$pvalue <- as.vector( mod[, 2] )
  }
  
} else if ( identical(test, testIndIGreg) ) {  ## Inverse Gaussian regression
  fit1 = glm(target ~ 1, family = inverse.gaussian(link = log), weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
      for ( i in 1:cols ) {
        fit2 = glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
        lik2[i] = as.numeric( logLik(fit2) )
        dof[i] = length( coef(fit2) ) - 1
      }
    stat = 2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
        fit2 = glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
        lik2 = as.numeric( logLik(fit2) )
        return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    parallel::stopCluster(cl)
    stat = as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, censIndCR) ) {  ## Cox regression
  stat = numeric(cols)
  dof = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = survival::coxph( target ~ dataset[, i], weights = wei)
      res <- anova(fit2)
      dof[i] <- res[2, 3]
      stat[i] <- res[2, 2]
    }
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 = survival::coxph( target ~ dataset[, i], weights = wei )
        res <- anova(fit2)
        return( c(res[2, 2], res[2, 3] ) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- pchisq(mod[, 1], mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, censIndWR) ) {  ## Weibull regression
  fit1 <- survival::survreg(target ~ 1, weights = wei)
  lik1 <- as.numeric( logLik(fit1) )
  lik2 <- numeric(cols)
  dof <- numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <- survival::survreg( target ~ dataset[, i], weights = wei, control=list(iter.max = 5000)  )
      lik2[i] <- as.numeric( logLik(fit2) )
      dof[i] <- length( coef(fit2) ) - 1
    }
    stat <- 2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::survreg( target ~ dataset[, i], weights = wei, control=list(iter.max = 5000)  )
        lik2 <- as.numeric( logLik(fit2) )
        return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    parallel::stopCluster(cl)
    stat <- as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, censIndLLR) ) {  ## Weibull regression
  fit1 <- survival::survreg(target ~ 1, weights = wei, dist = "loglogistic")
  lik1 <- as.numeric( logLik(fit1) )
  lik2 <- numeric(cols)
  dof <- numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <- survival::survreg( target ~ dataset[, i], weights = wei, control = list(iter.max = 5000), dist = "loglogistic" )
      lik2[i] <- as.numeric( logLik(fit2) )
      dof[i] <- length( coef(fit2) ) - 1
    }
    stat <- 2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
      fit2 <- survival::survreg( target ~ dataset[, i], weights = wei, control = list(iter.max = 5000), dist = "loglogistic" )
      lik2 <- as.numeric( logLik(fit2) )
      return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    parallel::stopCluster(cl)
    stat <- as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndTobit) ) {  ## Tobit regression
  fit1 = survival::survreg(target ~ 1, weights = wei, dist = "gaussian")
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = 2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
      fit2 = survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
      lik2 = as.numeric( logLik(fit2) )
      return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    parallel::stopCluster(cl)
    stat = as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndClogit) ) {  ## Conditional logistic regression
  case <- as.logical(target[, 1]);  ## case 
  subject <- target[, 2] #the patient id
  stat <- numeric(cols)
  dof <- numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <- survival::clogit( case ~ dataset[, i] + strata(subject) ) 
      dof[i] <- length( fit2$coefficients ) 
      stat[i] <- diff( fit2$loglik )
    }
    univariateModels$stat <- 2 * stat 
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::clogit(case ~ dataset[, i] + strata(subject) ) 
        return( c( diff(fit2$loglik) , length( fit2$coefficients ) ) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- 2 * as.vector( mod[, 1] )
    univariateModels$pvalue <- pchisq(mod[, 1], mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, censIndER) ) {  ## Exponential regression
  fit1 <- survival::survreg(target ~ 1, dist = "exponential", weights = wei)
  lik1 <- as.numeric( logLik(fit1) )
  lik2 <- numeric(cols)
  dof <- numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <- survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
      lik2[i] <- as.numeric( logLik(fit2) )
      dof[i] <- length( coef(fit2) ) - 1
    }
    stat <- 2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
        return( c(as.numeric( logLik(fit2) ), length( coef(fit2) ) - 1 ) )
    }
    parallel::stopCluster(cl)
    stat = as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, mod[ ,2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndQBinom)   ) {  ## Quasi Binomial regression

  fit1 <- glm(target ~ 1, family = quasibinomial(link = logit), weights = wei)
  lik1 <- fit1$deviance
  lik2 <- numeric(cols)
  dof <- numeric(cols)
  ina <- 1:cols
  phi <- numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in ina ) {
      fit2 <- glm( target ~ dataset[, i], family = quasibinomial(link = logit), weights = wei )
      lik2[i] <- fit2$deviance
      phi[i] <- summary(fit2)[[ 14 ]]
      dof[i] <- length( fit2$coefficients)
    }
    stat <- (lik1 - lik2) / (dof - 1) / phi
    univariateModels$stat <- stat
    univariateModels$pvalue <- pf(stat, dof - 1, rows - dof, lower.tail = FALSE, log.p = TRUE)
  
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
      fit2 <- glm( target ~ dataset[, i], family = quasibinomial(link = logit), weights = wei )
      return( c(fit2$deviance, length( fit2$coefficients ), summary(fit2)[[14]] ) )
    }
    parallel::stopCluster(cl)
    stat <- as.vector( lik1 - mod[, 1] ) / (mod[, 2] - 1) / mod[, 3] 
    univariateModels$stat <- stat
    univariateModels$pvalue <- pf(stat, mod[, 2] - 1, rows - mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndQPois)   ) {  ## Quasi Poisson regression
  
  fit1 <- glm(target ~ 1, family = quasipoisson(link = log), weights = wei)
  lik1 <- fit1$deviance
  lik2 <- numeric(cols)
  dof <- numeric(cols)
  ina <- 1:cols
  phi <- numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in ina ) {
      fit2 <- glm( target ~ dataset[, i], family = quasipoisson(link = log), weights = wei )
      lik2[i] <- fit2$deviance
      phi[i] <- summary(fit2)[[ 14 ]]
      dof[i] <- length( fit2$coefficients)
    }
    stat <- (lik1 - lik2) / (dof - 1) / phi
    univariateModels$stat <- stat
    univariateModels$pvalue <- pf(stat, dof - 1, rows - dof, lower.tail = FALSE, log.p = TRUE)
    
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
      fit2 <- glm( target ~ dataset[, i], family = quasipoisson(link = log), weights = wei )
      return( c(fit2$deviance, length( fit2$coefficients ), summary(fit2)[[14]] ) )
    }
    parallel::stopCluster(cl)
    stat <- as.vector( lik1 - mod[, 1] ) / (mod[, 2] - 1) / mod[, 3] 
    univariateModels$stat <- stat
    univariateModels$pvalue <- pf(stat, mod[, 2] - 1, rows - mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndSPML) ) {  ## Circular regression
  if ( !is.matrix(dataset) )  dataset <- as.matrix(dataset)
  mod <- Rfast::spml.regs(target, dataset, logged = TRUE, parallel = (ncores > 1) )
  univariateModels$stat <- mod[, 1]
  univariateModels$pvalue <- mod[, 2]
  
}  else   univariateModels <- NULL
  
  if ( !is.null(univariateModels) )  {
    if (targetID != - 1) {
      univariateModels$stat[targetID] <- 0
      univariateModels$pvalue[targetID] <- log(1)
    }
    if ( sum( id > 0 ) > 0 ) {
      univariateModels$stat[id] <- 0
      univariateModels$pvalue[id] <- log(1)
    }
  }
  
  univariateModels
}