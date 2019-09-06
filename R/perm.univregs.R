perm.univregs <- function(target, dataset, targetID = -1, test = NULL, user_test = NULL, wei = NULL, threshold = 0.05, R = 999, ncores = 1) {
  
  univariateModels <- list();
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]
  thres <- threshold * R + 1
  if (targetID != -1 ) {
    target <- dataset[, targetID]
    dataset[, targetID] <- rbinom(rows, 1, 0.5)
  }   
  id <- NULL
  if ( !identical(test, permFisher) ) {
    id <- Rfast::check_data(dataset)
    if ( sum(id > 0) )  dataset[, id] <- rbinom(rows * length(id), 1, 0.5 )
  }  
  la <- length( unique(target) )
  
if ( !is.null(user_test) ) {
  univariateModels <- perm.univariateScore(target, dataset, test, wei, targetID, threshold, R, ncores)

} else if ( identical(test, permMMFisher) )  { 
  univariateModels <- perm.univariateScore(target, dataset, test, wei, targetID, threshold, R, ncores)
  
} else if ( identical(test, permFisher) | ( identical(test, permReg) & is.null(wei) ) )  { ## Pearson's correlation 
  pvalue <- numeric(cols)
  a <- as.vector( cor(target, dataset) )
  dof <- rows - 3; #degrees of freedom
  wa <- abs( log( (1 + a) / (1 - a) ) )
  univariateModels$stat <- 0.5 * wa * sqrt(dof);
  dm <- dim(dataset)
  tb <- matrix(0, R, cols)
   for (j in 1:R) {
    yb <- sample(target, rows)
    tb[j, ] <- cor(yb, dataset)
  }  
  tb <- log( (1 + tb) / (1 - tb) )  ## the test statistic 
  for (i in 1:cols)  pvalue[i] <- ( sum( abs(tb[, i]) > wa[i] ) + 1 ) / (R + 1)  ## bootstrap p-value
  univariateModels$pvalue <- pvalue ;

} else if ( identical(test, permDcor) ) { ## Pearson's correlation 
  stat <- numeric(cols)
  pvalue <- numeric(cols)
  for (i in 1:cols) {
    mod <- energy::dcov.test(target, dataset[, i], R = R)
    stat[i] <- mod$statistic
    pvalue[i] <- mod$p.value
  }  
  univariateModels$stat <- stat ;
  univariateModels$pvalue <- pvalue ;

} else if ( identical(test, permgSquare) ) {  ## G^2 test
  z <- cbind(dataset, target)
  dc <- Rfast::colrange(z, cont = FALSE)
  a <- Rfast::g2tests_perm(data = z, x = 1:cols, y = cols + 1, dc = dc, nperm = R)
  univariateModels$stat <- a$statistic
  univariateModels$pvalue <- a$pvalue

} else if ( identical(test, permBeta) ) {  ## Beta regression
  mod <- perm.betaregs(target, dataset, wei, check = FALSE, logged = FALSE, R = R, threshold = threshold, ncores = ncores)
  univariateModels$stat <- mod[, 1]
  univariateModels$pvalue <- mod[, 2]

} else if ( identical(test, permMMReg)) {  ## M (Robust) linear regression
  fit1 <-MASS::rlm(target ~ 1, maxit = 2000, method = "MM")
  lik1 <-as.numeric( logLik(fit1) )
  lik2 <-numeric(cols)
  pval <-numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-MASS::rlm(target ~ dataset[, i], maxit = 2000, method = "MM" )
      lik2[i] = as.numeric( logLik(fit2) )
		  step <- 0
      j <- 1		
	   	x <- dataset[, i]
	   	while (j <= R & step < thres ) {
	   	  xb <- sample(x, rows)  
		    bit2 <- MASS::rlm(target ~ xb, maxit = 2000, method = "MM" )
        step <- step + ( as.numeric( logLik(bit2) ) > lik2[i] )
		    j <- j + 1
	   	}
		  pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- 2 * (lik2 - lik1)
    univariateModels$pvalue <- pval

  } else {
    
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "MASS") %dopar% {
      fit2 <-MASS::rlm(target ~ dataset[, i], maxit = 2000, method = "MM")
      lik2 <-as.numeric( logLik(fit2) )
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
		    bit2 <- MASS::rlm(target ~ xb, maxit = 2000, method = "MM")
        step <- step + ( as.numeric( logLik(bit2) ) > lik2[i] )
		    j <- j + 1
	    }
      return( c(lik2, step) ) 
    }  
    parallel::stopCluster(cl)
    univariateModels$stat <- 2 * abs(lik1 - mod[, 1])
    univariateModels$pvalue <- (mod[, 2] + 1) / (R + 1) 
  }   
  
} else if ( identical(test, permReg)  &  !is.null(wei) ) {  ## Weighted linear regression
  
  univariateModels = list();
  stat <-pval <-numeric(cols)
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
      tab <-anova(fit2)
     stat[i] <- tab[1, 4] 
		  step <- 0
      j <- 1		
		  x <- dataset[, i]
      while ( j <= R & step < thres ) {
		    xb <- sample(x, rows)  
		    bit2 <- lm(target ~ xb, weights = wei )
		    tab <- anova(bit2)[1, 4]
        step <- step + ( tab > stat[i] )
		    j <- j + 1
	   	}
      pval[i] <- (step + 1) / (R + 1) 
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
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while ( j <= R & step < thres ) {
		    xb <- sample(x, rows)  
		    bit2 <- lm(target ~ xb, weights = wei )
		    tab <- anova(bit2)[1, 4]
        step <- step + ( tab > stat )
		    j <- j + 1
		  }
      return( c(stat, step) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- (mod[, 2] + 1) / (R + 1)
  }   
  
} else if ( identical(test, permMVreg) ) {  ## Weighted linear regression
  
  univariateModels = list();
  stat <-pval <-numeric(cols)
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
     stat[i] <- anova(fit2)[2, 3]
      step <- 0
      j <- 1		
      x <- dataset[, i]
      while ( j <= R & step < thres ) {
        xb <- sample(x, rows)  
        bit2 <- lm(target ~ xb, weights = wei )
        step <- step + ( anova(bit2)[2, 3] > stat[i] )
        j <- j + 1
      }
      pvalue <- (step + 1) / (R + 1) 
    } 	  
    univariateModels$stat <- stat
    univariateModels$pvalue <- pval
    
  } else {
    
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
      fit2 <-lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
      stat <-anova(fit2)[2, 3]
      step <- 0
      j <- 1		
      x <- dataset[, i]
      while ( j <= R & step < thres ) {
        xb <- sample(x, rows)  
        bit2 <- lm(target ~ xb, weights = wei )
        step <- step + ( anova(bit2)[2, 3] > stat )
        j <- j + 1
      }
      return( c(stat, step) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- (mod[, 2] + 1) / (R + 1) 
  }   

} else if ( identical(test, permLogistic)  &  is.ordered(target) ) {  ## ordinal regression
  lik2 <- numeric(cols)
  pval <- numeric(cols)
  fit1 <- ordinal::clm(target ~ 1, weights = wei)
  lik1 <- as.numeric( logLik(fit1) )
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <- ordinal::clm.fit(target ~ dataset[, i], weights = wei)
      lik2[i] <- as.numeric( fit2$logLik )
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- ordinal::clm.fit(target ~ xb, weights = wei)            
		    step <- step + ( bit2$logLik > lik2[i] )
		    j <- j + 1
	    }
		  pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- 2 * (lik2 - lik1) 
    univariateModels$pvalue <- pval

  } else {
    
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "ordinal") %dopar% {
      fit2 <- ordinal::clm.fit(target ~ dataset[, i], weights = wei)
      lik2 <- as.numeric( fit2$logLik )
		  step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- ordinal::clm.fit(target ~ xb, weights = wei)  
        step <- step + ( bit2$logLik > lik2 )
		    j <- j + 1
	    }
		  pval <- (step + 1) / (R + 1) 
      return( c(lik2, pval ) )
    } 
    parallel::stopCluster(cl)
    univariateModels$stat <- as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$pvalue <- mod[, 2]
  }
  
} else if ( identical(test, permLogistic) == TRUE  &  la > 2  ) {  ## multinomial regression
  
  target = as.factor( as.numeric( target ) );
  lik2 <-numeric(cols)
  pval <-numeric(cols)
  fit1 <-nnet::multinom(target ~ 1, trace = FALSE, weights = wei)
  lik1 <-as.numeric( logLik(fit1) )

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-nnet::multinom(target ~ dataset[, i], trace = FALSE, weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
	    step <- 0
      j <- 1		
    	x <- dataset[, i]
      while (j <= R & step < thres ) {
	      xb <- sample(x, rows)  
        bit2 <- nnet::multinom(target ~ xb, trace = FALSE, weights = wei)  
        step <- step + ( as.numeric( logLik(bit2) ) > lik2[i] )
	      j <- j + 1
      }
	    pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- 2 * (lik2 - lik1)
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "nnet") %dopar% {
      fit2 <-nnet::multinom(target ~ dataset[, i], trace = FALSE, weights = wei)
      lik2 <-as.numeric( logLik(fit2 ) )
		  step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- nnet::multinom(target ~ xb, trace = FALSE, weights = wei)  
        step <- step + ( as.numeric( logLik(bit2 ) ) > lik2 )
		    j <- j + 1
	    }
		  pval <- (step + 1) / (R + 1) 
      return( c(lik2, pval ) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$pvalue <- pval
  }
  
} else if ( identical(test, permLogistic)  &  la == 2  &  is.matrix(dataset)  &  is.null(wei) ) {  ## logistic regression
  mod <- Rfast::univglms(target, dataset, oiko = "binomial", logged = FALSE)
  stat <- mod[, 1]
  bt <- matrix(0, R, cols)
  pval <- numeric(cols)
  for (j in 1:R) {
    y <- sample(target, rows)
    bt[j, ] <- Rfast::univglms(y, dataset, oiko = "binomial", logged = FALSE)[, 1]
  }
  for (i in 1:cols)  pval[i] <- (sum(bt[, i] > stat[i]) + 1) / (R + 1)
  univariateModels$stat <- stat
  univariateModels$pvalue <- pval

} else if ( identical(test, permLogistic)  &  la == 2  &  ( !is.null(wei)  ||  is.data.frame(dataset)  ) ) {  ## Logistic regression
  fit1 <-glm(target ~ 1, binomial, weights = wei)
  dev1 = fit1$deviance
  dev2 = numeric(cols)
  pval <-numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-glm( target ~ dataset[, i], binomial, weights = wei )
      dev2[i] = fit2$deviance
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, binomial, weights = wei)  
        step <- step + ( bit2$deviance < dev2[i] )
		    j <- j + 1
	    }
		  pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- dev1 - dev2
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
      fit2 <-glm( target ~ dataset[, i], binomial, weights = wei )
      dev2 = fit2$deviance
		  step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, binomial, weights = wei)  
        step <- step + ( bit2$deviance < dev2 )
		    j <- j + 1
	    }
		  pval <- (step + 1) / (R + 1) 
      return( c(dev2, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- as.vector( dev1 - mod[, 1])
    univariateModels$pvalue <- mod[, 2]
  }
  
} else if ( identical(test, permBinom) ) {  ## Binomial regression
  wei <- target[, 2] 
  y <- target[, 1] / wei
  fit1 <-glm(y ~ 1, binomial, weights = wei)
  dev1 = as.numeric( logLik(fit1) )
  dev2 = numeric(cols)
  pval <-numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-glm( y ~ dataset[, i], binomial, weights = wei )
      dev2[i] = fit2$deviance
		  step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, binomial, weights = wei)  
        step <- step + ( bit2$deviance < dev2[i] )
		    j <- j + 1
	    }
		  pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- dev1 - dev2
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    wei <- target[, 2] 
    y <- target[, 1] / wei
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
      fit2 <-glm( y ~ dataset[, i], binomial, weights = wei )
      dev2 = fit2$deviance
		  step <- 0
      j <- 1		
	    x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, binomial, weights = wei)  
        step <- step + ( bit2$deviance < dev2 )
		    j <- j + 1
	    }
		  pval <- (step + 1) / (R + 1) 
      return( c(dev2, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- as.vector( dev1 - mod[, 1])
    univariateModels$pvalue <- mod[, 2]
  }
  
} else if ( identical(test, permPois)  &  is.matrix(dataset)  &  is.null(wei) ) {  ## Poisson regression
  mod <- Rfast::univglms(target, dataset, oiko = "poisson", logged = FALSE)
  stat <- mod[, 1]
  bt <- matrix(0, R, cols)
  pval <- numeric(cols)
  for (j in 1:R) {
    y <- sample(target, rows)
    bt[j, ] <- Rfast::univglms(y, dataset, oiko = "poisson", logged = FALSE)[, 1]
  }
  for (i in 1:cols)  pval[i] <- (sum(bt[, i] > stat[i]) + 1) / (R + 1)
  univariateModels$stat <- stat
  univariateModels$pvalue <- pval

} else if ( identical(test, permPois)  &  ( !is.null(wei)  ||  is.data.frame(dataset)  ) ) {  ## Poisson regression
  fit1 <-glm(target ~ 1, poisson, weights = wei)
  dev1 = fit1$deviance
  dev2 = numeric(cols)
  pval <-numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-glm( target ~ dataset[, i], poisson, weights = wei )
      dev2[i] = fit2$deviance
		  step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, poisson, weights = wei)  
        step <- step + ( bit2$deviance < dev2[i] )
		    j <- j + 1
	    }
		  pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- dev1 - dev2
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
      fit2 <-glm( target ~ dataset[, i], poisson, weights = wei )
      dev2 = fit2$deviance
		  step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, poisson, weights = wei)  
        step <- step + ( bit2$deviance < dev2 )
		    j <- j + 1
	    }
		  pval <- (step + 1) / (R + 1) 
      return( c(dev2, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- as.vector(dev1 - mod[, 1])
    univariateModels$pvalue <- mod[, 2]
  }
  
} else if ( identical(test, permGamma) ) {  ## Gamma regression
  fit1 <-glm(target ~ 1, family = Gamma(link = log), weights = wei)
  dev1 = fit1$deviance
  dev2 = numeric(cols)
  pval <-numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-glm( target ~ dataset[, i], family = Gamma(link = log), weights = wei )
      dev2[i] = fit2$deviance
      step <- 0
      j <- 1		
      x <- dataset[, i]
      while (j <= R & step < thres ) {
        xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, family = Gamma(link = log), weights = wei)  
        step <- step + ( bit2$deviance < dev2[i] )
        j <- j + 1
      }
      pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- dev1 - dev2
    univariateModels$pvalue <- pval
    
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
      fit2 <-glm( target ~ dataset[, i], family = Gamma(link = log), weights = wei )
      dev2 = fit2$deviance
      step <- 0
      j <- 1		
      x <- dataset[, i]
      while (j <= R & step < thres ) {
        xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, family = Gamma(link = log), weights = wei)  
        step <- step + ( bit2$deviance < dev2 )
        j <- j + 1
      }
      pval <- (step + 1) / (R + 1) 
      return( c(dev2, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- as.vector(dev1 - mod[, 1])
    univariateModels$pvalue <- mod[, 2]
  }
  
} else if ( identical(test, permNormLog) ) {  ## Gaussian regression with a log link
  fit1 <-glm(target ~ 1, family = gaussian(link = log), weights = wei)
  dev1 = fit1$deviance
  dev2 = numeric(cols)
  pval <-numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
      dev2[i] = fit2$deviance
      step <- 0
      j <- 1		
      x <- dataset[, i]
      while (j <= R & step < thres ) {
        xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, family = gaussian(link = log), weights = wei)  
        step <- step + ( bit2$deviance < dev2[i] )
        j <- j + 1
      }
      pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- dev1 - dev2
    univariateModels$pvalue <- pval
    
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
      fit2 <-glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
      dev2 = fit2$deviance
      step <- 0
      j <- 1		
      x <- dataset[, i]
      while (j <= R & step < thres ) {
        xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, family = gaussian(link = log), weights = wei)  
        step <- step + ( bit2$deviance < dev2 )
        j <- j + 1
      }
      pval <- (step + 1) / (R + 1) 
      return( c(dev2, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- as.vector(dev1 - mod[, 1])
    univariateModels$pvalue <- mod[, 2]
  }
  
} else if ( identical(test, permNB) ) {  ## Negative binomial regression
  dev1 <- MASS::glm.nb( target ~ 1, weights = wei )$deviance

  if ( ncores <= 1 | is.null(ncores) ) {
    dev2 <- pval <- numeric(cols)
    for ( i in 1:cols ) {
      fit2 <-MASS::glm.nb( target ~ dataset[, i], weights = wei )
      dev2[i] = fit2$deviance
		  step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- MASS::glm.nb(target ~ xb, weights = wei)  
        step <- step + ( bit2$deviance < dev2[i] )
		    j <- j + 1
	    }
		  pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- dev1 - dev2
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "MASS") %dopar% {
      fit2 <-MASS::glm.nb( target ~ dataset[, i], weights = wei )
      dev2 = fit2$deviance
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, poisson, weights = wei)  
        step <- step + ( bit2$deviance < dev2 )
		    j <- j + 1
	    }
	    pval <- (step + 1) / (R + 1) 
      return( c(dev2, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- dev1 - mod[, 1]
    univariateModels$pvalue <- pval
  }
  
} else if ( identical(test, permZIP) ) {  ## Zero-inflated Poisson regression
  moda <- perm.zipregs(target, dataset, wei, check = FALSE, logged = FALSE, R = R, threshold = threshold, ncores = ncores) 
  univariateModels$stat <- moda[, 1]
  univariateModels$pvalue <- moda[, 2]

} else if ( identical(test, permRQ) ) {  ## Median (quantile) regression
  
  fit1 <-quantreg::rq(target ~ 1, weights = wei)
  stat <-pval <-numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-quantreg::rq(target ~ dataset[, i], weights = wei )
      ww = anova(fit1, fit2, test = "rank")
     stat[i] <- as.numeric( ww[[1]][3] )
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 = quantreg::rq(target ~ xb, weights = wei )
        ww = anova(fit1, bit2, test = "rank")
        step <- step + ( as.numeric( ww[[1]][3] ) > stat[i] )
		    j <- j + 1
	    }
		  pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- stat
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "quantreg") %dopar% {
      fit2 <-quantreg::rq(target ~ dataset[, i], weights = wei )
      ww = anova(fit1, fit2, test = "rank")
      stat <-as.numeric( ww[[1]][3] )
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 = quantreg::rq(target ~ xb, weights = wei )
        ww = anova(fit1, bit2, test = "rank")
        step <- step + ( as.numeric( ww[[1]][3] ) > stat )
		    j <- j + 1
	    }
		  pval <- (step + 1) / (R + 1) 
      return( c(stat, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- as.vector( mod[, 1] )
    univariateModels$pvalue <- as.vector( mod[, 2] )
  }
  
} else if ( identical(test, permIGreg) ) {  ## Inverse Gaussian regression
  fit1 <-glm(target ~ 1, family = inverse.gaussian(link = log), weights = wei)
  lik1 <-as.numeric( logLik(fit1) )
  lik2 <-numeric(cols)
  pval <-numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
		  step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, family = inverse.gaussian(link = log), weights = wei)  
        step <- step + ( as.numeric( logLik(bit2) ) > lik2[i] )
		    j <- j + 1
	    }
		  pval[i] <- (step + 1) / (R + 1) 
    }
    stat <-2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
      fit2 <-glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
      lik2 <-as.numeric( logLik(fit2) )
      step <- 0
      j <- 1		
	    x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- glm(target ~ xb, poisson, weights = wei)  
        step <- step + ( as.numeric( logLik(bit2) ) > lik2 )
		    j <- j + 1
	    }
		  pval <- (step + 1) / (R + 1) 
      return( c(lik2, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$pvalue <- mod[, 2]
  }
  
} else if ( identical(test, permCR) ) {  ## Cox regression
  stat <-numeric(cols)
  pval <-numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    fit2 <-survival::coxph( target ~ dataset[, i], weights = wei )
    stat[i] <- anova(fit2)[2, 2]
    step <- 0
    j <- 1		
	  x <- dataset[, i]
    while (j <= R & step < thres ) {
	    xb <- sample(x, rows)  
      bit2 = survival::coxph( target ~ xb, weights = wei )
      step <- step + ( anova(bit2)[2, 2] > stat[i] )
	    j <- j + 1
    }
	  pval[i] <- (step + 1) / (R + 1) 
    univariateModels$stat <- stat
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
      fit2 <-survival::coxph( target ~ dataset[, i], weights = wei )
      stat <- anova(fit2)[2, 2]
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 = survival::coxph( target ~ xb, weights = wei )
        step <- step + ( anova(bit2)[2, 2] > stat )
		    j <- j + 1
	    }
		  pval <- (step + 1) / (R + 1) 
      return( c(stat, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]
  }
  
} else if ( identical(test, permWR) ) {  ## Weibull regression
  fit1 <-survival::survreg(target ~ 1, weights = wei)
  lik1 <-as.numeric( logLik(fit1) )
  lik2 <-numeric(cols)
  pval <-numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-survival::survreg( target ~ dataset[, i], weights = wei, control = list(iter.max = 5000) )
      lik2[i] = as.numeric( logLik(fit2) )
		  step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- survival::survreg( target ~ xb, weights = wei, control = list(iter.max = 5000) ) 
        qa <- ( as.numeric( logLik(bit2) ) > lik2[i] )
        if ( is.na(qa) )  qa <- 0
        step <- step + qa
		    j <- j + 1
	    }
		  pval[i] <- (step + 1) / (R + 1) 
    }
    stat <-2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
      fit2 <-survival::survreg( target ~ dataset[, i], weights = wei, control = list(iter.max = 5000) )
		  lik2 <- as.numeric( logLik(fit2) )
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- survival::survreg( target ~ xb, weights = wei, control = list(iter.max = 5000) ) 
        qa <- ( as.numeric( logLik(bit2) ) > lik2[i] )
        if ( is.na(qa) )  qa <- 0
        step <- step + qa
		    j <- j + 1
	    }
		  pval <- (step + 1) / (R + 1) 
      return( c(lik2, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- 2 * (mod[, 1] - lik1) 
    univariateModels$pvalue <- mod[, 2]
  }

} else if ( identical(test, permLLR) ) {  ## Log-logistic regression
  fit1 <-survival::survreg(target ~ 1, weights = wei, dist = "loglogistic")
  lik1 <-as.numeric( logLik(fit1) )
  lik2 <-numeric(cols)
  pval <-numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-survival::survreg( target ~ dataset[, i], weights = wei, control = list(iter.max = 5000), dist = "loglogistic" )
      lik2[i] = as.numeric( logLik(fit2) )
      step <- 0
      j <- 1		
      x <- dataset[, i]
      while (j <= R & step < thres ) {
        xb <- sample(x, rows)  
        bit2 <- survival::survreg( target ~ xb, weights = wei, control = list(iter.max = 5000), dist = "loglogistic" ) 
        qa <- ( as.numeric( logLik(bit2) ) > lik2[i] )
        if ( is.na(qa) )  qa <- 0
        step <- step + qa
        j <- j + 1
      }
      pval[i] <- (step + 1) / (R + 1) 
    }
    stat <- 2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pval
    
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
      fit2 <- survival::survreg( target ~ dataset[, i], weights = wei, control = list(iter.max = 5000), dist = "loglogistic" )
      lik2 <- as.numeric( logLik(fit2) )
      step <- 0
      j <- 1		
      x <- dataset[, i]
      while (j <= R & step < thres ) {
        xb <- sample(x, rows)  
        bit2 <- survival::survreg( target ~ xb, weights = wei, control = list(iter.max = 5000), dist = "loglogistic" ) 
        qa <- ( as.numeric( logLik(bit2) ) > lik2[i] )
        if ( is.na(qa) )  qa <- 0
        step <- step + qa
        j <- j + 1
      }
      pval <- (step + 1) / (R + 1) 
      return( c(lik2, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- 2 * (mod[, 1] - lik1) 
    univariateModels$pvalue <- mod[, 2]
  }
  
  
} else if ( identical(test, permTobit) ) {  ## Tobit regression
  fit1 <-survival::survreg(target ~ 1, weights = wei, dist = "gaussian")
  lik1 <-as.numeric( logLik(fit1) )
  lik2 <-numeric(cols)
  pval <-numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
      lik2[i] = as.numeric( logLik(fit2) )
      step <- 0
      j <- 1		
      x <- dataset[, i]
      while (j <= R & step < thres ) {
        xb <- sample(x, rows)  
        bit2 <- survival::survreg( target ~ xb, weights = wei, dist = "gaussian" ) 
        qa <- ( as.numeric( logLik(bit2) ) > lik2[i] )
        if ( is.na(qa) )  qa <- 0
        step <- step + qa
        j <- j + 1
      }
      pval[i] <- (step + 1) / (R + 1) 
    }
    stat <-2 * (lik2 - lik1)
    univariateModels$stat <- stat
    univariateModels$pvalue <- pval
    
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
      fit2 <-survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
      lik2 <- as.numeric( logLik(fit2) )
      step <- 0
      j <- 1		
      x <- dataset[, i]
      while (j <= R & step < thres ) {
        xb <- sample(x, rows)  
        bit2 <- survival::survreg( target ~ xb, weights = wei, dist = "gaussian" ) 
        qa <- ( as.numeric( logLik(bit2) ) > lik2[i] )
        if ( is.na(qa) )  qa <- 0
        step <- step + qa
        j <- j + 1
      }
      pval <- (step + 1) / (R + 1) 
      return( c(lik2, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- 2 * (mod[, 1] - lik1) 
    univariateModels$pvalue <- mod[, 2]
  }  
  
} else if ( identical(test, permClogit) ) {  ## Conditional logistic regression
  subject = target[, 2] #the patient id
  case = as.logical(target[, 1]);  ## case 
  stat <-numeric(cols)
  pval <-numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-survival::clogit( case ~ dataset[, i] + strata(subject) ) 
     stat[i] <- 2 * abs( diff(fit2$loglik) )
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- survival::clogit( case ~ xb + strata(subject) ) 
        step <- step + ( 2 * abs( diff(bit2$loglik) ) > stat[i] )
		    j <- j + 1
	    }
		  pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- stat
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
      fit2 <-survival::clogit(case ~ dataset[, i] + strata(subject) ) 
		  stat <-2 * abs( diff(fit2$loglik) )
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- survival::clogit( case ~ xb + strata(subject) ) 
        step <- step + ( 2 * abs( diff(bit2$loglik) ) > stat )
		    j <- j + 1
	    }
		  pval <- (step + 1) / (R + 1) 
      return( c(stat, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[ ,2]
  }
  
} else if ( identical(test, permER) ) {  ## Exponential regression
  fit1 <-survival::survreg(target ~ 1, dist = "exponential", weights = wei)
  lik1 <-as.numeric( logLik(fit1) )
  lik2 <-numeric(cols)
  pval <-numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 <-survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
		  step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- survival::survreg( target ~ xb, dist = "exponential", weights = wei ) 
        qa <- ( as.numeric( logLik(bit2) ) > lik2[i] )
        if ( is.na(qa) )  qa <- 0
        step <- step + qa 
		    j <- j + 1
	    }
		  pval[i] <- (step + 1) / (R + 1) 
    }
    univariateModels$stat <- 2 * (lik2 - lik1)
    univariateModels$pvalue <- pval

  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
      fit2 <-survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
      lik2 <-as.numeric( logLik(fit2) )
      step <- 0
      j <- 1		
		  x <- dataset[, i]
      while (j <= R & step < thres ) {
		    xb <- sample(x, rows)  
        bit2 <- survival::survreg( target ~ xb, dist = "exponential", weights = wei ) 
        qa <- ( as.numeric( logLik(bit2) ) > lik2[i] )
        if ( is.na(qa) )  qa <- 0
        step <- step + qa
		    j <- j + 1
	    }
		  pval <- (step + 1) / (R + 1) 
      return( c(lik2, pval) )
    }
    parallel::stopCluster(cl)
    univariateModels$stat <- 2 * (mod[, 1] - lik1)
    univariateModels$pvalue <- mod[, 2]
  }
  
} else   univariateModels <- NULL
   
  if ( !is.null(univariateModels) )  {
    univariateModels$pvalue <- log( univariateModels$pvalue )
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
