ebic.regs <- function(target, dataset, xIndex, csIndex, gam = NULL, test = NULL, wei = NULL, ncores = 1) {
  
  if ( identical(csIndex, 0) ) {
    models <- MXM::ebic.univregs(target = target, dataset = dataset, test = test, wei = wei, ncores = ncores, gam = gam)
  } else {

    if ( length(xIndex) > 0 ) {
      oop <- options(warn = -1) 
      on.exit( options(oop) )
      dm <- dim(dataset)
      n <- dm[1]
      p <- dm[2]
      lik2 <- numeric(p)
      logn <- log(n)
      
      if ( is.null(gam) ) {
        con <- 2 - log(p) / logn
      } else con <- 2 * gam
      if ( (con) <= 0 )  con <- 0
      M <- length(csIndex) + 1
      common <- con * lchoose(p, M)
      
      if ( identical(test, testIndBeta) ) {  ## Beta regression
        for (i in xIndex) {
          fit2 <- beta.reg(target, dataset[, c(csIndex, xIndex[i] ) ])$loglik
          lik2[i] <-  -2 * fit2$loglik + (length(fit2$be) + 1) * logn 
        }

      } else if ( identical(test, testIndMMReg) ) {  ## M (Robust) linear regression

        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- MASS::rlm( target ~., data = dataset[, c(csIndex, i)], maxit = 2000, method = "MM" )
            lik2[i] <- BIC(fit2) 
          } 

        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "MASS") %dopar% {
            fit2 <- MASS::rlm( target ~., data = dataset[, c(csIndex, i)], maxit = 2000, method = "MM" )
            return( BIC(fit2) )
          }  
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )
        }   
        
      } else if ( identical(test, testIndReg) ) {  ## linear regression
        
        if ( ncores <= 1 | is.null(ncores) ) {

          for ( i in xIndex) {
            fit2 <- lm( target ~., data = dataset[, c(csIndex, i)], weights = wei, y = FALSE, model = FALSE )
            lik2[i] <- BIC(fit2)
          }
          
        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
            fit2 <- lm( target ~., data = dataset[, c(csIndex, i)], weights = wei, y = FALSE, model = FALSE )
            return( BIC(fit2) )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector(mod[, 1])
        }   
        
      } else if ( identical(test, testIndOrdinal) ) {  ## ordinal regression

        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- ordinal::clm( y ~., data = dataset[, c(csIndex, i)], weights = wei, model = FALSE )
            lik2[i] <-  BIC(fit2)
          }
 
        } else {
          
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "ordinal") %dopar% {
            fit2 <- ordinal::clm( y ~., data = dataset[, c(csIndex, i)], weights = wei, model = FALSE )
            return( BIC(fit2) )
          } 
          parallel::stopCluster(cl)
          lik2 <- as.vector(mod[, 1])
          
        }
        
      } else if ( identical(test, testIndMultinom) ) {  ## multinomial regression
        
        target <- as.factor( as.numeric( target ) );
        
        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- nnet::multinom( target ~., data = dataset[, c(csIndex, i)], weights = wei, trace = FALSE )
            lik2[i] <- BIC(fit2)
          }

        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "nnet") %dopar% {
            fit2 <- nnet::multinom( target ~., data = dataset[, c(csIndex, i)], weights = wei, trace = FALSE )
            return( BIC(fit2) )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector(mod[, 1])
        }
        
      } else if ( identical(test, testIndLogistic) ) {  ## Logistic regression
        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], binomial, weights = wei )
            lik2[i] <- BIC(fit2)
          }

        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
            fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], binomial, weights = wei )
            return( BIC(fit2) )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector(mod[, 1])
        }
        
      } else if ( identical(test, testIndBinom) ) {  ## Binomial regression
        wei <- target[, 2] 
        y <- target[, 1] / wei
        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- glm( y ~., data = dataset[, c(csIndex, i)], binomial, weights = wei )
            lik2[i] <- BIC(fit2)
          }

        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          wei <- target[, 2] 
          y <- target[, 1] / wei
          mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
            fit2 <- glm( y ~., data = dataset[, c(csIndex, i)], binomial, weights = wei )
            return( BIC(fit2) )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector(mod[, 1])
        }
        
      } else if ( identical(test, testIndPois) ) {  ## Poisson regression

        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], poisson, weights = wei )
            lik2[i] <- BIC(fit2)
          }

        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
            fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], poisson, weights = wei )
            return( BIC(fit2) )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )
        }
        
      } else if ( identical(test, testIndNB) ) {  ## Negative binomial regression

        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- MASS::glm.nb( target ~., data = dataset[, c(csIndex, i)], weights = wei )
            lik2[i] <- BIC(fit2)
          }

        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "MASS") %dopar% {
            fit2 <- MASS::glm.nb( target ~., data = dataset[, c(csIndex, i)], weights = wei )
            return( BIC(fit2) )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector(mod[, 1])
        }

      } else if ( identical(test, testIndNormLog) ) {  ## Normal log link regression

        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = gaussian(link = log), weights = wei )
            lik2[i] <- BIC(fit2)
          }
 
        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
            fit2 <- glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
            return( BIC(fit2) )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )

        }
        
      } else if ( identical(test, testIndGamma)   ) {  ## Gamma regression
        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = Gamma(link = log), weights = wei )
            lik2[i] <- BIC(fit2)
          }

        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
            fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = Gamma(link = log), weights = wei )
            return( BIC(fit2) )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )
        } 
        
      } else if ( identical(test, testIndIGreg) ) {  ## Inverse Gaussian regression

        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = inverse.gaussian(link = log), weights = wei )
            lik2[i] <- BIC(fit2)
          }

        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind) %dopar% {
            fit2 <- glm( target ~., data = dataset[, c(csIndex, i)], family = inverse.gaussian(link = log), weights = wei )
            return( BIC(fit2) )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )
        }
        
      } else if ( identical(test, testIndZIP) ) {  ## Zero-inflated Poisson regression
        for (i in xIndex) {
          fit2 <- zip.reg(target, dataset[, c(csIndex, xIndex[i] ) ])$loglik
          lik2[i] <-  -2 * fit2$loglik + (length(fit2$be) + 1) * logn
        }
        
      } else if ( identical(test, censIndCR) ) {  ## Cox regression
        
        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- survival::coxph( target ~., data = dataset[, c(csIndex, i)], weights = wei)
            lik2[i] <- BIC(fit2)
          }
          
        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
            fit2 <- survival::coxph( target ~., data = dataset[, c(csIndex, i)], weights = wei )
            return( BIC(fit2) )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )
        }
        
      } else if ( identical(test, censIndWR) ) {  ## Weibull regression
        
        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei )
            lik2[i] <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn
          }

        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
            fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei )
            return( - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )
        }
        
      } else if ( identical(test, censIndER) ) {  ## Exponential regression
        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], dist = "exponential", weights = wei )
            lik2[i] <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn
          }

        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
            fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], dist = "exponential", weights = wei )
            return( - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )
        }
        
      } else if ( identical(test, censIndLLR) ) {  ## Log-logistic regression
        
        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei, dist = "loglogistic" )
            lik2[i] <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn
          }
          
        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
            fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei, dist = "loglogistic" )
            return( - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )
        }
        
      } else if ( identical(test, testIndTobit) ) {  ## Tobit regression
        
        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei, dist = "gaussian" )
            dof[i] <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn
          }

          
        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
            fit2 <- survival::survreg( target ~., data = dataset[, c(csIndex, i)], weights = wei, dist = "gaussian" )
            return( - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )
        }
        
      } else if ( identical(test, testIndClogit) ) {  ## Conditional logistic regression
        subject <- target[, 2] #the patient id
        case <- as.logical(target[, 1]);  ## case 
        
        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- survival::clogit( case ~., data = dataset[, c(csIndex, i)] + strata(subject) ) 
            lik2[i] <- BIC(fit2)
          }

        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "survival") %dopar% {
            fit2 <- survival::clogit( case ~., data = dataset[, c(csIndex, i)] + strata(subject) ) 
            return( BIC(fit2) )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )
        }
        
      } else if ( identical(test, testIndSPML) ) {  ## Conditional logistic regression
        if ( ncores <= 1 | is.null(ncores) ) {
          for ( i in xIndex ) {
            fit2 <- Rfast::spml.reg( target, x = dataset[, c(csIndex, i)] ) 
            lik2[i] <-  -2 * fit2$loglik + length(fit2$be) * logn
          }
          
        } else {
          cl <- parallel::makePSOCKcluster(ncores)
          doParallel::registerDoParallel(cl)
          mod <- foreach::foreach(i = xIndex, .combine = rbind, .packages = "Rfast") %dopar% {
            fit2 <- Rfast::spml.reg( target, data = dataset[, c(csIndex, i)] ) 
            return( -2 * fit2$loglik + length(fit2$be) * logn )
          }
          parallel::stopCluster(cl)
          lik2 <- as.vector( mod[, 1] )
        }
        
      }  else  lik2 <- NULL  ## end of all if (test == )
      
      lik2 <- lik2 + common
      lik2[ - xIndex ] <- 0
      
    } else {
      lik2 <- NULL
    }  ##  end  if  ( length(xIndex) > 0 )  
    
  }  ## end if ( identical(csIndex, 0) )
  
  lik2
}