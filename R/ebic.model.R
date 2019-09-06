ebic.model <- function(target, test = NULL, wei = NULL) {
  
  n <- length(target)
  if ( identical(test, censIndCR) | identical(test, censIndWR) | identical(test, censIndER) | identical(test, testIndTobit) ) {
    n <- 0.5 * n 
  }
  logn <- log(n)

  if ( identical(test, testIndBeta) ) {  ## Beta regression
    fit2 <- beta.mod(target, dataset = NULL, wei = wei)
    ebic <-  -2 * fit2$loglik + 2 * logn
    
  } else if ( identical(test, testIndMMReg) ) {  ## M (Robust) linear regression
    fit2 <- MASS::rlm(target ~ 1, maxit = 2000, method = "MM" )
    ebic <- BIC(fit2) 
    
  } else if ( identical(test, testIndReg) ) {  ## linear regression
    fit2 <- lm( target ~ 1, weights = wei, y = FALSE, model = FALSE )
    ebic <- BIC(fit2)

  } else if ( identical(test, testIndOrdinal) ) {  ## ordinal regression
    fit2 <- ordinal::clm(target ~ 1, weights = wei)
    ebic <- BIC(fit2)

  } else if ( identical(test, testIndMultinom) ) {  ## multinomial regression
    fit2 <- nnet::multinom(target ~ 1, trace = FALSE, weights = wei )
    ebic <- BIC(fit2)

  } else if ( identical(test, testIndLogistic) ) {  ## logistic regression
    fit2 <- glm(y ~ 1, binomial, weights = wei)
    ebic <- BIC(fit2)

  } else if ( identical(test, testIndBinom) ) {  ## Binomial regression
    wei <- target[, 2] 
    y <- target[, 1] / wei
    fit2 <- glm( y ~ 1, binomial, weights = wei )
    ebic <- BIC(fit2)       

  } else if ( identical(test, testIndPois) ) {  ## Poisson regression
    fit2 <- glm( target ~ 1, poisson, weights = wei )
    ebic <- BIC(fit2)

  } else if ( identical(test, testIndNB) ) {  ## Negative binomial regression
    fit2 <- MASS::glm.nb( target ~ 1, weights = wei )
    ebic <- BIC(fit2)

  } else if ( identical(test, testIndNormLog)   ) {  ## Normal log link regression
    fit2 <- glm( target ~ 1, family = gaussian(link = log), weights = wei )
    ebic <- BIC(fit2)
    
  } else if ( identical(test, testIndGamma)   ) {  ## Gamma regression
    fit2 <- glm( target ~ 1, family = Gamma(link = log), weights = wei )
    ebic <- BIC(fit2)

  } else if ( identical(test, testIndZIP) ) {  ## Zero-inflated Poisson regression
    fit2 <- zip.mod(target, dataset = NULL, wei = wei) 
    ebic <-  -2 * fit2$loglik + 2 * logn
    
  } else if ( identical(test, testIndIGreg) ) {  ## Inverse Gaussian regression
    fit2 <- glm( target ~ 1, family = inverse.gaussian(link = log), weights = wei )
    ebic <- BIC(fit2)

  } else if ( identical(test, censIndCR) ) {  ## Cox regression
    fit2 <- survival::coxph( target ~ 1, weights = wei)
    ebic <-  - 2 * fit2$loglik

  } else if ( identical(test, censIndWR) ) {  ## Weibull regression
    fit2 <- survival::survreg( target ~ 1, weights = wei )
    ebic <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn 
    
  } else if ( identical(test, censIndER) ) {  ## Exponential regression
    fit2 <- survival::survreg( target ~ 1, weights = wei, dist = "exponential" )
    ebic <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn 
    
  } else if ( identical(test, censIndLLR) ) {  ## Log-logistic regression
    fit2 <- survival::survreg( target ~ 1, weights = wei, dist = "loglogistic" )
    ebic <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn 

  } else if ( identical(test, testIndTobit) ) {  ## Tobit regression
    fit2 <- survival::survreg( target ~ 1, weights = wei, dist = "gaussian" )
    ebic <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn

  } else if ( identical(test, testIndClogit) ) {  ## Conditional logistic regression
    subject <- target[, 2] #the patient id
    case <- as.logical(target[, 1]);  ## case 
    fit2 <- survival::clogit( case ~ 1 + strata(subject) ) 
    ebic <- BIC(fit2)
    
  } else if ( identical(test, testIndSPML) ) {
    if ( is.matrix(target) )  target <- ( atan(target[, 2]/target[, 1]) + pi * I(target[, 1] < 0) ) %% (2 * pi)
    ebic <-  -2 *Rfast::spml.mle(target)$loglik + 2 * logn
    
  }  else  ebic <- NULL
  
  ebic
  
}