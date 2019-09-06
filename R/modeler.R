modeler <- function(target,dataset = NULL, test = "testIndFisher") {
  
  n <- length(target)
  if ( survival::is.Surv(target) )  n <- 0.5 * n
  
  if ( is.null(dataset) ) {
    
    ## Linear regression
    if ( test == "testIndFisher" | test == "testIndReg" ) {
      mod <- lm(target ~ 1)
      res <- resid(mod)
      dev <- deviance(mod)
      bic <- dev
      
      ## MM regression
    } else if ( test == "testIndMMReg" ) {
      mod <- MASS::rlm(target ~ 1, maxit = 2000, method = "MM")
      res <- resid(mod)
      dev <-  - 2 * as.numeric( logLik(mod) )
      bic <- BIC(mod)  
      
      ## Quantile regression
    } else if ( test == "testIndRQ" ) {
      mod <- quantreg::rq(target ~ 1)
      res <- resid(mod)
      dev <-  - 2 * as.numeric( logLik(mod) )
      bic <- dev + 3 * log(n) 
      
      ## Gamma regression
    } else if ( test == "testIndGamma" ) {
      mod <- glm(target ~ 1, family = Gamma(log) ) 
      res <- resid(mod, type = "response")
      dev <-  deviance(mod)
      bic <- BIC(mod)
      
      ## Gaussian regression with log link
    } else if ( test == "testIndNormLog" ) {
      mod <- glm(target ~ 1, family = gaussian(log) ) 
      res <- resid(mod, type = "response")
      dev <-  deviance(mod)
      bic <- BIC(mod)
      
      ## binary logistic regression
    } else if ( test == "testIndLogistic" ) {
      mod <- glm(target ~ 1, binomial)
      res <- resid(mod, type = "response")
      dev <-  deviance(mod)
      bic <- BIC(mod)
      
      ## multinomial regression
    } else if ( test == "testIndMultinom" ) {
      mod <- nnet::multinom(target ~ 1, trace = FALSE)
      res <- resid(mod)
      dev <-  mod$deviance
      bic <- BIC(mod)
      
      ## ordinal regression
    } else if ( test == "testIndOrdinal" ) {
      mod <- MASS::polr(target ~ 1)
      res <- ord.resid(target, mod$fitted.values)
      dev <- deviance(mod)
      bic <- BIC(mod)
      
      ## Poisson regression
    } else if ( test == "testIndPois" ) {
      mod <- glm(target ~ 1, poisson)
      res <- resid(mod, type = "response")
      dev <-  deviance(mod)
      bic <- BIC(mod)
      
      ## quasi binomial regression
    } else if ( test == "testIndQBinom" ) {
      mod <- glm(target ~ 1, quasibinomial)
      res <- resid(mod, type = "response")
      dev <-  - 2 * as.numeric( logLik(mod) )
      bic <- NA
      
      ## quasi Poisson regression
    } else if ( test == "testIndQPois" ) {
      mod <- glm(target ~ 1, quasipoisson)
      res <- resid(mod, type = "response")
      dev <-  - 2 * as.numeric( logLik(mod) )
      bic <- NA
      
      ## negative binomial regression
    } else if ( test == "testIndNB" ) {
      mod <- MASS::glm.nb(target ~ 1 )
      res <- resid(mod, type = "response")
      dev <-  - 2 * as.numeric( logLik(mod) )
      bic <- BIC(mod)
      
      ## beta regression
    } else if ( test == "testIndBeta" ) {
      mod <- Rfast::beta.mle(target)
      res <- target - mod$param[1]/sum(mod$param)
      dev <-  - 2 * mod$loglik
      bic <- dev + 2 * log(n)
      
      ## Tobit regression
    } else if ( test == "testIndTobit" ) {
      mod <- survival::survreg(target ~ 1, dist = "gaussian" )
      res <- resid(mod)
      dev <-  - 2 * mod$loglik[2]
      bic <- dev + 2 * log(n)
      
      ## Cox proportional hazards
    } else if ( test == "censIndCR" ) {
      mod <- survival::coxph(target ~ 1)
      res <- mod$residuals
      dev <-  -2 * mod$loglik 
      bic <- dev 
      
      ## Weibull regression
    } else if ( test == "censIndWR" ) {
      mod <- survival::survreg(target ~ 1 )
      res <- resid(mod)
      dev <-  -2 * mod$loglik[2]
      bic <- dev + 2 * log(n)
      
      ## Log-logistic regression
    } else if ( test == "censIndLLR" ) {
      mod <- survival::survreg(target ~ 1, dist = "loglogistic" )
      res <- resid(mod)
      dev <-  -2 * mod$loglik[2]
      bic <- dev + 2 * log(n)
    }
    
    ##########################
    
  } else {  ## ifdataset is not NULL
    
    if ( test == "testIndFisher" | test == "testIndReg" ) {
      mod <- lm(target ~ dataset)
      res <- resid(mod)
      dev <- deviance(mod)
      bic <- dev
      
      ## MM regression
    } else if ( test == "testIndMMReg" ) {
      mod <- try( MASS::rlm(target ~ dataset, maxit = 2000, method = "MM"), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        res <- NA
        dev <- NA
        bic <- NA 
      } else { 
        res <- resid(mod)
        dev <-  - 2 * as.numeric( logLik(mod) )
        bic <- BIC(mod)
      }
      
      ## Quantile regression
    } else if ( test == "testIndRQ" ) {
      mod <- quantreg::rq(target ~ dataset)
      res <- resid(mod)
      dev <-  - 2 * as.numeric( logLik(mod) )
      bic <- dev + (length(mod$coefficients) + 1) * log(n)
      
      ## Gamma regression
    } else if ( test == "testIndGamma" ) {
      mod <- glm(target ~ dataset, family = Gamma(log) ) 
      res <- resid(mod, type = "response")
      dev <-  deviance(mod)
      bic <- BIC(mod)
      
      ## Gaussian regression with log link
    } else if ( test == "testIndNormLog" ) {
      mod <- glm(target ~ dataset, family = gaussian(log) ) 
      res <- resid(mod, type = "response")
      dev <-  deviance(mod)
      bic <- BIC(mod)
      
      ## binary logistic regression
    } else if ( test == "testIndLogistic" ) {
      mod <- glm(target ~ dataset, binomial)
      res <- resid(mod, type = "response")
      dev <-  deviance(mod)
      bic <- BIC(mod)

      ## multinomial regression
    } else if ( test == "testIndMultinom" ) {
      mod <- nnet::multinom(target ~ dataset, trace = FALSE)
      res <- resid(mod)
      dev <-  mod$deviance
      bic <- BIC(mod)
      
      ## ordinal regression
    } else if ( test == "testIndOrdinal" ) {
      mod <- MASS::polr(target ~ dataset)
      res <- ord.resid(target, mod$fitted.values)
      dev <-  deviance(mod)
      bic <- BIC(mod)
      
      ## Poisson regression
    } else if ( test == "testIndPois" ) {
      mod <- glm(target ~ dataset, poisson)
      res <- resid(mod, type = "response")
      dev <-  deviance(mod)
      bic <- BIC(mod)
      
      ## quasi binomial regression
    } else if ( test == "testIndQBinom" ) {
      mod <- glm(target ~ dataset, quasibinomial)
      res <- resid(mod, type = "response")
      dev <-  - 2 * as.numeric( logLik(mod) )
      bic <- NA
      
      ## quasi Poisson regression
    } else if ( test == "testIndQPois" ) {
      mod <- glm(target ~ dataset, quasipoisson)
      res <- resid(mod, type = "response")
      dev <-  - 2 * as.numeric( logLik(mod) )
      bic <- NA
      
      ## negative binomial regression
    } else if ( test == "testIndNB" ) {
      mod <- MASS::glm.nb(target ~ dataset )
      res <- resid(mod, type = "response")
      dev <-  - 2 * as.numeric( logLik(mod) )
      bic <- BIC(mod)
      
      ## beta regression
    } else if ( test == "testIndBeta" ) {
      mod <- try( beta.mod(target, dataset, xnew = dataset), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        res <- NA
        dev <- NA
        bic <- NA
      } else {
        res <- target - mod$est
        dev <-  - 2 * mod$loglik
        bic <- dev + (length(mod$be[, 1]) + 1) * log(n)
      }
      
      ## Tobit regression
    } else if ( test == "testIndTobit" ) {
      mod <- try( survival::survreg(target ~ dataset, dist = "gaussian" ), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        res <- NA
        dev <- NA
        bic <- NA
      } else {
        res <- resid(mod)
        dev <-  -2 * mod$loglik[2]
        bic <- dev + (length(mod$coefficients) + 1) * log(n)
      } 
      
      ## Cox proportional hazards
    } else if ( test == "censIndCR" ) {
      mod <- try( survival::coxph(target ~ dataset), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        res <- NA
        dev <- NA
      } else {
        res <- mod$residuals
        dev <-  -2 * mod$loglik[2]
        bic <- BIC(mod)
      }
      
      ## Weibull regression
    } else if ( test == "censIndWR" ) {
      mod <- try( survival::survreg(target ~ dataset ), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        res <- NA
        dev <- NA
        bic <- NA
      } else {
        res <- resid(mod)
        dev <-  -2 * mod$loglik[2]
        bic <- dev + (length(mod$coefficients) + 1) * log(n)
      }
      
      ## Log-logistic regression
    } else if ( test == "censIndLLR" ) {
      mod <- try( survival::survreg(target ~ dataset, dist = "loglogistic" ), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        res <- NA
        dev <- NA
        bic <- NA
      } else {
        res <- resid(mod)
        dev <-  -2 * mod$loglik[2]
        bic <- dev + (length(mod$coefficients) + 1) * log(n)
      }
      
    }
    
  } ## end if ( is.null(dataset) )
  
  list(mod = mod, dev = dev, bic = bic, res = res)  
}