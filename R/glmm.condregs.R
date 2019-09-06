glmm.condregs <- function(target, reps = NULL, id, dataset, xIndex, csIndex, test, wei = NULL, 
                                slopes = FALSE, ncores = 1) {
    
  if ( identical(csIndex, 0) ) {
    models <- glmm.univregs(target = target, reps = reps, id = id, dataset = dataset, targetID = -1, test = test, 
                                wei = wei, slopes = slopes, ncores = ncores)
  } else {
    
    if ( length(xIndex) > 0 ) {
      
    models <- list();
    cols <- dim(dataset)[2]
    lik2 <- numeric(cols)
    
    oop <- options(warn = -1) 
    on.exit( options(oop) )
	
    if ( identical(test, testIndLMM) | ( identical(test, testIndGLMMReg)  &  is.null(wei)  &  !slopes  &  is.null(reps) )  ) {
      stat <- numeric(cols)
      
        for (i in xIndex) {
          fit2 <- Rfast::rint.reg( target, dataset[, c(csIndex, i)], id )
          be <- fit2$be
          seb <- fit2$seb
          n <- length(target)
          p <- length(be)
          if  ( length( unique(round(fit2$be, 14) ) ) < p ) {  ## overloaded, nnon significant variable, probably collinear
            stat[i] <- 0
          } else   stat[i] <- be[p]^2/seb[p]^2    
        }  ##  end  for (i in xIndex) {

      models$stat <- stat
      models$pvalue <- pf(stat, 1, n - p - 2, lower.tail = FALSE, log.p = TRUE)

    ### Gaussian GLMM    
    } else if ( identical(test, testIndGLMMReg) ) {
      stat <- numeric(cols)
      cs <- dataset[, csIndex]
      dcs <- length(csIndex) + 1
      
      if ( is.null(reps) ) {
        for (i in xIndex) {
          fit2 <- lme4::lmer( target ~ (1|id) + cs + dataset[, i], weights = wei, REML = FALSE )
          if ( dcs < summary(fit2)[[ 3 ]]$dims[3] ) {
            mod <- anova(fit2)
            v2 <- as.numeric( summary(fit2)[[ 14 ]][5] )
            pr <- nrow(mod) 
            v1 <- mod[pr, 1]
            stat[i] <- mod[pr, 4]   
          }  
        }  ##  end  for (i in xIndex) {
        
      } else {
        reps <- reps 
        if ( slopes ) {
          for (i in xIndex) {
            fit2 <- lme4::lmer( target ~ reps + (reps|id) + cs + dataset[, i], weights = wei, REML = FALSE ) 	  
            if ( dcs < summary(fit2)[[ 3 ]]$dims[3] ) {
              mod <- anova(fit2)
              v2 <- as.numeric( summary(fit2)[[ 14 ]][5] )
              pr <- nrow(mod) 
              v1 <- mod[pr, 1]
              stat[i] <- mod[pr, 4]   
            } 
          }  ##  end for (i in xIndex)  
        } else {  ###  no slopes
          reps <- reps 
          for (i in xIndex) {
            fit2 <- lme4::lmer( target ~ reps + (1|id) + cs + dataset[, i], weights = wei, REML = FALSE )
            if ( dcs < summary(fit2)[[ 3 ]]$dims[3] ) {
              mod <- anova(fit2)
              v2 <- as.numeric( summary(fit2)[[ 14 ]][5] )
              pr <- nrow(mod) 
              v1 <- mod[pr, 1]
              stat[i] <- mod[pr, 4] 
            }      
          }  ##  end for (i in xIndex)
        }  ##  end if  slopes
      }  ##  end if ( is.null(reps) ) 
      
      models$stat <- stat
      models$pvalue <- pf(stat, v1, v2, lower.tail = FALSE, log.p = TRUE)
      
    ### Logistic GLMM    
    } else if ( identical(test, testIndGLMMLogistic) ) {
      cs <- dataset[, csIndex]
      
      if ( is.null(reps) ) {
        lik1 <- logLik( lme4::glmer( target ~ (1|id) + cs, weights = wei, family = binomial ) )
        for (i in xIndex) {
          fit2 <- lme4::glmer( target ~ (1|id) + cs + dataset[, i], weights = wei, family = binomial ) 
          lik2[i] <- logLik(fit2)
        }

      } else {
        reps <- reps 
        if ( slopes ) {
          lik1 <- logLik( lme4::glmer( target ~ reps + (reps|id) + cs, weights = wei, family = binomial ) )
          for (i in xIndex) {
            fit2 <- lme4::glmer( target ~ reps + (reps|id) + cs + dataset[, i], weights = wei, family = binomial ) 
            lik2[i] <- logLik(fit2)
          }  ##  end for (i in xIndex)  
        } else {  ###  yes slopes
          reps <- reps 
          lik1 <- logLik( lme4::glmer( target ~ reps + (1|id) + cs, weights = wei, family = binomial ) )
          for (i in xIndex) {
            fit2 <- lme4::glmer( target ~ reps + (1|id) + cs + dataset[, i], weights = wei, family = binomial ) 
            lik2[i] <- logLik(fit2)
          }  ##  end for (i in xIndex)
        }  ##  end if  slopes
      }  ##  end if ( is.null(reps) ) 
      
      models$stat <- 2 * lik2 - 2 * lik1
      models$pvalue <- pchisq(models$stat, 1, lower.tail = FALSE, log.p = TRUE)
    
      ### Poisson GLMM    
    } else if ( identical(test, testIndGLMMPois) ) {
      cs <- dataset[, csIndex]
      
      if ( is.null(reps) ) {
        lik1 <- logLik( lme4::glmer( target ~ (1|id) + cs, weights = wei, family = poisson ) )
        for (i in xIndex) {
          fit2 <- lme4::glmer( target ~ (1|id) + cs + dataset[, i], weights = wei, family = poisson ) 
          lik2[i] <- logLik(fit2)
        }
        
      } else {
        reps <- reps 
        if ( slopes ) {
          lik1 <- logLik( lme4::glmer( target ~ reps + (reps|id) + cs, weights = wei, family = poisson ) )
          for (i in xIndex) {
            fit2 <- lme4::glmer( target ~ reps + (reps|id) + cs + dataset[, i], weights = wei, family = poisson ) 
            lik2[i] <- logLik(fit2)
          }  ##  end for (i in xIndex)  
        } else {  ###  yes slopes
          reps <- reps 
          lik1 <- logLik( lme4::glmer( target ~ reps + (1|id) + cs, weights = wei, family = poisson ) )
          for (i in xIndex) {
            fit2 <- lme4::glmer( target ~ reps + (1|id) + cs + dataset[, i], weights = wei, family = poisson) 
            lik2[i] <- logLik(fit2)
          }  ##  end for (i in xIndex)
        }  ##  end if  slopes
      }  ##  end if ( is.null(reps) ) 
      
      models$stat <- 2 * lik2 - 2 * lik1
      models$pvalue <- pchisq(models$stat, 1, lower.tail = FALSE, log.p = TRUE)

      ### Gamma GLMM
    } else if ( identical(test, testIndGLMMGamma) ) {
      cs <- dataset[, csIndex]
      
      if ( is.null(reps) ) {
        lik1 <- logLik( lme4::glmer( target ~ (1|id) + cs, weights = wei, family = Gamma(log) ) )
        for (i in xIndex) {
          fit2 <- lme4::glmer( target ~ (1|id) + cs + dataset[, i], weights = wei, family = Gamma(log) ) 
          lik2[i] <- logLik(fit2)
        }
        
      } else {
        reps <- reps 
        if ( slopes ) {
          lik1 <- logLik( lme4::glmer( target ~ reps + (reps|id) + cs, weights = wei, family = Gamma(log) ) )
          for (i in xIndex) {
            fit2 <- lme4::glmer( target ~ reps + (reps|id) + cs + dataset[, i], weights = wei, family = Gamma(log) ) 
            lik2[i] <- logLik(fit2)
          }  ##  end for (i in xIndex)  
        } else {  ###  yes slopes
          reps <- reps 
          lik1 <- logLik( lme4::glmer( target ~ reps + (1|id) + cs, weights = wei, family = Gamma(log) ) )
          for (i in xIndex) {
            fit2 <- lme4::glmer( target ~ reps + (1|id) + cs + dataset[, i], weights = wei, family = Gamma(log) ) 
            lik2[i] <- logLik(fit2)
          }  ##  end for (i in xIndex)
        }  ##  end if  slopes
      }  ##  end if ( is.null(reps) ) 
      
      models$stat <- 2 * lik2 - 2 * lik1
      models$pvalue <- pchisq(models$stat, 1, lower.tail = FALSE, log.p = TRUE)
      
      ### Gaussian with log link GLMM
    } else if ( identical(test, testIndGLMMNormLog) ) {
      cs <- dataset[, csIndex]
      
      if ( is.null(reps) ) {
        lik1 <- logLik( lme4::glmer( target ~ (1|id) + cs, weights = wei, family = gaussian(log) ) )
        for (i in xIndex) {
          fit2 <- lme4::glmer( target ~ (1|id) + cs + dataset[, i], weights = wei, family = gaussian(log) ) 
          lik2[i] <- logLik(fit2)
        }
        
      } else {
        reps <- reps 
        if ( slopes ) {
          lik1 <- logLik( lme4::glmer( target ~ reps + (reps|id) + cs, weights = wei, family = gaussian(log) ) )
          for (i in xIndex) {
            fit2 <- lme4::glmer( target ~ reps + (reps|id) + cs + dataset[, i], weights = wei, family = gaussian(log) ) 
            lik2[i] <- logLik(fit2)
          }  ##  end for (i in xIndex)  
        } else {  ###  yes slopes
          reps <- reps 
          lik1 <- logLik( lme4::glmer( target ~ reps + (1|id) + cs, weights = wei, family = gaussian(log) ) )
          for (i in xIndex) {
            fit2 <- lme4::glmer( target ~ reps + (1|id) + cs + dataset[, i], weights = wei, family = gaussian(log) ) 
            lik2[i] <- logLik(fit2)
          }  ##  end for (i in xIndex)
        }  ##  end if  slopes
      }  ##  end if ( is.null(reps) ) 
      
      models$stat <- 2 * lik2 - 2 * lik1
      models$pvalue <- pchisq(models$stat, 1, lower.tail = FALSE, log.p = TRUE)
      
      ### Ordinal GLMM
    } else if ( identical(test, testIndGLMMOrdinal) ) {
      cs <- dataset[, csIndex]
      
      if ( is.null(reps) ) {
        lik1 <- logLik( ordinal::clmm( target ~ (1|id) + cs, weights = wei ) )
        for (i in xIndex) {
          fit2 <- ordinal::clmm( target ~ (1|id) + cs + dataset[, i], weights = wei ) 
          lik2[i] <- logLik(fit2)
        }
        
      } else {
        reps <- reps 
        if ( slopes ) {
          lik1 <- logLik( ordinal::clmm( target ~ reps + (reps|id) + cs, weights = wei ) )
          for (i in xIndex) {
            fit2 <- ordinal::clmm( target ~ reps + (reps|id) + cs + dataset[, i], weights = wei ) 
            lik2[i] <- logLik(fit2)
          }  ##  end for (i in xIndex)  
        } else {  ###  yes slopes
          reps <- reps 
          lik1 <- logLik( ordinal::clmm( target ~ reps + (1|id) + cs, weights = wei ) )
          for (i in xIndex) {
            fit2 <- ordinal::clmm( target ~ reps + (1|id) + cs + dataset[, i], weights = wei ) 
            lik2[i] <- logLik(fit2)
          }  ##  end for (i in xIndex)
        }  ##  end if  slopes
      }  ##  end if ( is.null(reps) ) 
      
      ### Mixed effects Cox regression
    } else if ( identical(test, testIndGLMMCR) ) {
      cs <- dataset[, csIndex]
      
      lik1 <- coxme::coxme( target ~ (1|id) + cs, weights = wei )$loglik[2]
      for (i in xIndex)  lik2[i] <- coxme::coxme( target ~ (1|id) + cs + dataset[, i], weights = wei )$loglik[2] 
      
      models$stat <- 2 * lik2 - 2 * lik1
      models$pvalue <- pchisq(models$stat, 1, lower.tail = FALSE, log.p = TRUE)
      
    }  else   models <- NULL  ## end of all if (test == )
    
  }  else {
     models <- list()
     models$stat <- NULL
     models$pvalue <- NULL
	 
  }  ##  end  if  ( length(xIndex) > 0 )  
  
  }  ## end if ( identical(csIndex, 0) )
  
  models
}