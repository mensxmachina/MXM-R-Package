glmm.univregs <- function(target, reps = NULL, id, dataset, targetID = -1, test, wei = NULL, 
                              slopes = FALSE, ncores = 1) {
  
  if ( identical(test, testIndGLMMReg) | identical(test, testIndLMM) ) {
      univariateModels <- univariateScore.glmm(target = target, reps = reps, group = id, dataset = dataset, test = test, wei = wei, 
                                          targetID = targetID, slopes = slopes, ncores = ncores)  
  } else {
      
    cols <- dim(dataset)[2]
    ind <- 1:cols
    univariateModels <- list();
    stat <- numeric(cols)
    poia <- Rfast::check_data(dataset)
    
    if ( sum(poia) > 0 )  {
      stat[poia] <- 0
      ind[poia] <- 0
    }
  
    ### Logistic GLMM    
    if ( identical(test, testIndGLMMLogistic) ) {
      
      if ( is.null(reps) ) {
        fit1 <- lme4::glmer( target ~ (1|id), weights = wei, family = binomial ) 
        for (i in ind) {
          fit2 <- lme4::glmer( target ~ (1|id) + dataset[, i], weights = wei, family = binomial ) 
          mod <- anova(fit1, fit2)
          stat[i] <- mod[2, 6]
        }

      } else {
        reps <- reps 
        if ( slopes ) {
          fit1 <- lme4::glmer( target ~ reps + (reps|id), weights = wei, family = binomial ) 
          for (i in ind) {
            fit2 <- lme4::glmer( target ~ reps + (reps|id) + dataset[, i], weights = wei, family = binomial ) 
            mod <- anova(fit1, fit2)
            stat[i] <- mod[2, 6]
          }  ##  end for (i in ind)  
        } else {  ###  yes slopes
          reps <- reps 
          fit1 <- lme4::glmer( target ~ reps + (1|id), weights = wei, family = binomial )  
          for (i in ind) {
            fit2 <- lme4::glmer( target ~ reps + (1|id) + dataset[, i], weights = wei, family = binomial ) 
            mod <- anova(fit1, fit2)
            stat[i] <- mod[2, 6]
          }  ##  end for (i in ind)
        }  ##  end if  slopes
      }  ##  end if ( is.null(reps) ) 

      ### Poisson GLMM    
    } else if ( identical(test, testIndGLMMPois) ) {
    
      if ( is.null(reps) ) {
        fit1 <- lme4::glmer( target ~ (1|id), weights = wei, family = poisson ) 
        for (i in ind) {
          fit2 <- lme4::glmer( target ~ (1|id) + dataset[, i], weights = wei, family = poisson ) 
          mod <- anova(fit1, fit2)
          stat[i] <- mod[2, 6]
        }

      } else {
        reps <- reps 
        if ( slopes ) {
          fit1 <- lme4::glmer( target ~ reps + (reps|id), weights = wei, family = poisson ) 
          for (i in ind) {
            fit2 <- lme4::glmer( target ~ reps + (reps|id) + dataset[, i], weights = wei, family = poisson ) 
            mod <- anova(fit1, fit2)
            stat[i] <- mod[2, 6]
          }  ##  end for (i in ind)  
        } else {  ###  yes slopes
          reps <- reps 
          fit1 <- lme4::glmer( target ~ reps + (1|id), weights = wei, family = poisson )  
          for (i in ind) {
            fit2 <- lme4::glmer( target ~ reps + (1|id) + dataset[, i], weights = wei, family = poisson ) 
            mod <- anova(fit1, fit2)
            stat[i] <- mod[2, 6]
          }  ##  end for (i in ind)
        }  ##  end if  slopes
      }  ##  end if ( is.null(reps) ) 

      ### Gamma GLMM
    } else if ( identical(test, testIndGLMMGamma) ) {
    
      if ( is.null(reps) ) {
        fit1 <- lme4::glmer( target ~ (1|id), weights = wei, family = Gamma(log) ) 
        for (i in ind) {
          fit2 <- lme4::glmer( target ~ (1|id) + dataset[, i], weights = wei, family = Gamma(log) ) 
          mod <- anova(fit1, fit2)
          stat[i] <- mod[2, 6]
        }
      
      } else {
        reps <- reps 
        if ( slopes ) {
          fit1 <- lme4::glmer( target ~ reps + (reps|id), weights = wei, family = Gamma(log) ) 
          for (i in ind) {
            fit2 <- lme4::glmer( target ~ reps + (reps|id) + dataset[, i], weights = wei, family = Gamma(log) ) 
            mod <- anova(fit1, fit2)
            stat[i] <- mod[2, 6]
          }  ##  end for (i in ind)  
        } else {  ###  yes slopes
          reps <- reps 
          fit1 <- lme4::glmer( target ~ reps + (1|id), weights = wei, family = Gamma(log) )  
          for (i in ind) {
            fit2 <- lme4::glmer( target ~ reps + (1|id) + dataset[, i], weights = wei, family = Gamma(log) ) 
            mod <- anova(fit1, fit2)
            stat[i] <- mod[2, 6]
          }  ##  end for (i in ind)
        }  ##  end if  slopes
      }  ##  end if ( is.null(reps) ) 

      ### Gaussian with log link GLMM
    } else if ( identical(test, testIndGLMMNormLog) ) {
      
      if ( is.null(reps) ) {
        fit1 <- lme4::glmer( target ~ (1|id), weights = wei, family = gaussian(log) ) 
        for (i in ind) {
          fit2 <- lme4::glmer( target ~ (1|id) + dataset[, i], weights = wei, family = gaussian(log) ) 
          mod <- anova(fit1, fit2)
          stat[i] <- mod[2, 6]
        }
        
      } else {
        reps <- reps 
        if ( slopes ) {
          fit1 <- lme4::glmer( target ~ reps + (reps|id), weights = wei, family = gaussian(log) ) 
          for (i in ind) {
            fit2 <- lme4::glmer( target ~ reps + (reps|id) + dataset[, i], weights = wei, family = gaussian(log) ) 
            mod <- anova(fit1, fit2)
            stat[i] <- mod[2, 6]
          }  ##  end for (i in ind)  
        } else {  ###  yes slopes
          reps <- reps 
          fit1 <- lme4::glmer( target ~ reps + (1|id), weights = wei, family = gaussian(log) )  
          for (i in ind) {
            fit2 <- lme4::glmer( target ~ reps + (1|id) + dataset[, i], weights = wei, family = gaussian(log) ) 
            mod <- anova(fit1, fit2)
            stat[i] <- mod[2, 6]
          }  ##  end for (i in ind)
        }  ##  end if  slopes
      }  ##  end if ( is.null(reps) )         

      ### Ordinal GLMM
    } else if ( identical(test, testIndGLMMOrdinal) ) {
      
      if ( is.null(reps) ) {
        fit1 <- ordinal::clmm( target ~ (1|id), weights = wei ) 
        for (i in ind) {
          fit2 <- ordinal::clmm( target ~ (1|id) + dataset[, i], weights = wei ) 
          mod <- anova(fit1, fit2)
          stat[i] <- mod[2, 6]
        }
        
      } else {
        reps <- reps 
        if ( slopes ) {
          fit1 <- ordinal::clmm( target ~ reps + (reps|id), weights = wei ) 
          for (i in ind) {
            fit2 <- ordinal::clmm( target ~ reps + (reps|id) + dataset[, i], weights = wei ) 
            mod <- anova(fit1, fit2)
            stat[i] <- mod[2, 6]
          }  ##  end for (i in ind)  
        } else {  ###  yes slopes
          reps <- reps 
          fit1 <- ordinal::clmm( target ~ reps + (1|id), weights = wei )  
          for (i in ind) {
            fit2 <- ordinal::clmm( target ~ reps + (1|id) + dataset[, i], weights = wei ) 
            mod <- anova(fit1, fit2)
            stat[i] <- mod[2, 4]
          }  ##  end for (i in ind)
        }  ##  end if  slopes
      }  ##  end if ( is.null(reps) )         
     
      ### Mixed effects Cox regression
    } else if ( identical(test, testIndGLMMCR) ) {
      
      fit1 <- coxme::coxme( target ~ (1|id), weights = wei ) 
      for (i in ind) {
        fit2 <- coxme::coxme( target ~ (1|id) + dataset[, i], weights = wei ) 
        mod <- anova(fit1, fit2)
        stat[i] <- mod[2, 2]
      }  

    }
    
    univariateModels$stat <- stat
    univariateModels$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)  
    
  }  ##  end  if ( identical(test, testIndGLMMReg) | identical(test, testIndLMM) ) 
  
  univariateModels
}
