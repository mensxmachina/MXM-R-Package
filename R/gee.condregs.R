gee.condregs <- function(target, reps = NULL, id, dataset, xIndex, csIndex, test, wei = NULL, 
                          correl = "echangeable", se = "jack", ncores = 1) {
  
  if ( identical(csIndex, 0) ) {
    models <- gee.univregs(target = target, reps = reps, id = id, dataset = dataset, targetID = -1, test = test, 
                                wei = wei, correl = correl, se = se, ncores = ncores)
  } else {
    
    if ( length(xIndex) > 0 ) {
      
    models <- list();
    cols <- dim(dataset)[2]
    stat <- numeric(cols)
    
    oop <- options(warn = -1) 
    on.exit( options(oop) )
    cs <- dataset[, csIndex]
    
    ### GEE normal regression  
    if ( identical(test, testIndGEEReg ) ) {
      
      for (i in xIndex) {
        if ( is.null(reps) ) {
          fit2 <- try( geepack::geeglm( target ~ cs + dataset[, i], family = gaussian, id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
        } else {
          fit2 <- try( geepack::geeglm( target ~ reps + cs + dataset[, i], family = gaussian, id = id, waves = reps, weights = wei, corstr = correl, std.err = se ), silent = TRUE) 
        }
        #calculate the p value and stat.
        if ( identical( class(fit2), "try-error" ) ) {
          stat[i] <- 0
        } else {
          mod <- summary(fit2)[[6]]
          nr <- dim(mod)[1]
          stat[i] <- mod[nr, 3]
        }	
      }  ##  end  for (i in xIndex) {
      
      models$stat <- stat
      models$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)

      ### GEE logistic regression  
      } else if ( identical(test, testIndGEELogistic ) ) {
        
        for (i in xIndex) {
          if ( is.null(reps) ) {
            fit2 <- try( geepack::geeglm( target ~ cs + dataset[, i], family = binomial(logit), id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
          } else {
            fit2 <- try( geepack::geeglm( target ~ reps + cs + dataset[, i], family = binomial(logit), id = id, waves = reps, weights = wei, corstr = correl, std.err = se ), silent = TRUE) 
          }
          #calculate the p value and stat.
          if ( identical( class(fit2), "try-error" ) ) {
            stat[i] <- 0
          } else  stat[i] <- anova(fit2)[2, 3]
        }  ##  end  for (i in xIndex) {
        
        models$stat <- stat
        models$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
   
        ### GEE Poisson regression  
      } else if ( identical(test, testIndGEEPois ) ) {
        
        for (i in xIndex) {
          if ( is.null(reps) ) {
            fit2 <- try( geepack::geeglm( target ~ cs + dataset[, i], family = poisson(log), id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
          } else {
            fit2 <- try( geepack::geeglm( target ~ reps + cs + dataset[, i], family = poisson(log), id = id, waves = reps, weights = wei, corstr = correl, std.err = se ), silent = TRUE) 
          }
          #calculate the p value and stat.
          if ( identical( class(fit2), "try-error" ) ) {
            stat[i] <- 0
          } else   stat[i] <- anova(fit2)[2, 3]
        }  ##  end  for (i in xIndex) {
        
        models$stat <- stat
        models$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      
        ### GEE Gamma regression  
      } else if ( identical(test, testIndGEEGamma ) ) {
        
        for (i in xIndex) {
          if ( is.null(reps) ) {
            fit2 <- try( geepack::geeglm( target ~ cs + dataset[, i], family = Gamma(log), id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
          } else {
            fit2 <- try( geepack::geeglm( target ~ reps + cs + dataset[, i], family = Gamma(log), id = id, waves = reps, weights = wei, corstr = correl, std.err = se ), silent = TRUE) 
          }
          #calculate the p value and stat.
          if ( identical( class(fit2), "try-error" ) ) {
            stat[i] <- 0
          } else   stat[i] <- anova(fit2)[2, 3]
        }  ##  end  for (i in xIndex) {
        
        models$stat <- stat
        models$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        
        ### GEE Gaussian regression with log link
      } else if ( identical(test, testIndGEENormLog ) ) {
        
        for (i in xIndex) {
          if ( is.null(reps) ) {
            fit2 <- try( geepack::geeglm( target ~ cs + dataset[, i], family = gaussian(log), id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
          } else {
            fit2 <- try( geepack::geeglm( target ~ reps + cs + dataset[, i], family = gaussian(log), id = id, waves = reps, weights = wei, corstr = correl, std.err = se ), silent = TRUE) 
          }
          #calculate the p value and stat.
          if ( identical( class(fit2), "try-error" ) ) {
            stat[i] <- 0
          } else   stat[i] <- anova(fit2)[2, 3]
        }  ##  end  for (i in xIndex) {
        
        models$stat <- stat
        models$pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        
      } else   models <- NULL  ## end of all if (test == )
    
    }  else {
      models <- list()
      models$stat <- NULL
      models$pvalue <- NULL
    }  ##  end  if  ( length(xIndex) > 0 )  
    
  }  ## end if ( identical(csIndex, 0) )
  
  models
}