gee.ci.mm <- function(ind1, ind2, cs = NULL, dat, group, se = "jack") {
  
  y <- dat[, ind1]
  x <- dat[, ind2]
  
  if ( is.null(cs) ) {
    if ( Rfast::sort_unique.length(y) == 2 ) {
      mod1 <- try( geepack::geeglm( y ~ x, family = binomial(logit), id = group, corstr = "exchangeable", std.err = se), silent = TRUE)
      if ( identical( class(mod1), "try-error" ) ) {
        t1 <- 0
        p1 <- log(1)
      } else {
        t1 <- anova(mod1)[1, 2]
        p1 <- pchisq(t1, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    } else if ( sum(round(y) - y) == 0 ) {
      mod1 <- try( geepack::geeglm( y ~ x, family = poisson(log), id = group, corstr = "exchangeable", std.err = se ), silent = TRUE)
      if ( identical( class(mod1), "try-error" ) ) {
        t1 <- 0
        p1 <- log(1)
      } else {
        t1 <- anova(mod1)[1, 2]
        p1 <- pchisq(t1, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    } else {
      mod1 <- try( geepack::geeglm(y ~ x, family = gaussian, id = group, corstr = "exchangeable", std.err = se ), silent = TRUE)
      if ( identical( class(mod1), "try-error" ) ) {
        t1 <- 0
        p1 <- log(1)
      } else {
        t1 <- anova(mod1)[1, 2]
        p1 <- pchisq(t1, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    }  ## end if ( Rfast::sort_unique.length(y) == 2 ) 
    
    if ( Rfast::sort_unique.length(x) == 2 ) {
      mod2 <- try( geepack::geeglm( x ~ y, family = binomial(logit), id = group, corstr = "exchangeable", std.err = se ), silent = TRUE)
      if ( identical( class(mod2), "try-error" ) ) {
        t2 <- 0
        p2 <- log(1)
      } else {
        t2 <- anova(mod2)[1, 2]
        p2 <- pchisq(t2, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    } else if ( sum(round(x) - x) == 0 ) {
      mod2 <- try( geepack::geeglm( x ~ y, family = poisson(log), id = group, corstr = "exchangeable", std.err = se ), silent = TRUE)
      if ( identical( class(mod2), "try-error" ) ) {
        t2 <- 0
        p2 <- log(1)
      } else {
        t2 <- anova(mod2)[1, 2]
        p2 <- pchisq(t2, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    } else {
      mod2 <- try( geepack::geeglm(x ~ y, family = gaussian, id = group, corstr = "exchangeable", std.err = se ), silent = TRUE)
      if ( identical( class(mod2), "try-error" ) ) {
        t2 <- 0
        p2 <- log(1)
      } else {
        t2 <- anova(mod2)[1, 2]
        p2 <- pchisq(t2, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    }  ## end if ( Rfast::sort_unique.length(x) == 2 ) 

    
  } else {  ### with conditioning set   
    z <- dat[, cs]
    
    if ( Rfast::sort_unique.length(y) == 2 ) {
      mod1 <- try( geepack::geeglm( y ~ z + x, family = binomial(logit), id = group, corstr = "exchangeable", std.err = se ), silent = TRUE)
      if ( identical( class(mod1), "try-error" ) ) {
        t1 <- 0
        p1 <- log(1)
      } else {
        t1 <- anova(mod1)[2, 2]
        p1 <- pchisq(t1, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    } else if ( sum(round(y) - y) == 0 ) {
      mod1 <- try( geepack::geeglm( y ~ z + x, family = poisson(log), id = group, corstr = "exchangeable", std.err = se ), silent = TRUE)
      if ( identical( class(mod1), "try-error" ) ) {
        t1 <- 0
        p1 <- log(1)
      } else {
        t1 <- anova(mod1)[2, 2]
        p1 <- pchisq(t1, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    } else {
      mod1 <- try( geepack::geeglm( y ~ z + x, family = gaussian, id = group, corstr = "exchangeable", std.err = se ), silent = TRUE)
      if ( identical( class(mod1), "try-error" ) ) {
        t1 <- 0
        p1 <- log(1)
      } else {
        t1 <- anova(mod1)[2, 2]
        p1 <- pchisq(t1, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    }  ## end if ( Rfast::sort_unique.length(y) == 2 )
    
    if ( Rfast::sort_unique.length(x) == 2 ) {
      mod2 <- try( geepack::geeglm( x ~ z + y, family = binomial(logit), id = group, corstr = "exchangeable", std.err = se ), silent = TRUE)
      if ( identical( class(mod2), "try-error" ) ) {
        t2 <- 0
        p2 <- log(1)
      } else {
        t2 <- anova(mod2)[2, 2]
        p2 <- pchisq(t2, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    } else if ( sum(round(x) - x) == 0 ) {
      mod2 <- try( geepack::geeglm( x ~ z + y, family = poisson(log), id = group, corstr = "exchangeable", std.err = se ), silent = TRUE)
      if ( identical( class(mod2), "try-error" ) ) {
        t2 <- 0
        p2 <- log(1)
      } else {
        t2 <- anova(mod2)[2, 2]
        p2 <- pchisq(t2, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    } else {  
      mod2 <- try( geepack::geeglm( x ~ z + y, family = gaussian, id = group, corstr = "exchangeable", std.err = se ), silent = TRUE)
      if ( identical( class(mod2), "try-error" ) ) {
        t2 <- 0
        p2 <- log(1)
      } else {
        t2 <- anova(mod2)[2, 2]
        p2 <- pchisq(t2, 1, lower.tail = FALSE, log.p = TRUE)
      }	
    }  # end if ( Rfast::sort_unique.length(x) == 2 )
	
  }   ## end if ( is.null(cs) ) 
  
  pval <- min( log(2) + min(p1, p2), max(p1, p2) )
  stat <-  max( 2 * max(t1, t2), min(t1, t2) )
  result <- c(stat, pval, 1)
  names(result) <- c('test', 'logged.p-value', 'df') 
  result
}