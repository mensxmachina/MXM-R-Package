glmm.ci.mm <- function(ind1, ind2, cs = NULL, dat, group) {
  
  y <- dat[, ind1]
  x <- dat[, ind2]
  n <- length(y)
  
  if ( is.null(cs) ) {
    dof <- 1
    if ( Rfast::sort_unique.length(y) == 2 ) {
      mod1 <- lme4::glmer(y ~ x + (1|group), family = "binomial")
      mod0 <- lme4::glmer(y ~ 1 + (1|group), family = "binomial")
      t1 <- anova(mod0, mod1)[2, 6]
      p1 <- pchisq(t1, 1, lower.tail = FALSE, log.p = TRUE)
    } else if ( sum(round(y) - y) == 0 ) {
      mod1 <- lme4::glmer(y ~ x + (1|group), family = "poisson")
      mod0 <- lme4::glmer(y ~ 1 + (1|group), family = "poisson")
      t1 <- anova(mod0, mod1)[2, 6]
      p1 <- pchisq(t1, 1, lower.tail = FALSE, log.p = TRUE)
    } else {
      dof <- n - 4
      mod1 <- Rfast::rint.reg(y, x, group)
      if ( is.infinite(mod1$info[5])  | length( unique(round(mod1$be, 14) ) ) < length(mod1$be) ) {
        p1 <- log(1)
        t1 <- 0
      } else {
        t1 <- ( mod1$be[2, ] / mod1$seb[2, ] )^2
        p1 <- pf(t1, 1, dof, lower.tail = FALSE, log.p = TRUE)
      }  ## end if ( is.infinite(mod1$info[5]) ) 
    }  ## end if ( Rfast::sort_unique.length(y) == 2 ) 
    
    if ( Rfast::sort_unique.length(x) == 2 ) {
      mod2 <- lme4::glmer(x ~ y + (1|group), family = "binomial")
      mod0 <- lme4::glmer(x ~ 1 + (1|group), family = "binomial")
      t2 <- anova(mod0, mod2)[2, 6]
      p2 <- pchisq(t2, 1, lower.tail = FALSE, log.p = TRUE)
    } else if ( sum(round(x) - x) == 0 ) {
      mod2 <- lme4::glmer(x ~ y + (1|group), family = "poisson")
      mod0 <- lme4::glmer(x ~ 1 + (1|group), family = "poisson")
      t2 <- anova(mod0, mod2)[2, 6]
      p2 <- pchisq(t2, 1, lower.tail = FALSE, log.p = TRUE)
    } else {
      dof <- n - 4
      mod2 <- Rfast::rint.reg(x, y, group)
	  if ( is.infinite(mod2$info[5])  | length( unique(round(mod2$be, 14) ) ) < length(mod2$be) ) {
        p2 <- log(1)
        t2 <- 0
      } else {
        t2 <- ( mod2$be[2, ] / mod2$seb[2, ] )^2
        p2 <- pf(t2, 1, dof, lower.tail = FALSE, log.p = TRUE)
      }  ## end if ( is.infinite(mod2$info[5]) ) 
    }  ## end if ( Rfast::sort_unique.length(x) == 2 ) 

    
  } else {  ### with conditioning set   
    z <- dat[, cs]
    dof <- 1
    if ( Rfast::sort_unique.length(y) == 2 ) {
      mod1 <- lme4::glmer(y ~ z + x + (1|group), family = "binomial")
      mod0 <- lme4::glmer(y ~ z + (1|group), family = "binomial")
      t1 <- anova(mod0, mod1)[2, 6]
      p1 <- pchisq(t1, 1, lower.tail = FALSE, log.p = TRUE)
    } else if ( sum(round(y) - y) == 0 ) {
      mod1 <- lme4::glmer(y ~ z + x + (1|group), family = "poisson")
      mod0 <- lme4::glmer(y ~ z + (1|group), family = "poisson")
      t1 <- anova(mod0, mod1)[2, 6]
      p1 <- pchisq(t1, 1, lower.tail = FALSE, log.p = TRUE)
    } else {
      mod0 <- Rfast::rint.reg(y, x, group)
      mod1 <- Rfast::rint.reg(y, cbind(x, z), group)
      if ( mod1$info[5] >= mod0$info[5] ) {
        p1 <- log(1)
        t1 <- 0
      } else {
        dof <- n - length(mod1$be) - 2
        t1 <- ( mod1$be[2, ] / mod1$seb[2, ] )^2
        p1 <- pf(t1, 1, dof, lower.tail = FALSE, log.p = TRUE)
      }  ## end if ( mod1$info[5] >= mod0$info[5] ) 
    }  ## end if ( Rfast::sort_unique.length(y) == 2 )
    
    if ( Rfast::sort_unique.length(x) == 2 ) {
      mod2 <- lme4::glmer(x ~ z + y + (1|group), family = "binomial")
      mod0 <- lme4::glmer(x ~ z + (1|group), family = "binomial")
      t2 <- anova(mod0, mod2)[2, 6]
      p2 <- pchisq(t2, 1, lower.tail = FALSE, log.p = TRUE)
    } else if ( sum(round(x) - x) == 0 ) {
      mod2 <- lme4::glmer(x ~ z + y + (1|group), family = "poisson")
      mod0 <- lme4::glmer(x ~ z + (1|group), family = "poisson")
      t2 <- anova(mod0, mod2)[2, 6]
      p2 <- pchisq(t2, 1, lower.tail = FALSE, log.p = TRUE)
    } else {  
      mod0 <- Rfast::rint.reg(x, y, group)
      mod2 <- Rfast::rint.reg(x, cbind(y, z), group)
      if ( mod2$info[5] >= mod0$info[5] ) {
        p2 <- log(1)
        t2 <- 0
      } else {
        dof <- n - length(mod2$be) - 2
        t2 <- ( mod2$be[2, ] / mod2$seb[2, ] )^2
        p2 <- pf(t2, 1, dof, lower.tail = FALSE, log.p = TRUE)
      }  ## end if ( mod2$info[5] >= mod0$info[5] ) 
    }  # end if ( Rfast::sort_unique.length(x) == 2 )
  }   ## end if ( is.null(cs) ) 
  
  pval <- min( log(2) + min(p1, p2), max(p1, p2) )
  stat <-  max( 2 * max(t1, t2), min(t1, t2) )
  result <- c(stat, pval, dof)
  names(result) <- c('test', 'logged.p-value', 'df') 
  result
}