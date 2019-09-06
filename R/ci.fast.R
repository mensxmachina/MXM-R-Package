ci.fast <- function(ind1, ind2, cs = NULL, dat, type, rob = FALSE, R = 1) {
  
  y <- dat[, ind1]
  x <- dat[, ind2]

  if ( is.null(cs) ) {
    ds <- data.frame(y = y, x = x)
    if ( length( unique(y) ) == 2 ) {  
      mod <- glm(y ~., data = data.frame(x), binomial)
      stat <- mod$null.deviance - mod$deviance
      dof <- length(mod$coefficients) - 1
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    } else if ( length( unique(x) ) == 2 ) {  
      mod <- glm(x ~., data = data.frame(y), binomial)
      stat <- mod$null.deviance - mod$deviance
      dof <- length(mod$coefficients) - 1
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.numeric(y) ) {
      mod <- lm(y ~., data = ds )  ## linear regression
      a1 <- anova(mod)
      stat <- a1[1, 4]
      dof <- a1[1, 1]
      pval <- pf(stat, dof, a1[2, 1], lower.tail = FALSE, log.p = TRUE)
    } else if ( is.numeric(x) ) {
      mod <- lm(x ~., data = ds )  ## linear regression
      a1 <- anova(mod)
      stat <- a1[1, 4]
      dof <- a1[1, 1]
      pval <- pf(stat, dof, a1[2, 1], lower.tail = FALSE, log.p = TRUE)
    } else if ( is.ordered(x) & is.ordered(y) ) {
      if ( length( Rfast::sort_unique(y) ) <= length( Rfast::sort_unique(x) ) ) {
         mod1 <- ordinal.reg(y ~., data = ds)  ## ordinal regression
         mod0 <- ordinal.reg(y ~ 1, ds)
      } else {
         mod1 <- ordinal.reg(x ~., data = ds)  ## ordinal regression
         mod0 <- ordinal.reg(x ~ 1, ds)
      }
      stat <-  mod0$devi - mod1$devi
      dof <- length( mod1$be ) - length( mod0$be )
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.factor(y) & !is.ordered(y) ) {
      mod1 <- nnet::multinom(y ~., data = ds, trace = FALSE)  ## multinomial regression
      mod0 <- nnet::multinom(y ~ 1., data = ds, trace = FALSE)
      a1 <- anova(mod0, mod1)
      stat <- a1[2, 6]
      dof <- a1[2, 5]
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.factor(x) & !is.ordered(x) ) {
      mod1 <- nnet::multinom(x ~., data = ds, trace = FALSE)  ## multinomial regression
      mod0 <- nnet::multinom(x ~ 1, data = ds, trace = FALSE)
      a1 <- anova(mod0, mod1)
      stat <- a1[2, 6]
      dof <- a1[2, 5]
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    }
    ### with conditioning set   
  } else {
    ########  
    ds0 <- dat[, cs, drop = FALSE]
    if ( length( unique(y) ) == 2) {
      ds <- data.frame(y = y, dat[, cs], x = x )
      mod1 <- glm(y ~., data = ds, binomial)
      mod0 <- glm(y ~ x, data = ds0, binomial)
      a1 <- anova(mod0, mod1)
      stat <- a1[2, 4]
      dof <- a1[2, 3]
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    } else if ( length( unique(x) ) == 2) {
      ds <- data.frame(x = x, dat[, cs], y = y )
      mod1 <- glm(x ~., data = ds, binomial)
      mod0 <- glm(x ~ y, data = ds0, binomial)
      a1 <- anova(mod0, mod1)
      stat <- a1[2, 4]
      dof <- a1[2, 3]
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.numeric(y) ) {
      ds <-  data.frame(y = y, dat[, cs], x = x)
      mod1 <- lm(y ~., data = ds)
      if ( any( is.na(mod1$coefficients) ) ) {
        pval <- log(1)
      } else {
        a1 <- anova(mod1)
        d1 <- dim(a1)[1] - 1
        stat <- a1[d1, 4]
        dof <- a1[d1, 1]   ;  df2 <- a1[d1 + 1, 1]
        pval <- pf(stat, dof, df2, lower.tail = FALSE, log.p = TRUE)
      }  
    } else if ( is.numeric(x) ) {
      ds <-  data.frame(x = x, dat[, cs], y = y)
      mod1 <- lm(x ~., data = ds)
      if ( any( is.na(mod1$coefficients) ) ) {
        pval <- log(1)
      } else {
        a1 <- anova(mod1)
        d1 <- dim(a1)[1] - 1
        stat <- a1[d1, 4]
        dof <- a1[d1, 1]   ;  df2 <- a1[d1 + 1, 1]
        pval <- pf(stat, dof, df2, lower.tail = FALSE, log.p = TRUE)
      }  
    } else if ( is.ordered(x) & is.ordered(y) ) {
      if ( length( Rfast::sort_unique(y) ) <= length( Rfast::sort_unique(x) ) ) {
        ds1 <- data.frame(y = y, dat[, cs], x = x)
        mod1 <- ordinal.reg(y ~., data = ds1)
        ds0 <- data.frame(y = y, dat[, cs])
        mod0 <- ordinal.reg(y ~., data = ds0 )
        stat <- mod0$devi - mod1$devi
        dof <- length( mod1$be ) - length( mod0$be )
        pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
      } else {
        ds1 <-  data.frame(x = x, dat[, cs], y = y)
        mod1 <- ordinal.reg(x ~., data = ds1)
        ds0 <- data.frame(x = x, dat[, cs])
        mod0 <- ordinal.reg(x ~., data = ds0 )
        stat <- mod0$devi - mod1$devi
        dof <- length( mod1$be ) - length( mod0$be )
        pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
      }  
    } else if ( is.factor(y)  &  !is.ordered(y) ) {
      ds1 <- data.frame(y = y, dat[, cs], x = x)
      mod1 <- nnet::multinom(y ~., data = ds1, trace = FALSE)
      ds0 <- data.frame(y = y, dat[, cs])
      mod0 <- nnet::multinom(y ~., data = ds0, trace = FALSE)
      stat <- mod0$deviance - mod1$deviance
      dof <- length( coef(mod1) ) - length( coef(mod0) )
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.factor(x)  &  !is.ordered(x) ) {
      ds1 <- data.frame(x = x, dat[, cs], y = y)
      mod1 <- nnet::multinom(y ~., data = ds1, trace = FALSE)
      ds0 <- data.frame(x = x, dat[, cs])
      mod0 <- nnet::multinom(x ~., data = ds0, trace = FALSE)
      stat <- mod0$deviance - mod1$deviance
      dof <- length( coef(mod1) ) - length( coef(mod0) )
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    }  
  }     
  result <- c(stat, pval, dof)
  names(result) <- c('test', 'logged.p-value', 'df') 
  result
}