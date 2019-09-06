ci.mm <- function(ind1, ind2, cs = NULL, dat, type, rob = FALSE, R = 1) {

  y <- dat[, ind1]
  x <- dat[, ind2]

  if ( is.null(cs) ) {
    ds <- data.frame(y = y, x = x)
    if ( length( unique(y) ) == 2 ) {  
      mod1 <- glm(y ~., data = ds, binomial) ## logistic regression
      t1 <- mod1$null.deviance - mod1$deviance
      dof1 <- length(mod1$coefficients) - 1
      p1 <- pchisq(t1, dof1, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.numeric(y) ) {
      mod1 <- lm(y ~., data = ds )  ## linear regression
      a1 <- anova(mod1)
      t1 <- a1[1, 4]
      dof1 <- a1[1, 1]
      p1 <- pf(t1, dof1, a1[2, 1], lower.tail = FALSE, log.p = TRUE)
    } else if ( is.ordered(y) ) {
      mod1 <- MXM::ordinal.reg(y ~., data = ds)  ## ordinal regression
      mod0 <- MXM::ordinal.reg(y ~ 1, data = ds)
      t1 <-  mod0$devi - mod1$devi
      dof1 <- abs( length( mod1$be ) - length( mod0$be ) )
      p1 <- pchisq(t1, dof1, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.factor(y)  &  !is.ordered(y) ) {  ## multinomial regression
      mod1 <- nnet::multinom(y ~., data = ds, trace = FALSE)
      mod0 <- nnet::multinom(y ~ 1, data = ds, trace = FALSE)
      a1 <- anova(mod1, mod0)
      t1 <- a1[2, 6]
      dof1 <- a1[2, 5]
      p1 <- pchisq(t1, dof1, lower.tail = FALSE, log.p = TRUE)
    } 
    ###### other direction now
    if ( length( unique(x) ) == 2) {
      mod2 <- glm(x ~., data = ds, binomial)  ## logistic regression
      t2 <- mod2$null.deviance - mod2$deviance
      dof2 <- length(mod2$coefficients) - 1
      p2 <- pchisq(t2, dof2, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.numeric(x) ) {
      mod2 <- lm(x ~., data = ds ) ## linear regression
      a2 <- anova(mod2)
      t2 <- a2[1, 4]
      dof2 <- a2[1, 1]
      p2 <- pf(t2, dof2, a2[2, 1], lower.tail = FALSE, log.p = TRUE)
    } else if ( is.ordered(x) ) {
      mod2 <- MXM::ordinal.reg(x ~., data = ds)  ## ordinal regression
      mod0 <- MXM::ordinal.reg(x ~ 1, data = ds)
      t2 <- mod0$devi - mod2$devi
      dof2 <- abs( length( mod2$be ) - length( mod0$be ) )
      p2 <- pchisq(t2, dof2, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.factor(x)  &  !is.ordered(x) ) {
      mod2 <- nnet::multinom(x ~., data = ds, trace = FALSE)  ## multinomial regression
      mod0 <- nnet::multinom(x ~ 1, data = ds, trace = FALSE)
      a2 <- anova(mod2, mod0)
      t2 <- a2[2, 6]
      dof2 <- a2[2, 5]
      p2 <- pchisq(t2, dof2, lower.tail = FALSE, log.p = TRUE)
    } 
  } else {
    ### with conditioning set   
    ds0 <- dat[, cs, drop = FALSE ]
    ds1 <- dat[, c(cs, ind2) ]
    if ( length( unique(y) ) == 2 ) { 
      mod1 <- glm(y ~., data = ds1, binomial)  ## logistic regression
      mod0 <- glm(y ~., data = ds0, binomial)
      a1 <- anova(mod0, mod1)
      t1 <- a1[2, 4]
      dof1 <- a1[2, 3]
      p1 <- pchisq(t1, dof1, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.numeric(y) ) {
      mod1 <- lm(y ~., data = ds1)
      a1 <- anova(mod1)
      if ( any( is.na(mod1$coefficients) ) ) {
        p1 <- log(1)
      } else {
        a1 <- anova(mod1)
        d1 <- dim(a1)[1] - 1
        t1 <- a1[d1, 4]
        dof1 <- a1[d1, 1]   ;  df2 <- a1[d1 + 1, 1]
        p1 <- pf(t1, dof1, df2, lower.tail = FALSE, log.p = TRUE)
      }  
    } else if ( is.ordered(y) ) {
      mod1 <- MXM::ordinal.reg(y ~., data = ds1)   ## ordinal regression
      mod0 <- MXM::ordinal.reg(y ~., data = ds0)
      t1 <- mod0$devi - mod1$devi
      if (t1 <0 )  t1 <- 0
      dof1 <- abs( length( mod1$be ) - length( mod0$be ) )
      p1 <- pchisq(t1, dof1, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.factor(y)  &  !is.ordered(y) ) {
      mod1 <- nnet::multinom(y ~., data = ds1, trace = FALSE)  ## multinomial regression
      mod0 <- nnet::multinom(y ~., data = ds0, trace = FALSE)
      t1 <- deviance(mod0) - deviance(mod1)
      dof1 <- length( coef(mod1) ) - length( coef(mod0) )
      p1 <- pchisq(t1, dof1, lower.tail = FALSE, log.p = TRUE)
    } 
    ###### other direction now
    ds2 <- dat[, c(cs, ind1) ]
    if ( length( unique(x) ) == 2 ) {
      mod2 <- glm(x ~., data = ds2, binomial)  ## binomial regression
      mod0 <- glm(x ~., data = ds0, binomial)
      a2 <- anova(mod0, mod2)
      t2 <- a2[2, 4]
      dof2 <- a2[2, 3]
      p2 <- pchisq(t2, dof2, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.numeric(x) ) {
      mod2 <- lm(x ~., data = ds2)  ## linear regression
      if ( any( is.na(mod2$coefficients) ) ) {
        p2 <- log(1)
      } else {
        a2 <- anova(mod2)
        d2 <- dim(a2)[1] - 1
        t2 <- a2[d2, 4]
        dof2 <- a2[d2, 1]   ;  df2 <- a2[d2 + 1, 1]
        p2 <- pf(t2, dof2, df2, lower.tail = FALSE, log.p = TRUE)
      }  
    } else if ( is.ordered(x) ) {
      mod2 <- MXM::ordinal.reg(x ~., data = ds2)  ## ordinal regression
      mod0 <- MXM::ordinal.reg(x ~., data = ds0)
      t2 <- mod0$devi - mod2$devi
      if (t2 <0 )  t2 <- 0
      dof2 <- abs( length( mod2$be ) - length( mod0$be ) )
      p2 <- pchisq(t2, dof2, lower.tail = FALSE, log.p = TRUE)
    } else if ( is.factor(x)  &  !is.ordered(x) ) {
      mod2 <- nnet::multinom(x ~., data = ds2, trace = FALSE)  ## multinomial regression
      mod0 <- nnet::multinom(x ~., data = ds0, trace = FALSE)
      t2 <- deviance(mod0) - deviance(mod2)
      dof2 <- length( coef(mod2) ) - length( coef(mod0) )
      p2 <- pchisq(t2, dof2, lower.tail = FALSE, log.p = TRUE)
    } 
  }   
  pval <- min( log(2) + min(p1, p2), max(p1, p2) )
  stat <-  max( 2 * max(t1, t2), min(t1, t2) )
  if (p1 < p2) { 
    dof <- dof1
  } else  dof <- dof2
  result <- c(stat, pval, dof)
  names(result) <- c('test', 'logged.p-value', 'df') 
  result
}
