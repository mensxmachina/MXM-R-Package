ci.mm2 <- function(ind1, ind2, cs = NULL, suffStat) {
  
  dataset <- suffStat$dataset
  y <- dataset[, ind1]
  x <- dataset[, ind2]

  if ( is.null(cs) ) {
    if ( is.numeric(y) ) {
      mod1 <- lm(y ~., data = data.frame(x) )
      a1 <- anova(mod1)
      p1 <- a1[1, 5]
    } else if ( is.ordered(y) ) {
      ds <- data.frame(y = y, x = x)
      mod1 <- ordinal.reg(y ~., data = ds)
      mod0 <- ordinal.reg(y ~ 1, ds)
      t1 <-  mod0$devi - mod1$devi
      dof <- length( mod1$be ) - length( mod0$be )
      p1 <- pchisq(t1, dof, lower.tail = FALSE)
    } 
    ###### other direction now
    if ( is.numeric(x) ) {
      mod2 <- lm(x ~., data = data.frame(y) )
      a2 <- anova(mod2)
      p2 <- a2[1, 5]
    } else if ( is.ordered(x) ) {
      ds <- data.frame(x = x, y = y)
      mod2 <- ordinal.reg(x ~., data = ds) 
      mod0 <- ordinal.reg(x ~ 1, ds)
      t2 <- mod0$devi - mod2$devi
      dof <- length( mod2$be ) - length( mod0$be )
      p2 <- pchisq(t2, dof, lower.tail = FALSE)
    } 
    ### with conditioning set   
  } else {
    ########  
    if ( is.numeric(y) ) {
      ds <-  data.frame(y = y, dataset[, cs], x = x)
      mod1 <- lm(y ~., data = ds)
      ds0 <-  data.frame(y = y, dataset[, cs])
      mod0 <- lm(y ~., data = ds0)
      a1 <- anova(mod0, mod1)
      p1 <- a1[2, 6]
    } else if ( is.ordered(y)) {
      ds1 <-  data.frame(y = y, dataset[, cs], x = x)
      mod1 <- ordinal.reg(y ~., data = ds1)
      ds0 <- data.frame(y = y, dataset[, cs])
      mod0 <- ordinal.reg(y ~., data = ds0 )
      t1 <- mod0$devi - mod1$devi
      dof <- length( mod1$be ) - length( mod0$be )
      p1 <- pchisq(t1, dof, lower.tail = FALSE)
    } 
    ###### other direction now
    if ( is.numeric(x) ) {
      ds <-  data.frame(x = x, dataset[, cs], y = y)
      mod2 <- lm(x ~., data = ds)
      ds0 <-  data.frame(x = x, dataset[, cs])
      mod0 <- lm(x ~., data = ds0)
      a2 <- anova(mod0, mod2)
      p2 <- a2[2, 6]
    } else if ( is.ordered(x) ) {
      ds2 <- data.frame(x = x, dataset[, cs], y = y)
      mod2 <- ordinal.reg(x ~., data = ds2)
      ds0 <- data.frame(x = x, dataset[, cs])
      mod0 <- ordinal.reg(x ~., data = ds0 )
      t2 <- mod0$devi - mod2$devi
      dof <- length( mod2$be ) - length( mod0$be )
      p2 <- pchisq(t2, dof, lower.tail = FALSE)
    } 
  }     
  min( 2 * min(p1, p2), max(p1, p2) )
}

