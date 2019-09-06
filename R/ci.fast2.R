ci.fast2 <- function(ind1, ind2, cs = NULL, suffStat) {
  
  dataset <- suffStat$dataset
  y <- dataset[, ind1]
  x <- dataset[, ind2]
  ########  
  if( is.numeric(y) ) {
    ds <- data.frame(y = y, dataset[, cs], x = x)
    mod1 <- lm(y ~., data = ds)
    ds0 <-  data.frame(y = y, dataset[, cs])
    mod0 <- lm(y ~., data = ds0)
    a1 <- anova(mod0, mod1)
    pval <- a1[2, 6]
  } else if ( is.numeric(x) )  {
    ds <-  data.frame(x = x, dataset[, cs], y = y)
    mod1 <- lm(x ~., data = ds)
    ds0 <-  data.frame(x = x, dataset[, cs])
    mod0 <- lm(x ~., data = ds0)
    a1 <- anova(mod0, mod1)
    pval <- a1[2, 6]
  } else if  ( is.ordered(y)  &  is.ordered(x) )  {
    if ( length( unique(y) ) <= length( unique(x) ) ) {
      ds1 <-  data.frame(y = y, dataset[, cs], x = x )
      mod1 <- ordinal.reg(y ~., data = ds1)
      ds0 <- data.frame(y = y, dataset[, cs])
      mod0 <- ordinal.reg(y ~ ., data = ds0 )
    } else {
      ds1 <-  data.frame(x = x, dataset[, cs], y = y )
      mod1 <- ordinal.reg(x ~., data = ds1)
      ds0 <- data.frame(x = x, dataset[, cs])
      mod0 <- ordinal.reg(x ~ ., data = ds0 )
    }  
    stat <- mod0$devi - mod1$devi
    dof <- length( mod1$be ) - length( mod0$be )
    pval <- pchisq(stat, dof, lower.tail = FALSE)
  } 
  pval
}
