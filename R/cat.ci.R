####################
#### G^2 (and X^2) test of (un)conditional independence
####################
cat.ci <- function(ind1, ind2, cs, dat, type, rob = FALSE, R = 1) {
  ## the ind1 and ind2 are two numbers, 1 and 2 for example
  ## indicating the two variables whose conditional independence 
  ## will be tested
  ## ind1, ind2 and cs must be different, non onverlapping numbers
  ## cs is one or more numbers indicating the conditioning variable(s)
  ## it is et to 0 by default. In this case an uncodntional test of 
  ## independence is  performed
  ## dat is the whole dat, and is expected to be a matrix
  ## type is set to NULL be default. This argument is not taken into consideration anywhere
  ## rob is FALSE by default, even if it is TRUE it is not taken into cosideration
  ## the type and rob arguments are put here so as to have the same signature as condi
  if ( R == 1 ) {
    if ( sum(cs == 0) > 0 ) {  ## There are no conditioning variables
      a1 <- Rfast::g2Test_univariate(dat[, c(ind1, ind2)], type)
      stat <- as.numeric( a1$statistic )
      dof <- as.numeric( a1$df )
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
      if ( dim(dat)[1] < 5 * dof ) {  ## condition to perform the test
        a1 <- Rfast::g2Test_univariate_perm(dat[, c(ind1, ind2)], type, 1000)           
        stat <- as.numeric( a1$statistic )
        dof <- as.numeric( a1$df )
        pval <- log(a1$pvalue)
      } 
    } else {   ## There are conditioning variables
      a1 <- Rfast::g2Test(dat, ind1, ind2, cs, type)           
      stat <- as.numeric( a1$statistic )
      dof <- as.numeric( a1$df )
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
      if ( dim(dat)[1] < 5 * dof ) {  ## condition to perform the test
        a1 <- Rfast::g2Test_perm(dat, ind1, ind2, cs, type, 1000)           
        stat <- as.numeric( a1$statistic )
        dof <- as.numeric( a1$df )
        pval <- log(a1$pvalue)
      } 
    }
  } else {
    if ( sum(cs == 0) > 0 ) {  ## There are no conditioning variables
      a1 <- Rfast::g2Test_univariate_perm(dat[, c(ind1, ind2)], type, R)           
      stat <- as.numeric( a1$statistic )
      dof <- NULL
      pval <- a1$pvalue
    } else {   ## There are conditioning variables
      a1 <- Rfast::g2Test_perm(dat, ind1, ind2, cs, type, R)           
      stat <- as.numeric( a1$statistic )
      dof <- as.numeric(a1$df)
      pval <- log(a1$pvalue)
    }
  }

  res <- c( stat, pval, dof )  
  names(res) <- c("Chi-squared test", "logged p-value", "df")
  res
}