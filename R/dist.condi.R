##### Conditional indpendence test for continuous variables using the distance correlation
#####
#####################
dist.condi <- function(ind1, ind2, cs, dat, type = NULL, rob = FALSE, R = 499) {
  ## ind1 and ind2 are the two indices of the two variables whose correlation is of interest
  ## cs is a vector with the indices of of variable(s),  over which the condition takes place
  ## dat is the data, a matrix form
  ## type is either "pearson" or "spearman"
  ## For a robust estimation of the PEarson correlation set rob = TRUE or FALSE otherwise
  d <- sum( cs>0 )  ## dimensionality of cs
  x1 <- dat[, ind1]
  x2 <- dat[, ind2 ]
    
  if ( d == 0 ) {  ## There are no conditioning variables
    mod <- energy::dcov.test(x1, x2, R = R)
    dof <- 1
    stat <- mod$statistic 
    pvalue <- log( mod$p.value )
  } else{  ## there are conditioning variables
    z <- dat[, cs]
    mod <- energy::pdcor.test(x1, x2, z, R)
    stat <- mod$statistic
    dof <- NCOL(z)
    pvalue <- log( mod$p.value )
  }
  #lets calculate the stat and p-value which are to be returned
  result <- c(stat, pvalue, dof)
  names(result) <- c('test', 'logged.p-value', 'df') 
 result
}
