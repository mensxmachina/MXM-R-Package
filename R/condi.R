#####################
#####################
##### Conditional indpendence test for continuous variables
#####
#####################
#####################
condi <- function(ind1, ind2, cs, dat, type = "pearson", rob = FALSE, R = 1) {
  ## ind1 and ind2 are the two indices of the two variables whose correlation is of interest
  ## cs is a vector with the indices of of variable(s),  over which the condition takes place
  ## dat is the data, a matrix form
  ## type is either "pearson" or "spearman"
  ## For a robust estimation of the Pearson correlation set rob = TRUE or FALSE otherwise
  n <- dim(dat)[1] ## sample size
  d <- sum( cs>0 )  ## dimensionality of cs
  if ( R == 1 ) {  ## no permutation
  ## NOTE: if you set test = "spearman", you must have the ranks of the data in 
  ## the dat argument and not the data themselves. This is to speed up the computations
  if (type == "spearman")  rob = FALSE   ## if spearman is chosen, robust is set to FALSE  

  if ( !rob ) {
    if ( d == 0 ) {
      r <- cor(dat[, ind1], dat[, ind2])
    } else {
      tmpm <- dat[, c(ind1, ind2, cs) ]
      tmpm <- as.matrix(tmpm)
      corrMatrix <- cor(tmpm)
      xyIdx <- 1:2
      csIdx <- 3:( d + 2 ) # or csIdx = 3
      residCorrMatrix <- corrMatrix[xyIdx, xyIdx] - as.matrix( corrMatrix[xyIdx, csIdx] ) %*% 
      ( solve( as.matrix( corrMatrix[csIdx, csIdx] ), rbind( corrMatrix[csIdx, xyIdx] ) ) ) 
      r <-  - residCorrMatrix[1, 2] / sqrt( residCorrMatrix[1, 1] * residCorrMatrix[2, 2]) 
      if ( abs(r) > 1 )   r <- 0.99999
    }
    
  } else {  ## robust estimation using M estimation
    if ( d == 0 ) {
      b1 <- coef( MASS::rlm( dat[, ind1] ~ dat[, ind2], maxit = 2000, method = "MM" ) )[2]
      b2 <- coef( MASS::rlm( dat[, ind2] ~ dat[, ind1], maxit = 2000, method = "MM" ) )[2]    
      r <- sqrt( abs(b1 * b2) )
    } else {
      e1 <- resid( MASS::rlm( dat[, ind1] ~.,  data = data.frame( dat[, c(ind2, cs) ] ), maxit = 2000, method = "MM") )
      e2 <- resid( MASS::rlm( dat[, ind2] ~., data = data.frame( dat[, c(ind1, cs) ] ), maxit = 2000, method = "MM" ) )
      r <- cor(e1, e2)
    }
  }
  
  if (type == "pearson") {
    stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(n - d - 3) )  ## absolute of the test statistic
  } else  stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(n - d - 3) ) /  1.029563  ## absolute of the test statistic
  
  pvalue <- log(2) + pt(stat, n - d - 3, lower.tail = FALSE, log.p = TRUE)  ## logged p-value
  result <- c(stat, pvalue, n - d - 3)
  names(result) <- c('test', 'logged.p-value', 'df')
  
  } else if ( R > 1 ) {  ## permutation based test
    
    x1 <- dat[, ind1]
    x2 <- dat[, ind2 ]
   
    if ( d == 0 ) {  ## There are no conditioning variables
      
      if ( rob ) { ## robust correlation

        e1 <- resid( MASS::rlm( x1 ~ 1, maxit = 2000, method = "MM" ) )
        e2 <- resid( MASS::rlm( x2 ~ 1, maxit = 2000, method = "MM" ) )       
        res <- Rfast::permcor( e1, e2, R )
        stat <- abs( res[1] )
        pvalue <- res[2]
        
      } else {
        res <- Rfast::permcor( x1, x2, R )
        stat <- abs( res[1] )
        pvalue <- res[2]
      }
      
    } else {  ## there are conditioning variables
      
      if ( rob ) { ## robust correlation
        e1 <- resid( MASS::rlm( x1 ~ ., data = data.frame(dat[, cs]), maxit = 2000, method = "MM") )
        e2 <- resid( MASS::rlm( x2 ~.,  data = data.frame(dat[, cs]), maxit = 2000, method = "MM") )
        res <- Rfast::permcor( e1, e2, R)
        stat <- abs( res[1] )
        pvalue <- res[2]
        
      } else {
        er <- .lm.fit(cbind(1, dat[, cs]), cbind( x1, x2 )  )$residuals
		res <- Rfast::permcor( er[, 1], er[, 2], R ) 
        stat <- abs( res[1] )
        pvalue <- (res[2])
      }
    }
    #lets calculate the stat and p-value which are to be returned
    dof <- n - d - 3; #degrees of freedom
    stat <- stat / dof
    pvalue <- log( pvalue )
    result <- c(stat, pvalue, dof)
    names(result) <- c('test', 'logged.p-value', 'df') 
  }
  
  result
}
