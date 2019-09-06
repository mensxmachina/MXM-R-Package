pearson_condis <-  function(ind1, ind2, cs1, cs2, Var, dat, type = "pearson", rob, max_k, R ) {
  
  d1 <- length(cs1)
  d2 <- length(cs2)
  
  if ( d1 == 0  &  d2 == 0) {
    pval <- condi(ind1, ind2, 0, dat = dat, type = type)[2]
    pval[2] <- condi(ind1, ind2, Var, dat = dat, type = type)[2]
    res <- cbind( c(0, 1), pval )
  } else {
  
    res1 <- matrix(nrow = 0, ncol = 2)   
    res2 <- matrix(nrow = 0, ncol = 2)   
  
    if (d1 > 0) {
      a <- list()
      pval <- list()
      cs <- c(cs1, Var)
      i <- 1
      a[[ i ]] <- combn(cs, i)
      while ( i < length(cs)  & i < max_k ) {
        i <- i + 1
        a[[ i ]] <- combn(cs, i)
      }  ## end while
    
      corrMatrix <- cor(dat)
      xyIdx <- c(ind1, ind2)
      n <- dim(dat)[1]
      ind <- list()
      for ( i in 1:length(a) ) {
        condset <- a[[ i ]]
        dm <- dim(condset)[2]
        a1 <- numeric(dm)
        a2 <- numeric(dm)
        dof <- numeric(dm)
        for ( j in 1:dm ) {
          dof <- n - j - 3
          csIdx <- condset[, j]
          a2[ j ] <- sum(csIdx == Var)
          residCorrMatrix <- corrMatrix[xyIdx, xyIdx] - as.matrix( corrMatrix[xyIdx, csIdx] ) %*% 
                            ( solve( as.matrix( corrMatrix[csIdx, csIdx] ), rbind( corrMatrix[csIdx, xyIdx] ) ) ) 
          r <-  residCorrMatrix[1, 2] / sqrt( residCorrMatrix[1, 1] * residCorrMatrix[2, 2]) 
          if ( abs(r) > 1 )   r <- 0.99999
          
          if (type == "pearson") {
            a1[j] <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(dof) )  ## absolute of the test statistic
          } else  a1[j] <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(dof) ) / 1.029563  ## absolute of the test statistic
        }  ## end for ( j in 1:dm )
        stat <- a1
        pval[[ i ]] <- log(2) + pt(stat, dof, lower.tail = FALSE, log.p = TRUE)  ## logged p-value
        ind[[ i ]] <- a2
      }  ## end for ( i in 1:length(a) ) 
      res1 <- cbind( unlist(ind), unlist(pval) )
    } ##  end if (d1 > 1)
  
    if (d2 > 0) {
      pval <- list()
      a <- list()
      cs <- c(cs2, Var)
      i <- 1
      a[[ i ]] <- combn(cs, i)
      while ( i < length(cs)  & i < max_k ) {
        i <- i + 1
        a[[ i ]] <- combn(cs, i)
      }  ## end while
    
      corrMatrix <- cor(dat)
      xyIdx <- c(ind1, ind2)
      n <- dim(dat)[1]
      ind <- list()
      for ( i in 1:length(a) ) {
        condset <- a[[ i ]]
        dm <- dim(condset)[2]
        a1 <- numeric(dm)
        a2 <- numeric(dm)
        dof <- numeric(dm)
        for ( j in 1:dm ) {
          dof <- n - j - 3
          csIdx <- condset[, j]
          a2[j] <- sum(csIdx == Var)
          residCorrMatrix <- corrMatrix[xyIdx, xyIdx] - as.matrix( corrMatrix[xyIdx, csIdx] ) %*% 
          ( solve( as.matrix( corrMatrix[csIdx, csIdx] ), rbind( corrMatrix[csIdx, xyIdx] ) ) ) 
          r <-  residCorrMatrix[1, 2] / sqrt( residCorrMatrix[1, 1] * residCorrMatrix[2, 2]) 
          if ( abs(r) > 1 )   r <- 0.99999
        
          if (type == "pearson") {
            a1[j] <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(dof) )  ## absolute of the test statistic
          } else  a1[j] <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(dof) ) / 1.029563  ## absolute of the test statistic
        }  ## end for ( j in 1:dm )
        stat <- a1
        pval[[ i ]] <- log(2) + pt(stat, dof, lower.tail = FALSE, log.p = TRUE)  ## logged p-value
        ind[[ i ]] <- a2
      }  ## end for ( i in 1:length(a) ) { 
      res2 <- cbind( unlist(ind), unlist(pval) )
    } ##  end if (d2 > 1)
    res <- rbind(res1, res2)
  
  }  ## end  if (d1 == 0  &  d2 == 0)  
  colnames(res) <- c("Var", "log.pvalue")
  res
}  
  
  