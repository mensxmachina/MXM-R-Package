cat_condis <-  function(ind1, ind2, cs1, cs2, Var, dat, type, rob = FALSE, max_k, R = 1 ){
  
  d1 <- sum(cs1 > 0 )
  d2 <- sum(cs2 > 0)
  
  if (d1 == 0  &  d2 == 0) {
    pval <- cat.ci(ind1, ind2, Var, 0, type = type, rob = FALSE, R = R)[2] 
    pval[2] <- cat.ci(ind1, ind2, Var, dat, type = type, rob = FALSE, R = R)[2] 
    res <- cbind( c(0, 1), pval )
  
  } else {  ## end if ( d1 == 0  &  d2 == 0)  and their else
    
    res1 <- matrix(nrow = 0, ncol = 2)   
    res2 <- matrix(nrow = 0, ncol = 2)   
    
    if (d1 > 0) {
      pval <- list()
      cs <- c(cs1, Var)
      i <- 1
      a[[ i ]] <- combn(cs, i)
      while ( i < length(cs)  & i < max_k ) {
        i <- i + 1
        a[[ i ]] <- combn(cs, i)
      }  ## end while
      ind <- list()
      for ( i in 1:length(a) ) {
        condset <- a[[ i ]]
        dm <- dim(condset)[2]
        a1 <- numeric(dm)
        a2 <- numeric(dm)
        for ( j in 1:dm ) {
          csIdx <- condset[, j]
          a2[j] <- sum(csIdx == Var)
          a1[j] <- cat.ci(ind1, ind2, csIdx, dat, type = type, rob = FALSE, R = R)[2] 
        }  ## end for ( j in 1:dm )
        pval[[ i ]] <- a1
        ind[[ i ]] <- a2
      }  ## end for ( i in 1:length(a) ) { 
      res1 <- cbind( unlist(ind), unlist(pval) )
    } ##  end if (d1 > 1)
    
    if (d2 > 0) {
      a <- list()
      pval <- list()
      cs <- c(cs2, Var)
      i <- 1
      a[[ i ]] <- combn(cs, i)
      while ( i < length(cs)  &  i < max_k ) {
        i <- i + 1
        a[[ i ]] <- combn(cs, i)
      }  ## end while
      ind <- list()
      for ( i in 1:length(a) ) {
        condset <- a[[ i ]]
        dm <- dim(condset)[2]
        a1 <- numeric(dm)
        a2 <- numeric(dm)
        for ( j in 1:dm ) {
          csIdx <- condset[, j]
          a2[j] <- sum(csIdx == Var)
          a1[j] <- cat.ci(ind1, ind2, csIdx, dat, type = type, rob = FALSE, R = R)[2] 
        }  ## end for ( j in 1:dm )
        pval[[ i ]] <- a1
        ind[[ i ]] <- a2
      }  ## end for ( i in 1:length(a) ) { 
      res2 <- cbind( unlist(ind), unlist(pval) )
    } ##  end if (d1 > 1)
    res <- rbind(res1, res2)
    
  }  ## end  if (d1 == 0  &  d2 == 0)  
  colnames(res) <- c("Var", "log.pvalue")
  res
}  

