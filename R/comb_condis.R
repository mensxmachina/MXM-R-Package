comb_condis <- function(ind1, ind2, cs1, cs2, Var, dat, type = "ci.mm", rob, max_k, R ) {
  
  d1 <- length(cs1)
  d2 <- length(cs2)
  if (type == "ci.mm") {
    test <- ci.mm
  } else  test <- ci.fast  
  
  if ( d1 == 0  &  d2 == 0) {
    pval <- test(ind1, ind2, cs = 0, dat = dat)[2]
    pval[2] <- test(ind1, ind2, cs = Var, dat = dat)[2]
    res <- cbind(c(0, 1), pval)
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
      ind <- list()
      for ( i in 1:length(a) ) {
        condset <- a[[ i ]]
        dm <- dim(condset)[2]
        a1 <- numeric(dm)
        a2 <- numeric(dm)
        for ( j in 1:dm ) {
          csIdx <- condset[, j]
          a2[j] <- sum(csIdx == Var)
          a1[j] <- test(ind1, ind2, csIdx, dat = dat)[2]
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
          a1[j] <- condi(ind1, ind2, csIdx, dat = dat)[2]
        }  ## end for ( j in 1:dm )
        pval[[ i ]] <- a1
        ind[[ i ]] <- a2
      }  ## end for ( i in 1:length(a) ) { 
      res2 <- cbind( unlist(ind), unlist(pval) )
    } ##  end if (d2 > 1)
    res <- rbind(res1, res2)
  
  }  ## end  if (d1 == 0  &  d2 == 0)  
  colnames(res) <- c("Var", "log.pvalue")
  res
}  

