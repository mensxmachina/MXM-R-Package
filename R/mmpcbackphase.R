mmpcbackphase <- function(target, dataset, max_k = 3, threshold = 0.05, test = NULL, wei = NULL, R = 1 ) {
  d <- NCOL(dataset)
  met <- 1:d
  counter <- 0
  
  if (R == 1) {
    
    threshold = log(threshold)
    
    if (d == 1) {
      tes <- test(target, dataset, 1, 0, wei = wei) 
      pv <- tes$pvalue
      counter <- 1
      if ( pv > threshold )  met <- 0
      pvalue <- pv
    
    } else  if ( d == 2  ) {
      tes1 <- test(target, dataset, 1, 2, wei = wei)
      tes2 <- test(target, dataset, 2, 1, wei = wei)
      pv1 <- tes1$pvalue
      pv2 <- tes2$pvalue
      counter <- 2
      if ( pv1 > threshold )  met[1] <- 0  
      if ( pv2 > threshold )  met[2] <- 0  
      pvalue <- c(pv1, pv2)
      
    } else {  
      pvalue <- 0
      for (i in 1:d) {
        k <- 0
        pval <-  -5
        b <- 0
        while ( k < d - 1  &  k < max_k  &  pval < threshold ) {
          k <- k + 1
          ze <- combn(d, k)
          rem <- which(ze == i,arr.ind = TRUE)[, 2]
          ze <- ze[, -rem, drop = FALSE]
          j <- 0
          while ( j < dim(ze)[2]  &  pval < threshold ) {
            j <- j + 1
            tes <- test(target, dataset, i, ze[, j], wei = wei)
            pval <- tes$pvalue
            b[j] <- pval
            counter <- counter + 1
          }  ## end while ( j < dim(ze)[2]  &  pval < threshold )
        }  ##end  while ( k < d - 1  &  k < max_k & pval < threshold )
        if ( pval > threshold )  met[i] <- 0
        pvalue[i] <- max(b)
      }  ## end for (i in 1:d)
    }  ## end if (d == 1) 
    
  } else {
    
    if (d == 1) {
      tes <- test(target, dataset, 1, 2, wei = wei, threshold = threshold, R = R)
      pv <- tes$pvalue
      counter <- 1
      if ( pv > threshold )  met <- 0
      pvalue <- pv
      
    } else if ( d == 2  ) {
      tes1 <- test(target, dataset, 1, 2, wei = wei, threshold = threshold, R = R)
      tes2 <- test(target, dataset, 2, 1, wei = wei, threshold = threshold, R = R)
      pv1 <- tes1$pvalue
      pv2 <- tes2$pvalue
      counter <- 2
      if ( pv1 > threshold )  met[1] <- 0  
      if ( pv2 > threshold )  met[2] <- 0  
      pvalue <- c(pv1, pv2)
      
    } else {  
      pvalue <- 0
      for (i in 1:d) {
        k <- 0
        pval <-  0
        b <- 0
        while ( k < d - 1  &  k < max_k & pval < threshold ) {
          k <- k + 1
          ze <- combn(d, k)
          rem <- which(ze == i,arr.ind = TRUE)[, 2]
          ze <- ze[, -rem, drop = FALSE]
          j <- 0
          while ( j < dim(ze)[2]  &  pval < threshold ) {
            j <- j + 1
            tes <- test(target, dataset, i, ze[, j], wei = wei, threshold = threshold, R = R)
            pval <- tes$pvalue
            b[j] <- pval
            counter <- counter + 1
          }  ## end while ( j < dim(ze)[2]  &  pval < threshold )
        }  ## end while ( k < d - 1  &  k < max_k & pval < threshold )
        if ( pval > threshold )   met[i] <- 0
        pvalue[i] <- max(b)
      }  ## end for (i in 1:d)
    }  ## end if (d == 1) 
  }  ## end if (R == 1)
  list(met = met, counter = counter, pvalues = pvalue)
}