internaliamb.bs <- function(target, dataset, threshold, test, wei, p) {
  if ( !is.null(dataset) |  p > 0 ) {
    ini <- test( target ~.,  data = dataset, weights = wei )
    dofini <- length( coef(ini) )
    likini <- logLik(ini) 
    stat <- dof <- numeric(p)
    
    if (p > 1) {
      for (j in 1:p) {
        mod <- test( target ~.,  data = dataset[, -j, drop = FALSE], weights = wei )
        stat[j] <- 2 * abs( likini - logLik(mod) )
        dof[j] <- dofini - length( coef(mod) ) 
      }
    } else {
      mod <- test( target ~ 1,  data = dataset, weights = wei )
      stat <- 2 * abs( likini - logLik(mod) )
      dof <- dofini - length( coef(mod) ) 
    }
    
    mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
    colnames(mat) <- c("variable", "log.p-values", "statistic" )
    rownames(mat) <- 1:p  
    info <- matrix( c(0, -10, -10) , ncol = 3 )
    sel <- which( mat[, 2] > threshold ) 
    
    if ( length(sel) == 0 ) {
      final <- ini 
      
    } else {
      
      info <- mat[sel, , drop = FALSE]
      mat <- mat[-sel, , drop = FALSE] 
      dat <- dataset[, -sel, drop = FALSE] 
      
      if ( p - length(sel) == 0 ) {
        final <- "No variables were selected"
        mat <- matrix(nrow = 0, ncol = 3)
      } else if ( p - length(sel) == 1 ) {
        mod1 <- test(target ~ ., data = dat, weights = wei )
        mod0 <- test(target ~ 1, weights = wei)
        stat <- 2 * abs( logLik(mod1) - logLik(mod0) )
        dof <- length( coef(mod1) ) - length( coef(mod0) ) 
        pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
        if (pval > threshold ) {
          final <- "No variables were selected"
          mat <- matrix(nrow = 0, ncol = 3)
        } else final <- mod1
      } else  final <- test(target ~ ., data = dat, weights= wei )
    }
    info <- info[ info[, 1] > 0, , drop = FALSE]
    
  } else { 
    info <- NULL  
    mat <- matrix(nrow = 0, ncol = 3) 
    final <- "No variables were selected"
  } 
  
  list(info = info, mat = mat, final = final) 
}  
