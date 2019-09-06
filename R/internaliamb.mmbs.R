######## Internal for linear regression
internaliamb.mmbs <- function(target, dataset, threshold, wei, p) {
  if ( !is.null(dataset) |  p > 0 ) {
    n <- length(target)
    if ( p > 1 ) {
      ini <- MASS::rlm( target ~., data = dataset, weights = wei, maxit = 2000, method = "MM" )
      lik1 <- as.numeric( logLik(ini) )
      dofini <- length( coef(ini) ) 
      for (i in 1:p) {
        fit2 <- MASS::rlm( target ~., data = dataset[, -i, drop = FALSE], weights = wei, maxit = 2000, method = "MM" )
        stat[i] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
        dof[i] <- dofini - length( coef(fit2) )
      }
      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
    } else {
      ini <- MASS::rlm( target ~., data = dataset, weights = wei, maxit = 2000, method = "MM" )
      mod0 <- MASS::rlm(target ~ 1, weights = wei, maxit = 2000, method = "MM" ) 
      stat <- 2 * abs( as.numeric( logLik(ini) ) - as.numeric( logLik(mod0) ) )
      mat <- cbind(1, pchisq( stat, length( coef(ini) ) - 1, lower.tail = FALSE, log.p = TRUE), stat )
    }
    
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
        ini <- MASS::rlm( target ~., data = dataset, weights = wei, maxit = 2000, method = "MM" )
        mod0 <- MASS::rlm( target ~ 1, weights = wei, maxit = 2000, method = "MM" ) 
        stat <- 2 * abs( as.numeric( logLik(ini) ) - as.numeric( logLik(mod0) ) )
        pval <- pchisq(stat, length( coef(ini) ) - 1, lower.tail = FALSE, log.p = TRUE)
        if (pval > threshold ) {
          final <- "No variables were selected"
          mat <- matrix(nrow = 0, ncol = 3)
        } else final <- ini
      } else {
        final <- MASS::rlm( target ~., data = dataset, weights = wei, maxit = 2000, method = "MM" )
      }
    }
    info <- info[ info[, 1] > 0, , drop = FALSE]
    
  } else { 
    info <- NULL  
    mat <- matrix(nrow = 0, ncol = 3) 
    final <- "No variables were selected"
  } 
  
  list(info = info, mat = mat, final = final) 
} 