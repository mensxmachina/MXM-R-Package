internaliamb.tobitbs <- function(target, dataset, threshold, wei, p) {
  if ( !is.null(dataset) |  p > 0 ) {
    if ( p > 1 ) {
        ma <- survival::survreg( target ~.,  data = dataset, weights = wei, dist = "gaussian" )
        ini <- 2 * logLik(ma)
        dofini <- length( coef(ma) )
        stat <- dof <- numeric(p)
        for (i in 1:p) {
          mod <- survival::survreg( target ~.,  data = dataset[, -i, drop = FALSE], dist = "gaussian" )
          stat[i] <- ini - 2 * logLik(ma)
          dof[i] <- dofini - length( coef(mod) ) 
        }
      
    } else {
      ini <- survival::survreg( target ~.,  data = dataset, weights = wei, dist = "gaussian" )
      mod0 <- survival::survreg( target ~ 1,  data = dataset, weights = wei, dist = "gaussian" )
      stat <- 2 * logLik(ini) - 2 * logLik(mod0) 
      dof <- length( coef(ini) ) - 1
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
        mod1 <- survival::survreg( target ~ .,  data = dat, weights = wei, dist = "gaussian" )
        mod0 <- survival::survreg( target ~ 1,  data = dat, weights = wei, dist = "gaussian" )
        stat <- 2 * logLik(mod1) - 2 * logLik(mod0)
        pval <- pchisq( stat, length( coef(mod1) ) - 1, lower.tail = FALSE, log.p = TRUE)
        if (pval > threshold ) {
          final <- "No variables were selected"
          mat <- matrix(nrow = 0, ncol = 3)
        } else final <- mod1
      } else  final <- survival::survreg( target ~ .,  data = dat, weights = wei, dist = "gaussian" )
    }
    info <- info[ info[, 1] > 0, , drop = FALSE]
    
  } else { 
    info <- NULL  
    mat <- matrix(nrow = 0, ncol = 3)
    final <- "No variables were selected"
  } 
  
  list(info = info, mat = mat, final = final) 
}  