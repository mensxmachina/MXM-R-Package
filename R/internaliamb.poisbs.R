######## Internal for Poisson regression
internaliamb.poisbs <- function(target, dataset, threshold, wei, p) {
  if ( !is.null(dataset) |  p > 0 ) {
    if ( p > 1 ) {
      ini <- glm( target ~.,  data = dataset, family = poisson(log), weights = wei, y = FALSE, model = FALSE )
      tab <- drop1( ini, test = "Chisq" )
      dof <- tab[-1, 1]
      stat <- tab[-1, 4]
    } else {
      ini <- glm( target ~.,  data = dataset, family = poisson(log), weights = wei, y = FALSE, model = FALSE )
      mod0 <- glm(target ~ 1, family = poisson(log), weights = wei, y = FALSE, model = FALSE)
      stat <- mod0$deviance - ini$deviance
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
        mat <- NULL
      } else if ( p - length(sel) == 1 ) {
        mod1 <- glm(target ~., data = dat, family = poisson(log), weights = wei, y = FALSE, model = FALSE)
        mod0 <- glm(target ~ 1, family = poisson(log), weights = wei, y = FALSE, model = FALSE)
        stat <- abs( mod1$deviance - mod0$deviance )
        pval <- pchisq( stat, length( coef(mod1) ) - 1, lower.tail = FALSE, log.p = TRUE)
        if (pval > threshold ) {
          final <- "No variables were selected"
          mat <- matrix(nrow = 0, ncol = 3)
        } else final <- mod1
      } else {
        final <- glm(target ~., data = dat, family = poisson(log), weights = wei, y = FALSE, model = FALSE)
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