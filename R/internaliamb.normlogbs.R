internaliamb.normlogbs <- function(target, dataset, threshold, wei, p) {
  n <- length(target)
  if ( !is.null(dataset) |  p > 0 ) {
    if ( p > 1 ) {
      ini <- glm( target ~.,  data = dataset, family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
      tab <- drop1( ini, test = "F" )
      dof <- tab[-1, 1]
      stat <- tab[-1, 4]    
    } else {
      ini <- glm( target ~.,  data = dataset, family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
      mod0 <- glm(target ~ 1, family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE)
      dof <- length( ini$coefficients ) - 1
      stat <- (mod0$deviance - ini$deviance) / dof / summary(ini)[[ 14 ]]
    }
    
    mat <- cbind(1:p, pf( stat, dof, n - dof, lower.tail = FALSE, log.p = TRUE), stat )
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
        mat <- mat <- matrix(nrow = 0, ncol = 3)
      } else if ( p - length(sel) == 1 ) {
        mod1 <- glm(target ~., data = dat, family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE)
        mod0 <- glm(target ~ 1, family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE)
        dof <- length( mod1$coefficients ) - 1
        stat <- abs( mod1$deviance - mod0$deviance )/dof/summary(mod1)[[ 14 ]]     
        pval <- pf( stat, dof, n - dof, lower.tail = FALSE, log.p = TRUE)
        if (pval > threshold ) {
          final <- "No variables were selected"
          mat <- matrix(nrow = 0, ncol = 3)
        } else final <- mod1
      } else {
        final <- glm(target ~., data = dat, family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE)
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