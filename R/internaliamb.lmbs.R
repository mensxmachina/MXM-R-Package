######## Internal for linear regression
internaliamb.lmbs <- function(target, dataset, threshold, wei, p) {
  if ( !is.null(dataset) |  p > 0 ) {
    n <- length(target)
    if ( p > 1 ) {
       ini <- lm( target ~., data = dataset, weights = wei, y = FALSE, model = FALSE )
       df2 <- n - length( coef(ini) )
       tab <- drop1( ini, test = "F" )
       df1 <- tab[-1, 1]
       stat <- tab[-1, 5]
       mat <- cbind(1:p, pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE), stat )
    } else {
      ini <- lm( target ~., data = dataset, weights = wei, y = FALSE, model = FALSE )
      df2 <- n - length( coef(ini) )
      tab <- drop1( ini, test = "F" )
      df1 <- tab[-1, 1]
      stat <- tab[-1, 5]
      mat <- cbind(1, pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE), stat )  
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
        ini <- lm( target ~., data = dataset, weights = wei, y = FALSE, model = FALSE )
        df2 <- n - length( coef(ini) )
        tab <- drop1( ini, test = "F" )
        df1 <- tab[-1, 1]
        stat <- tab[-1, 5]
        pval <- pf(stat, df1, df2, lower.tail = FALSE, log.p = TRUE)

        if (pval > threshold ) {
          final <- "No variables were selected"
          mat <- matrix(nrow = 0, ncol = 3)
        } else final <- ini
      } else {
        final <- lm( target ~., data = dataset, weights = wei, y = FALSE, model = FALSE )        
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