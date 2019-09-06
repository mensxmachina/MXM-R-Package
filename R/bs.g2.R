bs.g2 <- function(target, dataset, threshold = 0.05) {
  
  runtime <- proc.time() 
  dataset <- as.matrix(dataset)
  target <- as.numeric(target)
  z <- cbind(dataset, target)
  all.dc <- Rfast::colrange(z, cont = FALSE)
  p <- dim(z)[2] - 1
  ind <- 1:p
  stat <- rep( Inf, p )
  pval <- rep( -100, p )
  dof <- rep( 100, p )
  sig <- log(threshold)
  
  for ( i in 1:p )  {
    mod <- Rfast::g2Test(z, x = p + 1, y = i, cs = ind[ ind != i ], dc = all.dc) 
    dof[i] <- mod$df
    stat[i] <- mod$statistic
  }  ## end for (i in 2:p)  
  pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
  no <- which(pval == 0)
  if ( length(no) > 1 ) { 
    sel <- which.min( stat[no]/dof[no] )
  } else sel <- which.max(pval)
  
  mat <- cbind(1:p, pval, stat)
  colnames(mat) <- c("variable", "log.p-values", "statistic" )
  info <- matrix( c(0, -10, -10), ncol = 3 )
  colnames(info) <- c("Variables", "log.p-values", "statistic")
  
  if ( pval[sel] < sig ) {
    res <- list(info = matrix(0, 0, 3), runtime = proc.time() - runtime, mat = mat, ci_test = "gSquare") 
    
  } else {
    
    info[1, ] <- mat[sel, ]
    mat <- mat[-sel, , drop = FALSE] 
    dat <- z[, -sel, drop = FALSE] 
    dc2 <- all.dc[-sel]
    p <- p - 1
    ind <- 1:p
    stat <- rep( Inf, p )
    pval <- rep( -100, p )
    dof <- rep( 100, p )
    for ( i in 1:p )  {
      mod <- Rfast::g2Test(dat, x = p + 1, y = i, cs = ind[ ind != i ], dc = dc2)           
      stat[i] <- mod$statistic
      dof[i] <- mod$df
    }  ## end for (i in 2:p)  
    pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    mat[, 2:3] <- cbind(pval, stat)
    no <- which(pval == 0)
    if ( length(no) < 1 ) { 
      sel <- which.max( pval )
    } else sel <- which.min(stat[no]/dof[no])
    
    while ( pval[sel] > sig  &  p > 1 ) {
      info <- rbind(info, mat[sel, ])
      mat <- mat[-sel, , drop = FALSE] 
      dat <- z[, -info[, 1], drop = FALSE] 
      dc2 <- all.dc[ -info[, 1] ]
      p <- p - 1
      
      if ( p == 1 ) {
        mod <- Rfast::gchi2Test(target, dat[, -(p + 1)], logged = TRUE) 
        stat <- mod[2, 1]
        pval <- mod[2, 2]
        sel <- 1
        if ( pval > sig ) {
          info <- rbind(info, mat[sel, ])
          mat <- NULL
        }
        
      } else {
      
        ind <- 1:p
        stat <- rep( Inf, p )
        pval <- rep( -100, p )
        dof <- rep( 100, p )

        for ( i in 1:p )  {
          mod <- Rfast::g2Test(dat, x = p + 1, y = i, cs = ind[ ind != i ], dc = dc2)           
          stat[i] <- mod$statistic
          dof[i] <- mod$df
        }  ## end for (i in 2:p)     
        pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        mat[, 2:3] <- cbind(pval, stat)
        no <- which(pval == 0)
        if ( length(no) < 1 ) { 
          sel <- which.max( pval )
        } else sel <- which.min(stat[no]/dof[no])

      }  ## end else of if (p == 1)
      
    }  ## end  while ( pval[sel] > sig  &  p > 1 )
    runtime <- proc.time() - runtime
    info <- info[ info[, 1] > 0, , drop = FALSE]
    res <- list(runtime = runtime, info = info, mat = mat, ci_test = "gSquare" ) 
  }  ## end else of if ( pval[sel] < sig ) 
  
  res
}  
  
  
  