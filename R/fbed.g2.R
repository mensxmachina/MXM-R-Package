fbed.g2 <- function(y, x, alpha = 0.05, univ = NULL, K = 0, backward = TRUE) { 
  p <- dim(x)[2]
  ind <- 1:p
  sig <- log(alpha)
  sela <- NULL
  card <- 0
  sa <- NULL
  pva <- NULL
  runtime <- proc.time()
  x <- as.matrix(x)
  y <- as.numeric(y)
  z <- cbind(x, y)
  all.dc <- Rfast::colrange(z, cont = FALSE)
  
  if ( is.null(univ) ) {
    a <- Rfast::g2tests(data = z, x = 1:p, y = p + 1, dc = all.dc)
    stat <- a$statistic
    pval <- pchisq(stat, a$df, lower.tail = FALSE, log.p = TRUE)
    univ$stat <- stat
    univ$pvalue <- pval
    n.tests <- p
  } else {
    stat <- univ$stat
    pval <- univ$pvalue
    n.tests <- 0
  }
 
  s <- which(pval < sig)
  
  if ( length(s) > 0 ) {
    sel <- which.min(pval)
    sela <- sel
    s <- s[ - which(s == sel) ]
    sa <- stat[sel]
    pva <- pval[sel]
    stat <- numeric(p)
    pval <- numeric(p)
    #########
    while ( sum(s>0) > 0 ) {
      for ( i in ind[s] )  {
        zz <- cbind(y, x[, c(i, sela)] )
        k <- length(sela)
        dc <- all.dc[c(p + 1, i, sela)]
        mod <- Rfast::g2Test(zz, x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
        stat[i] <- mod$statistic
        pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
      }
      n.tests <- n.tests + length( ind[s] ) 

      s <- which(pval < sig) 
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel > 0] )
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        stat <- numeric(p)
        pval <- numeric(p)
      }  
    } ## end while ( sum(s > 0) > 0 )
    
    card <- sum(sela > 0)
    
    if (K == 1) {
      stat <- numeric(p)
      pval <- numeric(p)
      for ( i in ind[-sela] )  {
        zz <- cbind(y, x[, c(i, sela)] )
        k <- length(sela)
        dc <- all.dc[c(p + 1, i, sela)]
        mod <- Rfast::g2Test(zz, x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
        stat[i] <- mod$statistic
        pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
      }
      n.tests[2] <- length( ind[-sela] )
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel > 0] )
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        stat <- numeric(p)
        pval <- numeric(p)
      }  
      while ( sum(s>0) > 0 ) {
        for ( i in ind[s] )  {
          zz <- cbind(y, x[, c(i, sela)] )
          k <- length(sela)
          dc <- all.dc[c(p + 1, i, sela)]
          mod <- Rfast::g2Test(zz, x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
          stat[i] <- mod$statistic
          pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel > 0] )
        s <- s[ - which(s == sel) ]
        if (sel > 0) {
          stat <- numeric(p)
          pval <- numeric(p)
        }  
      } ## end while ( sum(s>0) > 0 ) 
      card <- c(card, sum(sela>0) )
    }  ## end if ( K == 1 ) 
    
    if ( K > 1) {
      
      for ( i in ind[-sela] )  {
        zz <- cbind(y, x[, c(i, sela)] )
        k <- length(sela)
        dc <- all.dc[c(p + 1, i, sela)]
        mod <- Rfast::g2Test(zz, x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
        stat[i] <- mod$statistic
        pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
      }
      n.tests[2] <- length( ind[-sela] ) 
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel > 0] )
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        stat <- numeric(p)
        pval <- numeric(p)
      }  
      while ( sum(s > 0) > 0 ) {
        for ( i in ind[s] )  {
          zz <- cbind(y, x[, c(i, sela)] )
          k <- length(sela)
          dc <- all.dc[c(p + 1, i, sela)]
          mod <- Rfast::g2Test(zz, x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
          stat[i] <- mod$statistic
          pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
        }
        
        n.tests[2] <- n.tests[2] + length( ind[s] )  
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel > 0] )
        s <- s[ - which(s == sel) ]
        if (sel > 0) {
          stat <- numeric(p)
          pval <- numeric(p)
        }  
      } ## end while ( sum(s>0) > 0 ) 
      
      card <- c(card, sum(sela > 0) )
      vim <- 1
      while ( vim < K  & card[vim + 1] - card[vim] > 0 ) {
        vim <- vim + 1
        for ( i in ind[-sela] )  {
          zz <- cbind(y, x[, c(i, sela)] )
          k <- length(sela)
          dc <- all.dc[c(p + 1, i, sela)]
          mod <- Rfast::g2Test(zz, x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
          stat[i] <- mod$statistic
          pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
        }
        n.tests[vim + 1] <- length( ind[-sela] )
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel > 0] )
        s <- s[ - which(s == sel) ]
        if (sel > 0) {
          stat <- numeric(p)
          pval <- numeric(p)
        }     
        while ( sum(s > 0) > 0 ) {
          for ( i in ind[s] )  {
            zz <- cbind(y, x[, c(i, sela)] )
            k <- length(sela)
            dc <- all.dc[c(p + 1, i, sela)]
            mod <- Rfast::g2Test(zz, x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
            stat[i] <- mod$statistic
            pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
          }
          n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
          s <- which(pval < sig)
          sel <- which.min(pval) * ( length(s)>0 )
          sa <- c(sa, stat[sel]) 
          pva <- c(pva, pval[sel])
          sela <- c(sela, sel[sel > 0] )
          s <- s[ - which(s == sel) ]
          if (sel > 0) {
            stat <- numeric(p)
            pval <- numeric(p)
          }  
        } ## end while ( sum(s > 0) > 0 ) 
        card <- c(card, sum(sela>0) )
      }  ## end while ( vim < K )
    } ## end if ( K > 1)
  } ## end if ( length(s) > 0 )
  
  runtime <- proc.time() - runtime
  len <- sum( sela > 0 )
  if (len > 0) {
    res <- cbind(sela[1:len], sa[1:len], pva[1:len] )
    info <- matrix(nrow = length(card), ncol = 2)
    info[, 1] <- card
    info[, 2] <- n.tests
  } else {
    res <- matrix(c(0, 0, 0), ncol = 3)
    info <- matrix(c(0, p), ncol = 2)
  }  
  colnames(res) <- c("Vars", "stat", "log p-value")
  rownames(info) <- paste("K=", 1:length(card)- 1, sep = "")
  colnames(info) <- c("Number of vars", "Number of tests")
  result <- list(univ = univ, res = res, info = info, runtime = runtime)
  
  result$back.rem <- 0
  result$back.n.tests <- 0

  if ( backward ) {
    
    if (result$info[1, 1] > 0) {
      a <- bs.g2(y, x[, result$res[, 1], drop = FALSE], threshold = alpha)
      
      if ( typeof(a) == "list" ) {
        result$back.rem <- result$res[a$info[, 1], 1]
        back.n.tests <- sum( dim(result$res)[1] : (dim(a$mat)[1] - 1) )
        sel <- result$res[a$mat[, 1], 1] 
        stat <- a$mat[, 3]
        pval <- a$mat[, 2]
        result$res <- cbind(sel, stat, pval)
        result$back.n.tests <- back.n.tests
        result$runtime <- result$runtime + a$runtime
      } else {
        back.rem <- 0
        back.n.tests <- 0
        result$back.rem <- back.rem
        result$back.n.tests <- back.n.tests
        result$runtime <- result$runtime 
      }  
    }   ## end if (result$info[1, 1] > 0)
  }  ## end if ( backward )
  
  result
}
