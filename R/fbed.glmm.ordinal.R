fbed.glmm.ordinal <- function(y, x, id, univ = NULL, alpha = 0.05, wei = NULL, K = 0) { 

  dm <- dim(x)
  p <- dm[2]
  n <- dm[1]
  ind <- 1:p
  sig <- log(alpha)
  lik1 <- logLik( ordinal::clmm(y ~ 1 + (1|id), weights = wei) )
  lik2 <- numeric(p)
  sela <- NULL
  card <- 0
  sa <- NULL
  pva <- NULL
  
  zevar <- Rfast::check_data(x)
  if ( sum( zevar > 0 ) > 0 )  x[, zevar] <- rnorm( n * length(zevar) )
  
  runtime <- proc.time()
  
  if ( is.null(univ) ) {
    for ( i in ind ) {
      fit2 <- ordinal::clmm( y ~ x[, i] + (1|id), weights = wei )
      lik2[i] <- logLik( fit2 )
    }
    n.tests <- p
    stat <- 2 * (lik2 - lik1)
    pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
    univ <- list()
    univ$stat <- stat
    univ$pvalue <- pval
	univ$stat[zevar] <- 0
	univ$pvalue[zevar] <- 0
	
  } else {
    n.tests <- 0
    stat <- univ$stat
    pval <- univ$pvalue
    lik2 <- 0.5 * stat + lik1
  }  
  s <- which(pval < sig)
  
  if ( length(s) > 0 ) {
    sel <- which.min(pval)
    sela <- sel
    s <- s[ - which(s == sel) ]
    lik1 <- lik2[sel] 
    sa <- stat[sel]
    pva <- pval[sel]
    lik2 <- rep( lik1, p )
    #########
    while ( sum(s>0) > 0 ) {
      for ( i in ind[s] )  {
        fit2 <- ordinal::clmm( y ~ x[, c(sela, i)] + (1|id), weights = wei )
        lik2[i] <- logLik( fit2 )
      }
      n.tests <- n.tests + length( ind[s] ) 
      stat <- 2 * (lik2 - lik1)
      pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig) 
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        lik1 <- lik2[sel] 
        lik2 <- rep(lik1, p)
      } 
    } ## end while ( sum(s > 0) > 0 )
    
    card <- sum(sela > 0)
    
    if (K == 1) {
      for ( i in ind[-c(sela, zevar)] )  {
        fit2 <- ordinal::clmm( y ~ x[, c(sela, i)] + (1|id), weights = wei )
        lik2[i] <- logLik( fit2 )
      }
      n.tests[2] <- length( ind[-c(sela, zevar)] )
      stat <- 2 * (lik2 - lik1)
      pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        lik1 <- lik2[sel] 
        lik2 <- rep(lik1, p)
      } 
      while ( sum(s>0) > 0 ) {
        for ( i in ind[s] )  {
          fit2 <- ordinal::clmm( y ~ x[, c(sela, i)] + (1|id), weights = wei )
          lik2[i] <- logLik( fit2 )
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )
        stat <- 2 * (lik2 - lik1)
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0) {
          lik1 <- lik2[sel] 
          lik2 <- rep(lik1, p)
        } 
      } ## end while ( sum(s>0) > 0 ) 
      card <- c(card, sum(sela>0) )
    }  ## end if ( K == 1 ) 
    
    if ( K > 1) {
      
      for ( i in ind[-c(sela, zevar)] )  {
        fit2 <- ordinal::clmm( y ~ x[, c(sela, i)] + (1|id), weights = wei )
        lik2[i] <- logLik( fit2 )
      }
      n.tests[2] <- length( ind[-c(sela, zevar)] ) 
      stat <- 2 * (lik2 - lik1)
      pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        lik1 <- lik2[sel] 
        lik2 <- rep(lik1, p)
      } 
      while ( sum(s > 0) > 0 ) {
        for ( i in ind[s] )  {
          fit2 <- ordinal::clmm( y ~ x[, c(sela, i)] + (1|id), weights = wei )
          lik2[i] <- logLik( fit2 )
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )  
        stat <- 2 * (lik2 - lik1)
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0) {
          lik1 <- lik2[sel] 
          lik2 <- rep(lik1, p)
        } 
      } ## end while ( sum(s>0) > 0 ) 
      
      card <- c(card, sum(sela > 0) )
      vim <- 1
      while ( vim < K  & card[vim + 1] - card[vim] > 0 ) {
        vim <- vim + 1
        for ( i in ind[-sela] )  {
          fit2 <- ordinal::clmm( y ~ x[, c(sela, i)] + (1|id), weights = wei )
          lik2[i] <- logLik( fit2 )
        }
        n.tests[vim + 1] <- length( ind[-sela] )
        stat <- 2 * (lik2 - lik1)
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0) {
          lik1 <- lik2[sel] 
          lik2 <- rep(lik1, p)
        }    
        while ( sum(s > 0) > 0 ) {
          for ( i in ind[s] )  {
            fit2 <- ordinal::clmm( y ~ x[, c(sela, i)] + (1|id), weights = wei )
            lik2[i] <- logLik( fit2 )
          }
          n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
          stat <- 2 * (lik2 - lik1)
          pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
          s <- which(pval < sig)
          sel <- which.min(pval) * ( length(s)>0 )
          sa <- c(sa, stat[sel]) 
          pva <- c(pva, pval[sel])
          sela <- c(sela, sel[sel>0])
          s <- s[ - which(s == sel) ]
          if (sel > 0) {
            lik1 <- lik2[sel] 
            lik2 <- rep(lik1, p)
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
  list(univ = univ, res = res, info = info, runtime = runtime)
}
