fbed.geelm <- function(y, x, id, univ = NULL, alpha = 0.05, wei = NULL, K = 0, correl = "exchangeable", se = "jack") { 
  
  dm <- dim(x)
  p <- dm[2]
  n <- dm[1]
  ind <- 1:p
  sig <- log(alpha)
  stat <- numeric(p)
  sela <- NULL
  card <- 0
  sa <- NULL
  pva <- NULL
  
  zevar <- Rfast::check_data(x)
  if ( sum( zevar > 0 ) > 0 )  x[, zevar] <- rnorm( n * length(zevar) )

  runtime <- proc.time()
  
  if ( is.null(univ) ) {
    for ( i in ind ) {
      fit2 <- geepack::geeglm( y ~ x[, i], family = gaussian, id = id, weights = wei, corstr = correl, std.err = se )
      stat[i] <- summary(fit2)[[ 6 ]][2, 3]
    }  
    n.tests <- p
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
  }  
  s <- which(pval < sig)
  
  if ( length(s) > 0 ) {
    sel <- which.min(pval)
    sela <- sel
    s <- s[ - which(s == sel) ]
    sa <- stat[sel]
    pva <- pval[sel]
    stat <- numeric(p )
    #########
    while ( sum(s>0) > 0 ) {
       
      for ( i in ind[s] )  {
        fit2 <- try( geepack::geeglm( y ~ x[, sela] + x[, i], family = gaussian, id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
        if ( identical( class(fit2), "try-error" ) ) {
          stat[i] <- 0
        } else {
          mod <- summary(fit2)[[ 6 ]]
          nr <- dim(mod)[1]
          stat[i] <- mod[nr, 3]
        }  
      } 
      n.tests <- n.tests + length( ind[s] ) 
      pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig) 
      sel <- which.min(pval) * ( length(s) > 0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0)  stat <- rep(0, p)
    } ## end while ( sum(s > 0) > 0 )
    
    card <- sum(sela > 0)
    
    if (K == 1) {
      for ( i in ind[-c(sela, zevar)] )  {
        fit2 <- try( geepack::geeglm( y ~ x[, sela] + x[, i], family = gaussian, id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
        if ( identical( class(fit2), "try-error" ) ) {
          stat[i] <- 0
        } else {
          mod <- summary(fit2)[[ 6 ]]
          nr <- dim(mod)[1]
          stat[i] <- mod[nr, 3]
        } 
      }
      n.tests[2] <- length( ind[-c(sela, zevar)] )
      pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0)   stat <- numeric(p)
      while ( sum(s>0) > 0 ) {
        for ( i in ind[s] )  {
          fit2 <- try( geepack::geeglm( y ~ x[, sela] + x[, i], family = gaussian, id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
          if ( identical( class(fit2), "try-error" ) ) {
            stat[i] <- 0
          } else {
            mod <- summary(fit2)[[ 6 ]]
            nr <- dim(mod)[1]
            stat[i] <- mod[nr, 3]
          } 
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0)  stat <- numeric(p)
      } ## end while ( sum(s>0) > 0 ) 
      card <- c(card, sum(sela>0) )
    }  ## end if ( K == 1 ) 
    
    if ( K > 1) {
      
      for ( i in ind[-c(sela, zevar)] )  {
        fit2 <- try( geepack::geeglm( y ~ x[, sela] + x[, i], family = gaussian, id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
        if ( identical( class(fit2), "try-error" ) ) {
          stat[i] <- 0
        } else {
          mod <- summary(fit2)[[ 6 ]]
          nr <- dim(mod)[1]
          stat[i] <- mod[nr, 3]
        } 
      }
      n.tests[2] <- length( ind[-c(sela, zevar)] )
      pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0)  stat <- numeric(p)
      
      while ( sum(s > 0) > 0 ) {
        for ( i in ind[s] )  {
          fit2 <- try( geepack::geeglm( y ~ x[, sela] + x[, i], family = gaussian, id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
          if ( identical( class(fit2), "try-error" ) ) {
            stat[i] <- 0
          } else {
            mod <- summary(fit2)[[ 6 ]]
            nr <- dim(mod)[1]
            stat[i] <- mod[nr, 3]
          } 
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )  
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0)  stat <- numeric(p)
      } ## end while ( sum(s>0) > 0 ) 
      
      card <- c(card, sum(sela > 0) )
      vim <- 1
      while ( vim < K  & card[vim + 1] - card[vim] > 0 ) {
        vim <- vim + 1
        for ( i in ind[-c(sela, zevar)] )  {
          fit2 <- try( geepack::geeglm( y ~ x[, sela] + x[, i], family = gaussian, id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
          if ( identical( class(fit2), "try-error" ) ) {
            stat[i] <- 0
          } else {
            mod <- summary(fit2)[[ 6 ]]
            nr <- dim(mod)[1]
            stat[i] <- mod[nr, 3]
          } 
        }
        n.tests[vim + 1] <- length( ind[-c(sela,zevar)] )
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0)  stat <- numeric(p)
        while ( sum(s > 0) > 0 ) {
          for ( i in ind[s] )  {
            fit2 <- try( geepack::geeglm( y ~ x[, sela] + x[, i], family = gaussian, id = id, weights = wei, corstr = correl, std.err = se ), silent = TRUE)
            if ( identical( class(fit2), "try-error" ) ) {
              stat[i] <- 0
            } else {
              mod <- summary(fit2)[[ 6 ]]
              nr <- dim(mod)[1]
              stat[i] <- mod[nr, 3]
            } 
          }
          n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
          pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
          s <- which(pval < sig)
          sel <- which.min(pval) * ( length(s)>0 )
          sa <- c(sa, stat[sel]) 
          pva <- c(pva, pval[sel])
          sela <- c(sela, sel[sel>0])
          s <- s[ - which(s == sel) ]
          if (sel > 0)  stat <- numeric(p)
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
