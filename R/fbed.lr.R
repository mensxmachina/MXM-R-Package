fbed.lr <- function(y, x, alpha = 0.05, univ = NULL, test = NULL, wei = NULL, K = 0) { 

  dm <- dim(x)
  p <- dm[2]
  ind <- 1:p
  sig <- log(alpha)
  sela <- NULL
  card <- 0
  sa <- NULL
  pva <- NULL
  
  test <- test.maker(test)
  
  runtime <- proc.time()
  if ( is.null(univ) ) {
    univ <- MXM::univregs(y, x, test = test)
    stat <- univ$stat
    pval <- univ$pvalue
    n.tests <- p
    zevar <- which(stat == 0)
  } else {
    stat <- univ$stat
    pval <- univ$pvalue
    zevar <- which(stat == 0)
    n.tests <- 0
  }  
  s <- which(pval < sig)
  
  if ( length(s) > 0 ) {
    if ( identical(test, testIndSPML) ) {
      if ( !is.matrix(y) )   y <- cbind( cos(y), sin(y) )
    }  
      
    sel <- which.min(pval)
    sela <- sel
    s <- s[ - which(s == sel) ]
    sa <- stat[sel]
    pva <- pval[sel]
    #########
    while ( sum(s>0) > 0 ) {
      mod <- MXM::cond.regs(y, x, xIndex = ind[s], csIndex = sela, test = test, wei = wei)
      n.tests <- n.tests + length( ind[s] ) 
      stat <- mod$stat
      pval <- mod$pvalue
      s <- which(pval < sig) 
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
    } ## end while ( sum(s > 0) > 0 )
    
    card <- sum(sela > 0)
    
    if (K == 1) {
      mod <- MXM::cond.regs(y, x, xIndex = ind[-c(sela, zevar)], csIndex = sela, test = test, wei = wei)
      n.tests[2] <- length( ind[-c(sela, zevar)] )
      stat <- mod$stat
      pval <- mod$pvalue
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      while ( sum(s>0) > 0 ) {
        mod <- MXM::cond.regs(y, x, xIndex = ind[s], csIndex = sela, test = test, wei = wei)
        n.tests[2] <- n.tests[2] + length( ind[s] )
        stat <- mod$stat
        pval <- mod$pvalue
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
      } ## end while ( sum(s>0) > 0 ) 
      card <- c(card, sum(sela>0) )
    }  ## end if ( K == 1 ) 
    
    if ( K > 1) {
      
      mod <- MXM::cond.regs(y, x, xIndex = ind[-c(sela, zevar)], csIndex = sela, test = test, wei = wei)
      n.tests[2] <- length( ind[-c(sela, zevar)] ) 
      stat <- mod$stat
      pval <- mod$pvalue
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
 
      while ( sum(s > 0) > 0 ) {
        mod <- MXM::cond.regs(y, x, xIndex = ind[s], csIndex = sela, test = test, wei = wei)
        n.tests[2] <- n.tests[2] + length( ind[s] )  
        stat <- mod$stat
        pval <- mod$pvalue
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
      } ## end while ( sum(s>0) > 0 ) 
      
      card <- c(card, sum(sela > 0) )
      vim <- 1
      while ( vim < K  & card[vim + 1] - card[vim] > 0 ) {
        vim <- vim + 1
        mod <- MXM::cond.regs(y, x, xIndex = ind[-c(sela, zevar)], csIndex = sela, test = test, wei = wei)
        n.tests[vim + 1] <- length( ind[-c(sela, zevar)] ) 
        stat <- mod$stat
        pval <- mod$pvalue
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
 
        while ( sum(s > 0) > 0 ) {
          mod <- MXM::cond.regs(y, x, xIndex = ind[s], csIndex = sela, test = test, wei = wei)
          n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
          stat <- mod$stat
          pval <- mod$pvalue
          s <- which(pval < sig)
          sel <- which.min(pval) * ( length(s)>0 )
          sa <- c(sa, stat[sel]) 
          pva <- c(pva, pval[sel])
          sela <- c(sela, sel[sel>0])
          s <- s[ - which(s == sel) ]
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
