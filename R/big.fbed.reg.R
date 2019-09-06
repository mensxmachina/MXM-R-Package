big.fbed.reg <- function(target = NULL, dataset, threshold = 0.01, ini = NULL, test = "testIndLogistic", K = 0, backward = FALSE) { 

  runtime <- proc.time()
    
  if ( is.null(target) ) {
    y <- dataset[, 1]
    dataset <- bigmemory::sub.big.matrix(dataset, firstCol = 2)
  } else   y <- target 
  
  dm <- dim(dataset)
  n <- dm[1]     ;     p <- dm[2]
  ind <- 1:p
  sig <- log(threshold)
  sela <- NULL
  card <- 0
  sa <- NULL
  pva <- NULL
  
  if ( is.null(ini) ) {
    
    if ( test == "testIndLogistic"  |  test == "testIndPois"  |  test == "testIndMultinom"  |  test == "censIndWR" ) {
      if (test == "testIndLogistic") {
        mod <- big.score.univregs(target = y, dataset = dataset, test = testIndLogistic)
      } else if (test == "testIndPois") {
        mod <- big.score.univregs(target = y, dataset = dataset, test = testIndPois)
      } else if (test == "testIndMultinom") {
        mod <- big.score.univregs(target = y, dataset = dataset, test = testIndMultinom)
      } else if (test == "censIndWR") {
        mod <- big.score.univregs(target = y, dataset = dataset, test = censIndWR) 
      }
      stat <- mod$stat
      pval <- mod$pvalue
      
    } else {
      if ( test == "testIndFisher" )   {
        r <- as.vector( cor(y, dataset[]) )
        stat <- 0.5 * log( (1 + r)/(1 - r) ) * sqrt(n - 3)
        pval <- log(2) + pt(abs(stat), n - 3, lower.tail = FALSE, log.p = TRUE)
      } else if (test == "testIndQPois") {
        stat <- numeric(p)
        phi <- rep(1, p)
        m <- colSums(dataset[])
        x2 <- colSums(dataset[]^2)
        s <- (x2 - m^2/n)
        ind[s == 0] <- 0
        for (i in ind)  {
          mod <- Rfast::qpois.reg(dataset[, i], y, maxiters = 1000)
          stat[i] <- mod$devi
          phi[i] <- mod$phi
        }  
        stat <- ( 2 * sum( y * log(n * y/sum(y) ), na.rm = TRUE) - stat )/phi
        pval <- pf(stat, 1, n - 2, lower.tail = FALSE, log.p = TRUE)        
      } 
       # else if ( test == "testIndQBinom" ) {
       # stat <- numeric(p)
      #  for (i in 1:p)    stat[i] <- Rfast::prop.reg(y, x[, i], varb = "glm")$info[2, 3]
      #  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
    }
    ini$stat <- stat
    ini$pvalue <- pval
    n.tests <- p
  } else {
    stat <- ini$stat
    pval <- ini$pvalue
    n.tests <- 0
  }  
  
  if ( test == "censIndWR" )   y <- survival::Surv(y, rep(1, n) )  
  
  s <- which(pval < sig)
  
  if ( length(s) > 0 ) {
    sel <- which.min(pval)
    sela <- sel
    s <- s[ - which(s == sel) ]
    lik1 <- big.model(y, dataset[, sel], test)$devi
    sa <- stat[sel]
    pva <- pval[sel]
    lik2 <- rep( lik1, p )
    phi <- numeric(p)
    #########
    while ( sum(s>0) > 0 ) {
      for ( i in ind[s] )  {
        fit2 <- big.model( y, dataset[, c(sela, i)], test )
        lik2[i] <- fit2$devi
        phi[i] <- fit2$phi
      }    
      if ( any( is.na(phi) ) ) {
        stat <- abs(lik2 - lik1)
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      } else {
        stat <- abs(lik2 - lik1)/phi
        pval <- pf(stat, 1, n - length(sela) - 2, lower.tail = FALSE, log.p = TRUE)
      }  
      n.tests <- n.tests + length( ind[s] ) 
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
      for ( i in ind[-sela] )  {
        fit2 <- big.model( y, dataset[, c(sela, i)], test )
        lik2[i] <- fit2$devi
        phi[i] <- fit2$phi
      }    
      if ( any( is.na(phi) ) ) {
        stat <- abs(lik2 - lik1)
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      } else {
        stat <- abs(lik2 - lik1)/phi
        pval <- pf(stat, 1, n - length(sela) - 2, lower.tail = FALSE, log.p = TRUE)
      }  
      n.tests[2] <- length( ind[-sela] )
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
          fit2 <- big.model( y, dataset[, c(sela, i)], test )
          lik2[i] <- fit2$devi
          phi[i] <- fit2$phi
        }    
        if ( any( is.na(phi) ) ) {
          stat <- abs(lik2 - lik1)
          pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        } else {
          stat <- abs(lik2 - lik1)/phi
          pval <- pf(stat, 1, n - length(sela) - 2, lower.tail = FALSE, log.p = TRUE)
        }  
        n.tests[2] <- n.tests[2] + length( ind[s] )
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
      
      for ( i in ind[-sela] )  {
        fit2 <- big.model( y, dataset[, c(sela, i)], test )
        lik2[i] <- fit2$devi
        phi[i] <- fit2$phi
      }    
      if ( any( is.na(phi) ) ) {
        stat <- abs(lik2 - lik1)
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      } else {
        stat <- abs(lik2 - lik1)/phi
        pval <- pf(stat, 1, n - length(sela) - 2, lower.tail = FALSE, log.p = TRUE)
      }  
      n.tests[2] <- length( ind[-sela] ) 
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
        for ( i in ind[s] ) {  
          fit2 <- big.model( y, dataset[, c(sela, i)], test )
          lik2[i] <- fit2$devi
          phi[i] <- fit2$phi
        }    
        if ( any( is.na(phi) ) ) {
          stat <- abs(lik2 - lik1)
          pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        } else {
          stat <- abs(lik2 - lik1)/phi
          pval <- pf(stat, 1, n - length(sela) - 2, lower.tail = FALSE, log.p = TRUE)
        } 
        n.tests[2] <- n.tests[2] + length( ind[s] )  
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
        for ( i in ind[-sela] ) {
          fit2 <- big.model( y, dataset[, c(sela, i)], test )
          lik2[i] <- fit2$devi
          phi[i] <- fit2$phi
        }    
        if ( any( is.na(phi) ) ) {
          stat <- abs(lik2 - lik1)
          pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        } else {
          stat <- abs(lik2 - lik1)/phi
          pval <- pf(stat, 1, n - length(sela) - 2, lower.tail = FALSE, log.p = TRUE)
        } 
        n.tests[vim + 1] <- length( ind[-sela] )
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
          for ( i in ind[s] ) {  
            fit2 <- big.model( y, dataset[, c(sela, i)], test )
            lik2[i] <- fit2$devi
            phi[i] <- fit2$phi
          }    
          if ( any( is.na(phi) ) ) {
            stat <- abs(lik2 - lik1)
            pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
          } else {
            stat <- abs(lik2 - lik1)/phi
            pval <- pf(stat, 1, n - length(sela) - 2, lower.tail = FALSE, log.p = TRUE)
          } 
          n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
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
  result <- list(univ = ini, res = res, info = info, runtime = runtime)
  
  result$back.rem <- 0
  result$back.n.tests <- 0
  
  if ( backward ) {
    
    if (result$info[1, 1] > 0) {
      a <- bsreg.big(target, dataset[, result$res[, 1], drop = FALSE], threshold = threshold, test = test)
      
      if ( typeof(a) == "list" ) {
        result$back.rem <- result$res[a$info[, 1], 1]
        back.n.tests <- sum( dim(result$res)[1] : dim(a$mat)[1] )
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
