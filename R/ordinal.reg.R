ordinal.reg <- function(formula, data) {
   call <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf$drop.unused.levels <- TRUE
   mf[[1]] <- as.name("model.frame")
   mf <- eval.parent(mf)
   mt <- attr(mf, "terms")
   y <- model.response(mf)
   x <- model.matrix(mt, mf, contrasts)
   k <- length( unique(y) )

  if ( k == 2 ) {
    mod <- glm(formula, data = data, binomial)
    res <- list(be = mod$coefficients, devi = mod$deviance)
    
  } else if ( sum(x) == 2 * length(y) ) {
    mod <- Rfast::ordinal.mle(y)
    res <- list(be = mod$param, devi = - 2 * mod$loglik) 
    
  } else {
    yd <- model.matrix( ~ y - 1)
    y <- as.numeric(y)
    m <- dim(x)[2]
    b <- matrix(0, m, k - 1)
    oop <- options(warn = -1) 
    on.exit( options(oop) )
    for (i in 1:(k - 1) ) {
      y1 <- y
      y1[y <= i ] <- 0
      y1[ y > i ] <- 1
      b[, i] <-  - glm.fit(x, y1, family = binomial(logit) )$coefficients 
    }
    ia <- which( is.na( rowSums(b) ) )
    if ( length(ia) > 0 ) {
      u <- x[, -ia] %*% b[-ia, ]
    } else  u <- x %*% b
    cump <- 1 / (1 + exp(-u) )
    p <- cbind(cump[, 1], Rfast::coldiffs(cump), 1 - cump[, k - 1] )
    mess <- NULL

    if ( any( p < 0) ) {
      poia <- which(p < 0, arr.ind = TRUE)[, 1]
      a <- p[poia, , drop = FALSE]
      a <- abs( a )  
      a <- a / Rfast::rowsums(a)   
      p[poia, ] <- a
      mess <- "problematic region"
    }
  
    rownames(b) <- colnames(x)
    colnames(b) <- paste("Y", 1:(k-1), sep = "" )
    devi <-  - 2 * sum( log( p[yd > 0] ) )
    res <- list(message = mess, be = b, devi = devi)
  }
  res
}
 
 
