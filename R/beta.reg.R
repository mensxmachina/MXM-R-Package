beta.reg <- function(target, dataset, wei = NULL) {
  
  n <- length(target)
  if ( NCOL(dataset) == 0 ) {
    if ( is.null(wei) ) {
      mod <- Rfast::beta.mle(target)
    } else mod <- betamle.wei(target, wei)
    res <- list(phi = sum(mod$param), be = mod$param[1]/mod$param[2], loglik = mod$loglik)
    
  } else {
    x <- model.matrix(target ~ ., data.frame(dataset) )
    iniphi <- log( sum( target * (1 - target) ) / Rfast::Var(target) / n )
    oop <- options(warn = -1) 
    on.exit( options(oop) )
  
    if ( is.null(wei) ) {
      ly1 <- log(1 - target)     ;    ly <- log(target) - ly1           
      sly1 <- sum(ly1)      ;    sly2 <- sum( log(target) ) + sly1   
      mod1 <- nlm(regbeta, c( iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
      mod2 <- nlm(regbeta, mod1$estimate, ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
    } else {
      w <- wei / sum(wei)
      ly1 <- w * log(1 - target)     ;    ly <- w * log(target) - ly1
      sly1 <- sum(ly1)    ;    sly2 <- sum( w * log( target ) ) + sly1   
      mod1 <- nlm(regbetawei, c(iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
      mod2 <- nlm(regbetawei, mod1$estimate, ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
    }
    res <-list(phi = exp(mod2$estimate[1]), be = mod2$estimate[-1], loglik = - mod2$minimum - sly2)
  }
  res
}