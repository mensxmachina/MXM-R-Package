betamle.wei <- function(y, wei) {
  
  n <- length(y)
  w <- wei / sum(wei)
  ly1 <- sum( w * log(y) )      ;      ly2 <- sum( w * log(1 - y) )  
  
  betawei <- function(pa, ly1, ly2) {
    a <- exp(pa[1])     ;     b <- exp(pa[2])
    lbeta(a, b) - (a - 1) * ly1 - (b - 1) * ly2
  } 
  
  iniphi <- sum( y * (1 - y) ) / Rfast::Var(y) / n 
  a1 <- sum(y) * iniphi / n        ;        a2 <- iniphi - a1
  oop <- options(warn = -1) 
  on.exit( options(oop) )
  lik <- nlm( betawei, c( log(a1), log(a2) ), ly1 = ly1, ly2 = ly2, iterlim = 10000 )
  lik2 <- nlm( betawei, lik$estimate, ly1 = ly1, ly2 = ly2, iterlim = 10000 )  
  list(iters = lik$iterations + lik2$iterations, param = exp(lik2$estimate), loglik = -lik2$minimum )
}