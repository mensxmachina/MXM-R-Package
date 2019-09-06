zipmle.wei <- function(y, wei) {
  n <- length(y)
  poia <- which(y == 0)
  n0 <- length(poia)   ;    n1 <-  n - n0
  y1 <- y[ -poia ]    
  wei <- wei / sum(wei)
  w0 <- wei[poia]     ;    w1 <- wei[-poia]
  pini <- n0 / n1
  pini <- log(pini / (1 - pini) )
  lik <- nlm( zipwei, c(pini, log( sum(y1) / n ) ), y1 = y1, w1 = w1, w0 = w0, iterlim = 10000 )
  lik2 <- nlm( zipwei, lik$estimate, y1 = y1, w1 = w1, w0 = w0, iterlim = 10000 )  
  prop <- exp(lik2$estimate[1]) / ( 1 + exp(lik2$estimate[1]) )
  list(iters = lik$iterations + lik2$iterations, prop = prop, lam = exp(lik2$estimate[2]),
       loglik = -lik2$minimum - sum( w1 * lgamma(y1 + 1) ) )
}
