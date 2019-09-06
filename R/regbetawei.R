regbetawei <- function(pa, ly, sly1, x, w) {
  phi <- exp(pa[1])    ;    b <- pa[-1]
  m <- exp( tcrossprod(b, x) )
  m <- m / ( 1 + m )
  a1 <- m * phi   ;   a2 <- phi - a1
  - lgamma(phi) + sum( w * lgamma(a1) ) + sum( w * lgamma(a2) ) - sum(a1 * ly) - phi * sly1
} 