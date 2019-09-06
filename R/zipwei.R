zipwei <- function(pa, y1, w1, w0) {
  a <- exp(pa[1]) / ( 1 + exp(pa[1]) )      
  lam <- exp(pa[2])
  - sum(w0) * log( a + (1 - a) * exp( - lam ) ) - sum(w1) * log( 1 - a ) - sum(w1 * y1) * pa[2] + sum(w1) * lam 
}  
