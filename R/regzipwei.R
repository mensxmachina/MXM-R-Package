regzipwei <- function(pa, y1, x, n1, w1, w0, poia) {
  a <- exp(pa[1]) / ( 1 + exp(pa[1]) )    ;    b <- pa[-1]
  es <- tcrossprod(b, x)
  - sum( w0 * log( a + (1 - a) * exp( - exp( es[poia]) ) ) ) - sum(w1 * log( 1 - a ) ) - sum( w1 * y1 * es[ -poia ] - w1 * exp(es[ -poia ]) ) 
}  
