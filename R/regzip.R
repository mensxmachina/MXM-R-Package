regzip <- function(pa, y1, x, n1, poia) {
  a <- exp(pa[1]) / ( 1 + exp(pa[1]) )    ;    b <- pa[-1]
  es <- tcrossprod(b, x)
  - sum( log( a + (1 - a) * exp( - exp( es[poia]) ) ) ) - n1 * log( 1 - a ) - sum( y1 * es[ -poia ] - exp(es[ -poia ]) ) 
}  