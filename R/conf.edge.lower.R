conf.edge.lower <- function(p) {
  
  p <- sort(p)
  fn <- c(0, cumsum(p <= p) / length(p) )
  dp <- diff( c(0, sort(p), 1) )
  scou <- function(tau)  sum( abs(fn - tau) * dp )
  tau <- optimise(scou, c(0, 1) )$minimum
  p[ which(p>tau)[1] ]
}

  
  

