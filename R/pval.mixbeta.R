pval.mixbeta <- function(p) {
  
  ell <- function(a, prob, p)  sum( log(prob + (1 - prob) * a * p^(a-1) ) )
  pi0 <- pi0est(p)
  if ( pi0 < 1 ) {
    ksi <- optimise(ell, prob = pi0, p = p, c(0, 1), maximum = TRUE)$maximum
  } else ksi <- NULL
  res <- c(pi0, ksi)
  names(res) <- c("pi0", "ksi")
  res  
}