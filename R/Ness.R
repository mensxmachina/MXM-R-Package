Ness <- function(propNt, N, K = 10000) {
  
  stats <- numeric(K)
  dm <- length(propNt)
  f <- round(N * propNt)
  if ( sum(f) < N)   f[1] <- f[1] + 1  
  if ( sum(f) > N)   f[1] <- f[1] + 1  
  
  Nt <- rep(0:(dm - 1), f)
  dat <- matrix(0, nrow = N, ncol = 3)
  dat[, 3] <- Nt
  
  for (i in 1:K) {
    dat[, 1:2] <- matrix( sample(0:1, 2 * N, replace = TRUE), ncol = 2 )
    stats[i] <- Rfast::g2Test(data = dat, 1, 2, 3,  dc = c(2, 2, dm) )$statistic
  }
  
  stats <- 0.5 * stats / N
  stats <- sort(stats)
  statsdot <- qchisq( seq(0, K - 1)/ (K - 1) , 1)
  round( median(statsdot/stats) * 0.5 )
}