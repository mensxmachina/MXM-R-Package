pi0est <- function(p, lambda = seq(0.05, 0.95, by = 0.01), dof = 3) {
  p <- sort(p)
  len <- length(lambda)
  prob <- numeric(len)
  for (i in 1:len)   prob[i] <- mean( p > lambda[i] ) / (1 - lambda[i])
  spi0 <- smooth.spline(lambda, prob, df = dof)
  min( predict(spi0, x = lambda[len])$y, 1 )
}

