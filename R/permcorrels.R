permcorrels <- function(y, x, R = 999) {
  dm <- dim(x)
  tb <- matrix(nrow = dm[2], ncol = R)
  n <- dm[1]
  my <- sum(y)        ;     my2 <- sum(y^2)
  mx <- Rfast::colsums(x)     ;     mx2 <- Rfast::colsums(x^2)
  up <-  my * mx / n
  down <- sqrt( (my2 - my^2 / n) * (mx2 - mx^2 / n) )
  r <- ( Rfast::colsums(y * x) - up) / down
  test <- abs( log( (1 + r) / (1 - r) ) )  ## the test statistic
  for (i in 1:R) {
    yb <- sample(y, n)
    tb[, i] <- Rfast::colsums(yb * x)
  }  
  tb <- (tb - up) / down
  tb <- log( (1 + tb) / (1 - tb) )  ## the test statistic
  pvalue <- ( Rfast::rowsums( abs(tb) > test ) + 1 ) / (R + 1)  ## bootstrap p-value
  res <- cbind(r, test, log( pvalue) )
  colnames(res) <- c("correlation", "statistic", "log p-value")
  res
}
