################################
#### Permutation based hypothesis testing 
#### for a zero correlation coefficient 
################################
permcor <- function(x1, x2, R = 999) {
  Rfast::permcor(x1, x2, R = R)
}


# permcor <- function(x1, x2, R = 999) {
  # x is a 2 column matrix containing the data
  # type can be either "pearson" or "spearman"
  # R is the number of permutations
  # n <- length(x1)
  # m1 <- sum(x1)     ;     m12 <- sum(x1^2)
  # m2 <- sum(x2)     ;     m22 <- sum(x2^2)
  # up <-  m1 * m2 / n
  # down <- sqrt( (m12 - m1^2 / n) * (m22 - m2^2 / n) )
  # r <- ( sum(x1 * x2) - up) / down
  # test <- log( (1 + r) / (1 - r) )  ## the test statistic
  # sxy <- numeric(R)
  # B <- round( sqrt(R) ) 
  # xp <- matrix(0, n, B)
  # yp <- matrix(0, n, B)
  # for (i in 1:B) {
    # xp[, i] <- sample(x1, n)
    # yp[, i] <- sample(x2, n) 
  # }
  # sxy <- crossprod(xp, yp)
  # rb <- (sxy - up) / down
  # tb <- log( (1 + rb) / (1 - rb) )  ## the test statistic
  # pvalue <- ( sum( abs(tb) > abs(test) ) + 1 ) / (B^2 + 1)  ## bootstrap p-value
  # res <- c( r, pvalue )
  # names(res) <- c('correlation', 'p-value')
  # res
# }