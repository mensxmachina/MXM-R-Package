cor.drop1 <- function(y, x, logged = FALSE) {
  if ( is.matrix(x) )    x <- data.frame(x)
  mod <- lm(y ~., data = x)
  a <- 1:dim(x)[2]
  bna <- which( is.na(mod$coefficients[-1]) )
  if ( length(bna) > 0 ) {
    a[bna] <- 0
    a[-bna] <- summary(mod)[[4]][-1, 3]^2
  } else  a <- summary(mod)[[4]][-1, 3]^2
  dof <- summary(mod)[[10]][3]
  r <- sqrt( a / (a + dof) )
  stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) ) * sqrt(dof - 1) 
  if ( logged ) {
    pval <- log(2) + pt(stat, dof - 1, lower.tail = FALSE, log.p = TRUE)
  } else pval <- 2 * pt(stat, dof - 1, lower.tail = FALSE)
  cbind(stat, pval)
}

