auc <- function(group, preds, roc = FALSE, cutoffs = NULL) {
  
  ri <- rank(preds)
  n <- length(preds)
  n1 <- sum( group )
  n0 <- n - n1
  s1 <- sum( ri[ group == 1 ] )
  auc <- ( s1 - 0.5 * n1 * (n1 + 1) ) / n0 / n1
  result <- list(auc = auc)
  
  if ( roc ) {
    if ( is.null(cutoffs) ) {
      boun <- quantile(preds, probs = c(0.01, 0.99))
      lena <- seq(from = boun[2], to = boun[1], length = 100)
    } else  lena <- cutoffs  
    nu <- length(lena)
    sens <- spec <- numeric(nu)
    npos <- n1
    nnegs <- n0
    for ( i in 1:nu ) {
      estp <- as.numeric( preds >= lena[i] )      
      sens[i] <-  sum( group == 1  &  estp == 1 ) / npos
      spec[i] <-  sum( group == 0  &  estp == 0 ) / nnegs
    }
    
    spec <- c(spec, 0)
    sens<- c(0, sens)
    plot( 1 - spec, sens, pch = 16, xlab = "1- specificity", ylab = "Sensitivity", col = 3, 
          xlim = c(0, 1), ylim = c(0, 1), type = "l")
    abline(a = 0, b = 1, lty = 2, col = 2)
    qa <- which.max( sens + spec - 1 )
    points( 1 - spec[qa], sens[qa], lwd = 2, col = 4, pch = 3)
    youden <- c( 1 - spec[qa], sens[qa], spec[qa] + sens[qa] - 1 )
    names(youden) <- c("1-specificity", "sensitivity", "Youden's J")
    result <- list(cutoffs = cutoffs, sensitivity = sens, specificity = spec, youden = youden, auc = auc)
  }
  
  result
}