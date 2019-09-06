bn.skel.utils2 <- function(mod, G = NULL, roc = TRUE, alpha = 0.01) {
  
  ell <- function(a, prob, p)  sum( log(prob + (1 - prob) * a * p^(a-1) ) )
  preds <- mod$pvalue
  theseis <- which(row(preds) < col(preds), arr.ind = TRUE)
  theseis <- cbind(theseis, preds[theseis] )
  p <- exp( theseis[, 3] )
  m <- length(p)
  pi0 <- pi0est(p)
  if ( pi0 < 1 ) {
    ksi <- optimise(ell, prob = pi0, p = p, c(0, 1), maximum = TRUE)$maximum
    po <- pi0/( ksi * p^(ksi - 1) * (1 - pi0) )
    pxy <- 1 / (1 + po)
  } else pxy <- numeric(m)
  pxy <- cbind(theseis, pxy )
  
  pxy <- pxy[order( - pxy[, 4]), ]
  colnames(pxy)[3:4] <- c( "p-value", "confidence")
  pxy[, 3] <- exp(pxy[, 3]) 
  area <- NULL
  if ( !is.null(G) ) {
    group <- G[ pxy[, 1:2] ] 
    area <- auc(group,  pxy[, 4], roc = roc)
  }
  lower <- conf.edge.lower(pxy[, 4])
  list(area = area, pxy = pxy, lower = lower)
}