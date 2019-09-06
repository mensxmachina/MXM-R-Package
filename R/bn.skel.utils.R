bn.skel.utils <- function(mod, G = NULL, roc = TRUE, alpha = 0.01) {
  area <- NULL
  preds <- mod$pvalue
  preds <- preds[ upper.tri(preds) ]
  if ( !is.null(G) ) {
    group <- G[ upper.tri(G) ] 
    area <- auc(group,  -preds, roc = roc)
  }  
  p <- exp(preds)
  sig.p <- p[ which(p <= alpha) ]
  fdr <- min(1, length(p) / length(sig.p) * max(sig.p) )
  theseis <- which( mod$pvalue < log(alpha), arr.ind = TRUE )
  sig.p <- cbind(theseis, mod$pvalue[theseis])
  sig.p <- Rfast::rowSort(sig.p)
  sig.p <- unique(sig.p)
  sig.p <- cbind(sig.p[, 2:3], sig.p[, 1] )
  sig.p <- sig.p[order(sig.p[, 3]), ]
  sig.p[, 3] <- exp(sig.p[, 3])
  theseis <- NULL
  list(area = area, fdr = fdr, sig.pvalues = sig.p)
}