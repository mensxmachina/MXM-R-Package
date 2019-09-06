supervised.pca <- function(target, dataset, indices, center = TRUE, scale = TRUE, colours = NULL, graph = TRUE) {
  
  mod.all <- prcomp(dataset, center = center, scale = scale)
  mod.sel <- prcomp(dataset[, indices], center = center, scale = scale)
  rat <- sum(mod.sel$sdev^2) / sum(mod.all$sdev^2)
  
  if ( graph ) {
    if ( is.null(colours) )  target <- as.numeric( as.factor(target) )
    plot( mod.all$x[, 1:2], col = target, main = "Scores using all variables" )
    dev.new()
    plot( mod.sel$x[, 1:2], col = target, main = "Scores using the selected variables" )
  }
  
  list(mod.all = mod.all, mode.sel = mod.sel, var.percent = rat)
}