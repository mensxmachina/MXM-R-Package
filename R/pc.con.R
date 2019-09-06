#######################
#######################
##### Fast skeleton of the PC algorithm for continuous data only using pearson or spearman
#####
#######################
#######################
pc.con <- function(dataset, method = "pearson", alpha = 0.01) {
  ## dataset contains the data, it must be a matrix 
  title <- deparse( substitute(dataset) )
  res <- Rfast::pc.skel(dataset = dataset, method = method, alpha = alpha)
  info <- summary( Rfast::rowsums(res$G) )
  n <- dim(dataset)[2]
  density <- sum(res$G) / n / ( n - 1 ) 
  res$density <- density
  res$info <- info
  res$title <- title 
  res  
}