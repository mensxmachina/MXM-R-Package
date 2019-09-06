## finds descenadants of some given node(s)
findDescendants <- function(G, node = NULL, graph = FALSE) {

  n <- dim(G)[1]
  dag <- matrix(0, n, n)
  dag[ G == 3  &  t(G) == 2 ] <- 1
  isDesc <- transitiveClosure(dag)
  
  if ( is.null(node) )  {
  
    res <- list( isDesc = isDesc )
	
  } else {
    desc <- as.vector( isDesc[, node] )
    desc <- which( desc > 0 )
    Gdesc <- G[c(node, desc), c(node, desc)]
    if ( graph ) {
      if ( length(desc) > 0 )   plotnetwork(Gdesc, titlos = paste("Completed partially directed graph with descendants of ", node ) )
    }
    res <- list(isDesc = isDesc, Gdesc = Gdesc, desc = desc)     
  }
  
  res
}


