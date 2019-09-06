## finds ancestors of some given node(s)
findAncestors <- function(G, node = NULL, graph = FALSE) {

  n <- dim(G)[1]
  dag <- matrix(0, n, n )
  dag[ G == 2  &  t(G) == 3 ] <- 1
  isAnc <- transitiveClosure(dag)
  
  if ( is.null(node) )  {
    res <- list( isAnc = isAnc )
	
  } else {
    anc <- as.vector( isAnc[, node] )
    anc <- which( anc > 0 )
    Ganc <- G[c(node, anc), c(node, anc)]
    
    if ( graph ) {
      if ( length(anc) > 0 )   plotnetwork(Ganc, titlos = paste("Completed partially directed graph with ancestors of ", node ) )
    }

	  res <- list(isAnc = isAnc, Ganc = Ganc, anc = anc)   
	
  } 
  
  res
}


