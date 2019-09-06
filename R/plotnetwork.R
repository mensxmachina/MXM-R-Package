#####################
##### Plot a graph using the adjacency matrix
#####
####################
plotnetwork <- function(G, titlos = NULL) {

  if ( sum(G) > 0 ) {
  
    if ( is.null( colnames(G) ) ) {
      nam <- colnames(G) <- rownames(G) <- paste("X", 1:dim(G)[2])
    }  else  nam <- colnames(G)
  
    nodes <- data.frame(id = nam)
  
    if ( sum( abs( G - t(G) ) ) == 0 ) { ## no orientations
      ed <- which(G > 0, arr.ind = TRUE) 
      ed <- unique( t( Rfast::colSort( t(ed) ) ) )
      edges <- data.frame( from = nam[ed[, 1]], to = nam[ed[, 2]], type = numeric( dim(ed)[1] ) )
    
    } else {
      ed1 <- which(G == 2, arr.ind = TRUE)
      ed0 <- which(G == 1, arr.ind = TRUE)
      ed0 <- unique( t( Rfast::colSort( t(ed0) ) ) )
      ed <- rbind(ed1, ed0)
      edges <- data.frame( from = nam[ed[, 1]], to = nam[ed[, 2]], type = c(numeric( dim(ed1)[1] ) + 1, numeric( dim(ed0)[1] )) )
    }

    id = paste(nodes[1:dim(nodes)[1], 1])

    nodesprint <- data.frame( id = id,
                           title = id, #on hover
                           label = id, #current label
                           shape = c("dot"),
                           mass = 1, 
                           size = 40,
                           shadow = TRUE,
                           #color = c("#1abc9c"),
                           physics = TRUE, font = "25px" )

    from <- edges$from
    to <- edges$to
    type <- edges$type
    a <- type
    a[type == 1] = "to"
    a[type == 0] = ""

    edgesprint <- data.frame(from = from, 
                         to = to,
                         #title = weight, #on hover
                         smooth = FALSE,
                         width = 3,
                         shadow = TRUE,
                         arrows = a, #determine the type of the edges here with a vector of length as the number of edges
                         color = list(highlight = "red", hover = "green"),
                         physics = FALSE
    )

  } else  edges <- NULL      

if(is.null(edges)){
  return();
}else{
  
    if ( is.null(titlos) )   titlos <- "Network" 
    visNetwork::visNetwork(nodesprint, edgesprint, width = "100%", height = "600px", main = titlos) %>% 
    visNetwork::visPhysics(minVelocity = 1, maxVelocity = 50,repulsion = list(nodeDistance = 200, centralGravity = 10000), barnesHut = list(gravitationalConstant = -16000)) %>%
    visNetwork::visInteraction(navigationButtons = TRUE,hover = TRUE, hoverConnectedEdges = TRUE) %>%
    visNetwork::visOptions(nodesIdSelection = list(enabled = TRUE))
}

}
