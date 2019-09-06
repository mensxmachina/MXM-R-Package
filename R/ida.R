ida <- function(x, y, G, dataset) {
  parx <- which(G[, x] == 2)
  poss.parx <- which(G[, x] == 1)
  t1 <- length(parx)
  t2 <- length(poss.parx)
  mess<- NULL
  
  if (G[x, y] == 2  &  t1 > 0 & t2 == 0) {  ## clear situation
    tc <- summary( lm(dataset[, y] ~ dataset[, x] + dataset[, parx]) )[[ 4 ]][2, ] 
  
  } else if ( G[x, y] == 3 ) {  ## clear situation, take cases
    tc <- matrix(0, 1, 4)
    
  } else if ( sum(findDescendants(G, x)$desc == y) == 1 ) {  ## clear situation, take cases
    tc <- matrix(0, 1, 4)
  
  } else if ( G[x, y] == 1 ) {  ## not clear situation, take cases
    if ( length( G[, y] == 2 ) > 0 ) {
      tc <- matrix(0, 1, 4)

    } else if ( t1 == 0  &  t2 == 0) {
      tc <- summary( lm(dataset[, y] ~ dataset[, x]) )[[ 4 ]][2, ] 
      mess <- "No parents of x"
      
    } else if ( t1 > 0 ) {
      tc <- summary( lm(dataset[, y] ~ dataset[, x] + dataset[, parx]) )[[ 4 ]][2, ] 
      
    } else if ( t1 == 0  &  t2 > 0) {
      tc <- matrix(0, 0, 4)
      for (i in 1:t2) {
        tc <- rbind(tc, summary( lm(dataset[, y] ~ dataset[, x] + dataset[, poss.parx(i)]) )[[ 4 ]][2, ] )
      }
    }
    
  } else if ( G[x, y] == 2  |  G[x, y] == 0 ) {  ## not clear situation, take cases

    if ( t1 == 0  &  t2 == 0) { 
      tc <- summary( lm(dataset[, y] ~ dataset[, x]) )[[ 4 ]][2, ] 
      mess <- "No parents of x"
      
    } else if ( t1 > 0 ) {  
      tc <- summary( lm(dataset[, y] ~ dataset[, x] + dataset[, parx]) )[[ 4 ]][2, ] 
      
    } else if ( t1 == 0  &  t2 > 0) {  
      tc <- matrix(0, 0, 4)
      for (i in 1:t2) {
        tc <- rbind(tc, summary( lm(dataset[, y] ~ dataset[, x] + dataset[, poss.parx(i)]) )[[ 4 ]][2, ] )
      }
    }
  }  
  list(tc = tc, mess = mess) 
}