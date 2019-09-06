equivdags <- function(g1, g2) {

  if (  sum( dim(g1) - dim(g2) ) != 0 ) {
    mes <- paste("The dimensions do not match")
    apo <- FALSE

  } else {
    
    a1 <- g1     ;       a2 <- g2
    a1[ a1 > 0 ] <- 1
    a2[ a2 > 0 ] <- 1    
    a <- sum(a1 - a2)
    
    if (a != 0 ) {
      mes <- paste("The two DAgs do not share the same adjacencies")
      apo <- FALSE

    } else {
    
      b1 <- which( colSums(g1 == 2) > 1 )
      b2 <- which( colSums(g2 == 2) > 1 )
      
      if ( length(b1) != length(b2) ) {
        mes <- paste("The two DAgs do not share the same number of unshilded colliders")
        apo <- FALSE

      } else if ( sum(b1 - b2) != 0 ) {
        mes <- paste("The two DAgs do not share the same number of unshilded colliders")
        apo <- FALSE
       
      } else {
        
        d1 <- as.matrix( g1[, b1] )   ;    d2 <- as.matrix( g2[, b2] )
        d1[ d1 != 2 ] <- 0     ;    d2[ d2 != 2 ] <- 0 
        
        if ( sum(d1 - d2) != 0 ) {
          mes <- paste("The two DAgs do not share the same number of unshilded colliders")
          apo <- FALSE 

        } else {
          mes <- paste("The two DAgs are Markov equivalent")
          apo <- TRUE
        }

      }

     }

   } 
   
   list(apofasi = apo, mes = mes)

}
