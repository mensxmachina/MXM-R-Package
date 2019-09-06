dag2eg <- function(dag, type = NULL) {

  ## dag is a square matrix corresponding to a DAG
  if ( nrow(dag) != ncol(dag) ) {
    essential <- paste("This matrix is not square")
    
  } else {

    if ( sum( dag >= 2 ) > 0 ) {
      typos = 1  
      dag[ dag == 3 ] <- 0
      dag[ dag == 2 ] <- 1
    } else  typos = 2
    
     ## essential <- ggm::essentialGraph(dag)
      essential <- dag_to_eg(dag)$eg 
      essential[ essential == 2 ] <- 1
 
      if ( typos == 1 ) { 
        eg <- essential + t(essential)
        a <- which(eg == 2) 
        b <- which(eg == 1, arr.ind = TRUE)   
        b <- t( Rfast::colSort( t(b) ) )   ## t( apply(b, 1, sort ) )         
        b <- unique( b )
        if ( nrow(b) > 0 ) {  
          eg[cbind(b[, 2], b[, 1]) ] <- 3
          eg[cbind(b[, 1], b[, 2]) ] <- 2
        }
        eg[ a ] <- 1
        essen <- eg
         
      } else if (typos == 2)    essen <- essential 
  }    
 
  essen
} 
       
     
  
  