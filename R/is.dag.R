is.dag <- function(dag) {
  if ( any( dag == 2 ) ) {
    dag[ dag != 2 ] <- 0
    dag[ dag == 2 ] <- 1
  } 
  a <- Rfast::topological_sort(dag)
  any( !is.na(a) )
} 
       
     
  
  