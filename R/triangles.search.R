triangles.search <- function(G) {
  n <- dim(G)[1]
  trigwna <- matrix(0, 0, 3)
  a <- vector("list", n)
  for (i in 1:n) {
    
    a[[ i ]] <- which( G[i, ] > 0 )
    if ( length( a[[ i ]] > 0 ) ) {
       for (j in 1:length( a[[ i ]]) ) {
         b1 <- intersect(which( G[a[[ i ]][j],  ] > 0 ), a[[ i ]]) 
         if ( length(b1) > 0 )  {
           for ( k in 1:length(b1) )   trigwna <- rbind(trigwna, c(i, a[[ i ]][j], b1[k]) )
         }  ## end if ( length(b1) > 0 )   
       }   ## end for (j in 1:length( a[[ i ]]) )  
    }   ## end if ( length( a[[ i ]] > 0 ) )
  }  ## end for (i in 1:n) 
  trigwna <- Rfast::rowSort(trigwna)
  unique(trigwna)
}