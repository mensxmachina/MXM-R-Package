generatefolds <- function(target, nfolds = 10, stratified = TRUE, seed = FALSE) {
  
  a <- paste("Fold", 1:nfolds)
  runs <- sapply(a, function(x) NULL)
  if ( seed )  set.seed(1234)
  
  if ( !stratified ) {
    
    oop <- options(warn = -1) 
    mat <- matrix( sample( length(target) ), ncol = nfolds )
	on.exit( options(oop) )
    for ( i  in 1:c(nfolds - 1) )  runs[[ i ]] <- mat[, i]
    a <- prod( dim(mat) ) - length(target)
    runs[[ nfolds  ]] <- mat[ 1:c(nrow(mat) - a), nfolds ]

  } else {
    
   labs <- unique(target)
   run <- list()
     
   for ( i in 1:length(labs) ) {
     a <- which( target == labs[i] )
     run[[ i ]] <- sample(a) 
   }

   run <- unlist(run) 

    for ( i in 1:length(target) ) {
        k <- i %% nfolds
        if ( k == 0 )  k <- nfolds
        runs[[ k ]] <- c( runs[[ k ]], run[i] )
    }

  }

  runs
}


