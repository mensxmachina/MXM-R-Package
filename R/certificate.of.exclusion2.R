certificate.of.exclusion2 <- function(xIndex, mmpc2object) {
    threshold <- mmpc2object$threshold
    pval.ini <- mmpc2object$univ$pvalue
    pval <- mmpc2object$kapa_pval
   
    info <- list()
    for ( i in 1:length(xIndex) ) {
      if ( sum( mmpc2object$selectedVars == xIndex[i]) ) {
        info[[ i ]] <- paste("Variable", xIndex[i], "has been selected", sep = " ")
        
      } else if ( pval.ini[ xIndex[i] ] > threshold ) {
        info[[ i ]] <- c(0, mmpc2object$univ$pvalue[ xIndex[i] ] )
        names(info[[ i ]])[2] <- "p-value"
      } else {
        j <- 0
        ok <- FALSE
        while ( !ok  & j < length(pval) ) {
          j <- j + 1
          poies <- which( pval[[ j ]][2, ] == xIndex ) 
          if  ( length(poies) > 0 ) {
            a <- which( pval[[ j ]][1, poies] > threshold )
            if ( length(a) > 0 ) {
              info[[ i ]] <- c( pval[[ j ]][-c(1,2), a], pval[[ j ]][1, a] ) 
              names(info[[ i ]])[ length( info[[ i ]] ) ] <- "p-value"
              ok <- TRUE 
            }  ##  end if ( length(a) > 0 ) {
          }  ##  end  if  ( length(poies) > 0 ) {
        }  ##  end  while ( !ok  & j < length(pval) ) {
      }  ##  end  if ( sum( mmpc2object$selectedVars == xIndex[i]) ) {
    }  ## end  for ( i in 1:length(xIndex) ) {
     
  names(info) <- xIndex
  info  
}



