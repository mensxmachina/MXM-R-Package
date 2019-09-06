certificate.of.exclusion <- function(xIndex, sesObject = NULL, mmpcObject = NULL) {
   
  if ( !is.null(sesObject) ) {  ## SES
    threshold <- sesObject@threshold
    pval <- Rfast::hash2list( as.list.environment( sesObject@hashObject$pvalue_hash ) )
    stat <- Rfast::hash2list( as.list.environment( sesObject@hashObject$stat_hash ) )
    info <- list()
    for ( i in 1:length(xIndex) ) {
      if ( sum( sesObject@signatures == xIndex[i]) ) {
         info[[ i ]] <- paste("Variable", xIndex[i], "has been selected", sep = " ")
         
      } else if ( sesObject@univ$pvalue[ xIndex[i] ] > threshold ) {
          info[[ i ]] <- c(0, sesObject@univ$stat[ xIndex[i] ], sesObject@univ$pvalue[ xIndex[i] ] )
          names(info[[ i ]])[2:3] <- c("statistic", "p-value")
      } else {
        for ( j in 1:length(pval) ) {
          if ( pval[[ j ]][1] == xIndex[i] ) {
            d <- length(pval[[ j ]])
            if ( pval[[ j ]][d] > threshold ) {             
              info[[ i ]] <- c( stat[[ j ]][-1], pval[[ j ]][ d ] )
              names(info[[ i ]])[c(d - 1, d)] <- c("statistic", "p-value")
            }
          }  ## end if ( pval[[ j ]][1] == xIndex[i] )
        }  ## end for ( j in 1:length(pval) )    
      }  ## end else
    }  ## end for ( i in 1:length(xIndex) )
    
  } else {  ## MMPC
    threshold <- mmpcObject@threshold
    pval <- Rfast::hash2list( as.list.environment( mmpcObject@hashObject$pvalue_hash ) )
    stat <- Rfast::hash2list( as.list.environment( mmpcObject@hashObject$stat_hash ) )
    info <- list()
    for ( i in 1:length(xIndex) ) {
      if ( sum( mmpcObject@selectedVars == xIndex[i]) ) {
        info[[ i ]] <- paste("Variable", xIndex[i], "has been selected", sep = " ")
        
      } else if ( mmpcObject@univ$pvalue[ xIndex[i] ] > threshold ) {
        info[[ i ]] <- c(0, mmpcObject@univ$stat[ xIndex[i] ], mmpcObject@univ$pvalue[ xIndex[i] ] )
        names(info[[ i ]])[2:3] <- c("statistic", "p-value")
      } else {
        for ( j in 1:length(pval) ) {
		      if ( pval[[ j ]][1] == xIndex[i] ) {
            d <- length(pval[[ j ]])
            if ( pval[[ j ]][d] > threshold) {
              info[[ i ]] <- c( stat[[ j ]][-1], pval[[ j ]][ d ] )
              names(info[[ i ]])[c(d - 1, d)] <- c("statistic", "p-value")
            }
          }  ## end if ( pval[[ j ]][1] == xIndex[i] )
        }  ## end for ( j in 1:length(pval) )    
      }  ## end else
    }  ## end for ( i in 1:length(xIndex) )
    
  }  ## end if ( !is.null(sesObject) )
  
  names(info) <- xIndex
  info  
}