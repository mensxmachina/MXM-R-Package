gSquare = function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels=NULL, hash = FALSE, 
stat_hash=NULL, pvalue_hash=NULL) {
  #Conditional Independence test based on the G test of independence (log likelihood ratio  test)
  csIndex[which(is.na(csIndex))] = 0;
  
  if( hash ) {
    csIndex2 = csIndex[which(csIndex!=0)]
    csindex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if (is.null(stat_hash[[key]]) == FALSE ) {
      stat = stat_hash[[key]];
      pvalue = pvalue_hash[[key]];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  
  xIndex = as.integer(xIndex);
  csIndex = as.integer(csIndex);
  if (length(csIndex) == 1)  {
    if (xIndex == csIndex)  {
      if ( hash ) {    #update hash objects
        stat_hash[[key]] <- 0;     #.set(stat_hash , key , 0)
        pvalue_hash[[key]] <- log(1);     #.set(pvalue_hash , key , 1)
      }
      pvalue = log(1);
      stat = 0;
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
    
  } else {
    csIndex = csIndex[csIndex!=0]
    if (length(csIndex) == 0) {
      csIndex = NULL
    }
  }
  
  if ( is.null(csIndex) ) {
    zz <- cbind(target, dataset[, xIndex] )
    dc <- Rfast::colrange(zz, cont = FALSE)  ##  as.numeric( apply(zz, 2, function(x) { length(unique(x)) } ) )
    mod <- cat.ci(1, 2, 0, zz, type = dc )
    stat <- mod[1]
    pvalue <- mod[2]
  } else {
    zz <- cbind(target, dataset[, c(xIndex, csIndex)] )
    dc <- Rfast::colrange(zz, cont = FALSE)  ##  as.numeric( apply(zz, 2, function(x) { length(unique(x)) } ) )
    #levels for each variable
    mod <- cat.ci(1, 2, c(1:dim(zz)[2])[-c(1:2)], zz, type = dc )
    stat <- mod[1]
    pvalue <- mod[2]
  }
  
  if ( hash ) {
    stat_hash[[key]] <- stat;      #.set(stat_hash , key , stat)
    pvalue_hash[[key]] <- pvalue;     #.set(pvalue_hash , key , pvalue)
  }
  
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
}








