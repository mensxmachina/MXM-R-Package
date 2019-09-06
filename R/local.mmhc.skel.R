local.mmhc.skel = function(dataset, node, max_k = 3, threshold = 0.05, test = "testIndFisher") {
  
  durat <- proc.time()  
  
  d <- dim(dataset)[2]
  G <- matrix(0, d, d)
  if ( is.null( colnames(dataset) ) ) {
    nam  <- colnames(dataset) <- paste("X", 1:d, sep ="")
    rownames(G) <- paste("X", 1:d, sep ="")
    colnames(G) <- paste("X", 1:d, sep ="")
  } else  {
    nam <- colnames(dataset)
    colnames(G) <- rownames(G) <- nam  
  }
  
  a <- MMPC(node, dataset, max_k = max_k, threshold = threshold, test = test, backward = TRUE)
  pct <- a@selectedVars;
  G[node, pct] <- 1
  ntests <- list()
  lista <- list()
  lista[[ node ]] <- pct
  ntests[[ node ]] <- a@n.tests
  
  if ( length(pct) > 0 ) {
    
    for ( i in pct) {
      res <- MMPC(i, dataset, max_k = max_k, threshold = threshold, test = test, backward = TRUE) 
      ntests[[ i ]] <- res@n.tests
      poies <- res@selectedVars
      lista[[ i ]] <- poies
      G[i, poies] <- 1
    }
    s1 <- which( Rfast::colsums(G) > 0 )
    Gloc <- G[s1, s1]
    ntests <- unlist(ntests)
    names(ntests) <- nam[c(node, pct)]
    res <- list()
    b <- c(node, pct)
    k <- length(pct) + 1
    for (i in 1:k )  res[[ i ]] <- lista[[ b[i] ]]
    names(res) <- nam[ c(node, pct)] 
    
  } else {
    Gloc <- matrix(0, nrow = 0, ncol = 0)
    ntests <- ntests
    names(ntests) <- nam[node]
    res <- list()
    res[[ 1 ]] <- numeric(0)
    names(res) <- nam[node]
  }
  
  runtime <- proc.time() - durat   
  list(runtime = runtime, ntests = ntests, res = res, Gloc = Gloc )
}