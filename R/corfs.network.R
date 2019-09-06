corfs.network <- function(x, threshold = 0.05, tolb = 2, tolr = 0.02, stopping = "BIC", symmetry = TRUE, nc = 1) {

  dm <- dim(x)
  n <- dm[1]
  D <- dm[2]
  G <- matrix(0, D, D)
  counter <- 0
  x <- Rfast::standardise(x)
  if ( nc <= 1  || is.null(nc) ) {
    pa <- proc.time()
       for (i in 1:D) {
         id <- 1:D
         id <- id[-i] 
         a <- Rfast::cor.fsreg(x[, i], x[, -i], ystand = FALSE, xstand = FALSE, threshold = threshold, tolb = tolb, tolr = tolr, stopping = stopping)[, 1]
         sel <- id[a]        
         G[i, sel] <- 1 
		     counter <- counter + sum(D - 0:a)
       } 
       runtime <- proc.time() - pa
   
  } else {
    pa <- proc.time() 
    cl <- makePSOCKcluster(nc)
    registerDoParallel(cl)
    sel <- numeric(n)
    mod <- foreach(i = 1:D, .combine = rbind, .export = "cor.fsreg", .packages = "Rfast" ) %dopar% {
	  id <- 1:D
      id <- id[-i] 
      sela <- numeric(D)  
      a <- Rfast::cor.fsreg(x[, i], x[, -i], ystand = FALSE, xstand = FALSE, threshold = threshold, tolb = tolb, tolr = tolr, stopping = stopping)[, 1]
      sel <- id[a]  
      sela[sel] <- 1
      return( sum(D - 0:a), sela)
    }    
    stopCluster(cl)
    G <- as.matrix(mod[, -1])
	  counter <- sum(mod[, 1])
    runtime <- proc.time() - pa
  }	   
  
  diag(G) <- 0	   
  if ( symmetry ) {
    a <- which( G == 1  &  t(G) == 1 ) 
    G[ -a ] <- 0
  } else {
    G <- G + t(G)
    G[ G > 0 ] <- 1
  }

  info <- summary( Rfast::rowsums(G) )
  density <- sum(G) / D / ( D - 1 )  
  if ( is.null( colnames(x) ) ) {
    colnames(G) <- rownames(G) <- paste("X", 1:D, sep = "")
  } else  colnames(G) <- rownames(G) <- colnames(x) 
  
  list(runtime = runtime, density = density, info = info, ntests = counter, G = G)
}