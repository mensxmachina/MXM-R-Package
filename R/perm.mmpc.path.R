perm.mmpc.path <- function(target, dataset, wei = NULL, max_ks = NULL, alphas = NULL, test = NULL, user_test = NULL, R = 999, ncores = 1) {

  if( is.null(alphas) )  alphas <- c(0.1, 0.05, 0.01)
  if( is.null(max_ks) )   max_ks <- c(4, 3, 2)  
  
  alphas <- sort(alphas, decreasing = TRUE)
  max_ks <- sort(max_ks, decreasing = TRUE)
  
  nalpha <- length(alphas);
  nmaxk <- length(max_ks);
  
  size <- matrix(0, nalpha, nmaxk)
  bic <- matrix(0, nalpha, nmaxk)
  vars <- list()
  
  iniset <- NULL
  inihash <- NULL
  maxj <- length(max_ks)
  
  tic <- proc.time()
  
  for (i in 1:nalpha) {
    for (j in 1:nmaxk) {

      results <- perm.mmpc(target, dataset, max_k = max_ks[j], threshold = alphas[i], test = test, ini = iniset, wei = wei, 
                           hash = TRUE, hashObject = inihash, R = R, ncores = ncores) 
      iniset <- results@univ
      inihash <- results@hashObject;
      
      a <- mmpc.model(target, dataset, wei = wei, results)$signature 
      
      if ( !is.null(a) ) {
        bic[i, j] <- a[length(a)]    
      } else bic[i, j] <- NULL
      
      size[i, j] <- length( results@selectedVars );    
      k <- (i - 1) * maxj + j
      vars[[ k ]] <-  results@selectedVars
      names(vars)[[ k ]] <- paste("alpha=", alphas[i], " & max_k=", max_ks[j], sep = "")
    }
  }
  
  runtime <- proc.time() - tic
  rownames(bic) <- paste("alpha=", alphas, sep = "")
  colnames(bic) <- paste("max_k=", max_ks, sep = "")
  rownames(size) <- paste("alpha=", alphas, sep = "")
  colnames(size) <- paste("max_k=", max_ks, sep = "")
  
  list(bic = bic, size = size, variables = vars, runtime = runtime)
  
}
