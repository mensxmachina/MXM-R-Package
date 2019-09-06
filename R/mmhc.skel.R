mmhc.skel <- function(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher", type = "MMPC", hash = FALSE, 
                      backward = TRUE, symmetry = TRUE, nc = 1, ini.pvalue = NULL) {
  ## dataset is either conitnuous or categorical data  
  ## max_k is the maximum number of variables upon which to condition
  ## threshold is the level of significance to reject the independence
  ## test can be either testIndFisher (default) or testIndSpearman for continuous data
  ## OR gSquare (default) for categorical data
  ## nc is the number of cores to use, set to 1 by default
  dm <- dim(dataset)
  n <- dm[2]
  m <- dm[1]
  G <- matrix(0, n, n)
  pvalue <- matrix(1, n, n)
  ntests <- numeric(n)
  initial.tests <- 0
  ini <- NULL
  nam <- colnames(dataset)
  initial.tests <- 0
  ini.pval <- matrix(0, n, n)
  
  if ( !is.null(test) ) {
    if ( is.null(ini.pvalue)  & ( test == "testIndSpearman" | test == "testIndFisher" | test == "gSquare") ) {
	  oop <- options(warn = -1) 
	  on.exit( options(oop) ) 
      if ( !is.matrix(dataset) )   dataset <- as.matrix(dataset)
      initial.tests <- 0.5 * n * (n - 1)
      if ( test == "testIndSpearman" ) {
        dataset <- Rfast::colRanks(dataset)
        R <- Rfast::cora(dataset)
        stat <- 0.5 * log( (1 + R)/( (1 - R) ) ) * sqrt(m - 3) / 1.029563
        ini.pvalue <- log(2) + pt( abs(stat), m - 3, lower.tail = FALSE, log.p = TRUE)
        diag(ini.pvalue) <- 0
        R <- NULL
        stat <- NULL
      } else if ( test == "testIndFisher" ) {
        R <- Rfast::cora(dataset)
        stat <- 0.5 * log( (1 + R)/( (1 - R) ) ) * sqrt(m - 3) 
        ini.pvalue <- log(2) + pt( abs(stat), m - 3, lower.tail = FALSE, log.p = TRUE)
        diag(ini.pvalue) <- 0
        R <- NULL
        stat <- NULL
      } else  if ( test == "gSquare" ) {
        dc <- Rfast::colrange(dataset, cont = FALSE)
        stat <- Rfast::g2Test_univariate(dataset, dc)
        ini.pvalue <- pchisq(stat$statistic, stat$df, lower.tail = FALSE, log.p = TRUE)
        ini.pvalue <- Rfast::squareform(ini.pvalue)
      }  
    }  ## end if ( is.null(ini.pvalue)  & ( test == "testIndSpearman" | test == "testIndFisher" | test == "gSquare") )
  }  ## end if ( !is.null(test) )
  ############
  #### MMPC
  ############
  if ( type == "MMPC" ) {
    ms <- NULL
    if ( !is.null(ini.pvalue) ) {
      ini <- list()
      ini$stat <- rnorm(n)
    }  
    
    if (nc <= 1  ||  is.null(nc) ) {
      pa <- proc.time()
      for (i in 1:n) {
        if ( !is.null(ini.pvalue) )   ini$pvalue <- ini.pvalue[i, ]
        a <- MMPC(i, dataset, max_k = max_k, test = test, threshold = threshold, hash = hash, backward = backward, ini = ini )
        ini.pval[i, ] <- a@univ$pvalue
        pvalue[i, ] <- a@pvalues
        sel <- a@selectedVars
        G[i, sel] <- 1
        ntests[i] <- a@n.tests
      } 
      runtime <- proc.time() - pa
    
    } else {
    
      pa <- proc.time() 
      cl <- makePSOCKcluster(nc)
      registerDoParallel(cl)
      mod <- foreach(i = 1:n, .combine = rbind, .export = c("MMPC") ) %dopar% {
        sel <- numeric(n)
        if ( !is.null(ini.pvalue) )   ini$pvalue <- ini.pvalue[i, ]
        a <- MMPC(i, dataset, max_k = max_k, test = test, threshold = threshold, hash = hash, backward = backward, ini = ini )
        sel[a@selectedVars] <- 1
        return( c(a@n.tests, sel, a@pvalues, a@univ$pvalue) )
      }
      
      stopCluster(cl)
      G <- as.matrix(mod[, 2:(n + 1) ])
      pvalue <- as.matrix( mod[ , (n + 2):(2 * n + 1) ] )
      ini.pval <- as.matrix( mod[, (2 * n + 2):(3 * n + 1) ] )
      ntests <- as.vector(mod[, 1])
      runtime <- proc.time() - pa
    }
    
  } else if (type == "SES") {
    ms <- numeric(n)
    if (nc <= 1  ||  is.null(nc) ) {
      if ( !is.null(ini.pvalue) )   ini$stat <- rnorm(n)
      
      pa <- proc.time()
      for (i in 1:n) {
        if ( !is.null(ini.pvalue) )   ini$pvalue <- ini.pvalue[i, ]
        a <- SES(i, dataset, max_k = max_k, test = test, threshold = threshold, hash = hash, ini = ini)
        ini.pval[i, ] <- a@univ$pvalue
        poies <- a@signatures
        ms[i] <- dim(poies)[1]
        sel <- unique(poies)
        G[i, sel] <- 1 
        ntests[i] <- a@n.tests
        pvalue[i, ] <- a@pvalues
      } 
      runtime <- proc.time() - pa
        
    } else {
        
      pa <- proc.time() 
      cl <- makePSOCKcluster(nc)
      registerDoParallel(cl)
      if ( !is.null(ini.pvalue) )   ini$stat <- rnorm(n)
      mod <- foreach(i = 1:n, .combine = rbind, .export = c("SES") ) %dopar% {
        ## arguments order for any CI test are fixed
        sel <- numeric(n)
        if ( !is.null(ini.pvalue) )   ini$pvalue <- ini.pvalue[i, ]
        a <- SES(i, dataset, max_k = max_k, test = test, threshold = threshold, hash = hash, ini = ini)
        poies <- a@signatures
        sel[ unique(poies) ] <- 1
        return( c(a@n.tests, dim(poies)[1], sel, a@pvalues, a@univ$pvalue) ) 
      }
        
      stopCluster(cl)
      G <- as.matrix(mod[, 3:(n + 2) ])
      pvalue <- as.matrix( mod[ , (n + 3):(2 * n + 2) ] )
      ini.pval <- as.matrix( mod[, (2 * n + 3):(3 * n + 2) ] )
      ntests <- as.vector(mod[, 1])
      ms <- as.vector(mod[, 2])
      runtime <- proc.time() - pa
    }
    ms <- rbind( which(ms>1), ms[ms>1] )
    rownames(ms) <- c("Nodes", "No of equivalent signatures")  
  }  ## end if (type == "MMPC")    
  
  diag(G) <- 0
  
  if ( symmetry ) {
    a <- which( G == 1  &  t(G) == 1 ) 
    G[ -a ] <- 0
    pvalue <- pmax( pvalue, t(pvalue) )
  } else {
    G <- G + t(G)
    G[ G > 0 ] <- 1
    pvalue <- pmin( pvalue, t(pvalue) )
  }
  
  info <- summary( Rfast::rowsums(G) )
  density <- sum(G) / n / ( n - 1 ) 
  
  if ( is.null( nam) ) {
    colnames(G) <- rownames(G) <- paste("X", 1:n, sep = "")
    colnames(pvalue) <- rownames(pvalue) <- paste("X", 1:n, sep = "")
    names(ntests) <- paste("X", 1:n, sep = "")
  } else  {
    colnames(G) <- rownames(G) <- nam
    colnames(pvalue) <- rownames(pvalue) <- nam
    names(ntests) <- nam
  } 
  
  ntests <- c(initial.tests, ntests)
  names(ntests) <- c("univariate tests", nam)
  if ( is.null(ini.pvalue) )  ini.pvalue <- ini.pval
  ini.pval <- NULL
  list(runtime = runtime, density = density, info = info, ms = ms, ntests = ntests, ini.pvalue = ini.pvalue, pvalue = pvalue, G = G)
}