mmpc.or <- function(x, max_k = 5, threshold = 0.01, test = "testIndFisher", backward = TRUE, 
                   symmetry = TRUE, ini.pvalue = NULL) {
  dm <- dim(x)
  p <- dm[2]
  n <- dm[1]
  e <- list()
  nam <- colnames(x)
  if ( is.null(nam) )  nam <- paste("Node", 1:p)
  e <- sapply(nam, function(x) NULL)
  G <- matrix(0, p, p)
  ntests <- numeric(p)
  if ( !is.matrix(x) )    x <- as.matrix(x)
  runtime <- proc.time()
  initial.tests <- 0
  ini <- NULL
    
  if ( is.null(ini.pvalue)  &  !is.null(test) ) {
    oop <- options(warn = -1) 
	on.exit( options(oop) )   
    initial.tests <- 0.5 * p * (p - 1)
    if ( test == "testIndSpearman" ) {
      x <- apply(x, 2, rank)
      R <- Rfast::cora(x)
      stat <- 0.5 * log( (1 + R)/( (1 - R) ) ) * sqrt(n - 3) / 1.029563
      ini.pvalue <- log(2) + pt( abs(stat), n - 3, lower.tail = FALSE, log.p = TRUE)
	  diag(ini.pvalue) <- 0
      R <- NULL
      stat <- NULL
    } else if ( test == "testIndFisher" ) {
      R <- Rfast::cora(x)
      stat <- 0.5 * log( (1 + R)/( (1 - R) ) ) * sqrt(n - 3) 
      ini.pvalue <- log(2) + pt( abs(stat), n - 3, lower.tail = FALSE, log.p = TRUE)
	  diag(ini.pvalue) <- 0
      R <- NULL
      stat <- NULL
    } else  if ( test == "gSquare" ) {
      dc <- Rfast::colrange(x, cont = FALSE)
      stat <- Rfast::g2Test_univariate(x, dc)
      ini.pvalue <- pchisq(stat$statistic, stat$df, lower.tail = FALSE, log.p = TRUE)
      ini.pvalue <- Rfast::squareform(ini.pvalue)
    } 
  }
  
  if ( !is.null(ini.pvalue) ) {
    ini <- list()
    ini$stat <- rnorm(p)
  } 
  
  for (j in 1:p) {
    if ( !is.null(ini.pvalue) )   ini$pvalue <- ini.pvalue[j, ]
    a <- MMPC(j, x, max_k = max_k, threshold = threshold, test = test, backward = backward, hash = TRUE, ini = ini)
    if ( !is.null(ini.pvalue) )   ini.pvalue[j, ] <- a@univ$pvalue
    sel <- a@selectedVars
    G[j, sel] <- 1
    ntests[j] <- a@n.tests
    if ( is.null(a@hashObject$pvalue_hash) | length(a@hashObject$pvalue_hash) == 0 ) {
      e[[ j ]] <- NA
    } else  e[[ j ]] <- Rfast::hash2list( as.list.environment( a@hashObject$pvalue_hash ) ) 
  }
  
  dm <- numeric(p)
  for (i in 1:p)   dm[i] = sum( !is.na(e[[ i ]]) )
  poia <- (1:p)[dm>0]
  for (i in poia)   for (j in 1:dm[i])  e[[ i ]][[ j ]] <- c(i, e[[ i ]][[ j ]] )     
  
  a <- unlist(e, recursive = FALSE)
  if ( sum( is.na(a) ) < p ) {
    b <- list()
    lg <- log(threshold)
    kapa <- max_k
    for (i in 1:kapa) b[[ i ]] = list()
    dom <- numeric(kapa)

    for (i in 1:kapa) {
      for ( j in 1:length(a) ) {
         if ( length(a[[ j ]]) - 3 == i ) {
           d <- length(a[[ j ]])
           if ( a[[ j ]][d] > lg ) {
             dom[i] <- dom[i] + 1 
             b[[ i ]][[ dom[i] ]] = c( a[[ j ]][1:(d-1)], 0, a[[ j ]][d])
           }
         }
      }
    }
    
    k <- kapa  
    while ( length( b[[ k ]] ) == 0  & k > 1) {
      b[[ k ]] <- NULL
      k <- k - 1
    }
    
    if (k == 1) {
      if ( length( b[[ k ]] ) == 0 )   b[[ k ]] <- NULL
    }
    
    kapa <- k
    a <- list()
    if ( length(b) > 0 ) { 
      for (i in 1:kapa) {
        d <- length( b[[ i ]] )
        mat <- matrix(0, d, i + 4)
        for (j in 1:d) {
          ela <- as.numeric( b[[ i ]][[ j ]] ) 
          ela[1:2] <- sort(ela[1:2])
          mat[j, ] <- ela 
        }
        a[[ i ]] <- mat  
      }
    }  
  }  else  {
     a <- list()
     kapa <- 0
  }
  
  if ( symmetry ) {
    G <- G + t(G)
    G[ G == 1 ] <- 0
    G[G > 1] <- 1
  } else {
    G <- G + t(G)
    G[ G > 0 ] <- 1
  }
  
  runtime <- proc.time() - runtime
  colnames(G) <- rownames(G) <- nam
  ntests <- c(initial.tests, ntests)
  names(ntests) <- c("univariate tests", nam)
  
  mod <- list()
  mod$ini.pvalue <- ini.pvalue
  mod$kapa <- kapa
  mod$ntests <- ntests
  mod$info <- summary( Rfast::rowsums(G) )
  mod$density <- sum(G) / n / ( n - 1 ) 
  mod$runtime <- runtime
  mod$G <- G
  mod$sepset <- a
  tic <- proc.time()
  mod$G <- pc.or(mod)$G
  mod$runtime.or <- proc.time() - tic
  mod$Gini <- G
  mod
}
 


