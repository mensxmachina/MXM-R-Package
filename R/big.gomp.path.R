big.gomp.path <- function(target = NULL, dataset, tol = c(8, 9, 10), 
                     test = "testIndFisher", method = "ar2") {
  
  tic <- proc.time()
  mod <- MXM::big.gomp(target, dataset, tol = min(tol), test = test, method = method) 
  sel <- mod$res[, 1]
  crit <- mod$res[, 2]
  dm <- dim(mod$res)[1]
  res <- matrix(0, nrow = dm[1], ncol = length(tol) + 1 )
  colnames(res) <- c( paste("tol=", tol, sep = ""), "Deviance")
  res[, 1] <- sel
  res[, dim(res)[2]] <- crit
  
  if ( test == "testIndGamma" | test == "testIndNormLog" | test == "testIndQBinom" | test == "testIndQPois") {
    b <- abs( diff(crit) / mod$phi[-1] )
  } else  b <- abs( diff(crit) )
  
  n <- length(crit)
  for ( i in 2:length(tol) ) {
    ep1 <- which(b < tol[i]) 
    if ( length(ep1) == 0 ) {
      ep <- 1:n
    } else   ep <- 1:min( ep1 )
    res[ep, i] <- sel[ep]     
  }
  runtime <- proc.time() - tic
  list(runtime = runtime, res = res)
}
