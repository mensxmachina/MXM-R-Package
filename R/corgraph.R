corgraph <- function(dataset, test = "testIndFisher", threshold = 0.01) {
  runtime <- proc.time()
  dm <- dim(dataset)
  n <- dm[1]    ;     p <- dm[2]
  oop <- options(warn = -1) 
  on.exit( options(oop) )
	
  if ( test == "testIndSpearman" ) {
    dataset <- Rfast::colRanks(dataset)
    R <- Rfast::cora(dataset)
    stat <- 0.5 * log( (1 + R)/( (1 - R) ) ) * sqrt(n - 3) / 1.029563
    pvalue <- log(2) + pt( abs(stat), n - 3, lower.tail = FALSE, log.p = TRUE)
    diag(pvalue) <- 0
    R <- NULL
  } else if ( test == "testIndFisher" ) {
    R <- Rfast::cora(dataset)
    stat <- 0.5 * log( (1 + R)/( (1 - R) ) ) * sqrt(n - 3) 
    pvalue <- log(2) + pt( abs(stat), n - 3, lower.tail = FALSE, log.p = TRUE)
    diag(pvalue) <- 0
    R <- NULL
  } else  if ( test == "gSquare" ) {
    dc <- Rfast::colrange(dataset, cont = FALSE)
    stat <- Rfast::g2Test_univariate(dataset, dc)
    pvalue <- pchisq(stat$statistic, stat$df, lower.tail = FALSE, log.p = TRUE)
    pvalue <- Rfast::squareform(pvalue)
  } ## end if ( test == "testIndFisher" ) 
  
  G <- matrix(0, nrow = p, ncol = p)
  G[ pvalue < log(threshold) ] <- 1
  runtime <- proc.time() - runtime
  
  list(runtime = runtime, stat = stat, pvalue = pvalue, G = G)    
}  