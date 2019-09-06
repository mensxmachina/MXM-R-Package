mmpc2 <- function(target, dataset, max_k = 3, threshold = 0.05, test = "testIndLogistic", ini = NULL, wei = NULL, ncores = 1, backward = FALSE) {
  
  runtime <- proc.time()

 #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any(is.na(dataset)) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
      for ( i in poia )  {
        xi <- dataset[, i]
        if ( is.numeric(xi) ) {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) ) {
          xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        }
        dataset[, i] <- xi
      }
    }
  }

 av_tests = c("testIndReg", "testIndBeta", "censIndCR", "censIndWR", "testIndClogit", "testIndOrdinal",
              "testIndLogistic", "testIndPois", "testIndNB", "testIndBinom", "auto", "testIndZIP", 
              "testIndRQ", "testIndGamma", "testIndNormLog", "testIndTobit", "testIndQPois", "censIndCR", 
              "censIndWR", "censIndER", "censIndLLR", "testdIndQBinom", "testIndMMReg", "testIndMultinom", 
              "testIndIGreg", "testIndSPML", NULL);

 ci_test <- test
 if (length(test) == 1) {      #avoid vectors, matrices etc
   test <- match.arg(test, av_tests, TRUE);
   if (test == "testIndFisher") {
     test <- testIndFisher;
     
   } else if (test == "testIndMMFisher") {
     test <- testIndMMFisher;
     
   } else if (test == "testIndMMReg") {
     test <- testIndMMReg;
     
   } else if (test == "testIndSpearman")  {
     target <- rank(target)
     dataset <- Rfast::colRanks(dataset)  
     test <- testIndSpearman;  ## Spearman is Pearson on the ranks of the data
     
   } else if (test == "testIndReg")  {  ## It uMMPC the F test
     test <- testIndReg
     
   }  else if (test == "testIndMVreg") {
     if ( min(target) > 0 & sd( Rfast::rowsums(target) ) == 0 )  target = log( target[, -1]/target[, 1] ) 
     test <- testIndMVreg;
     
   } else if (test == "testIndBeta") {
     test <- testIndBeta;
     
   } else if (test == "testIndRQ") {  ## quantile regression
     #an einai posostiaio target
     test <- testIndRQ;
     
   } else if (test == "testIndIGreg") { ## Inverse Gaussian regression
     test <- testIndIGreg;
     
   } else if (test == "testIndPois") { ## Poisson regression
     test <- testIndPois;
     
   } else if (test == "testIndNB") { ## Negative binomial regression
     test <- testIndNB;
     
   } else if (test == "testIndGamma") {  ## Gamma regression
     test <- testIndGamma;
     
   } else if (test == "testIndNormLog") { ## Normal regression with a log link
     test <- testIndNormLog;
     
   } else if (test == "testIndZIP") { ## Zero inflated Poisson regression
     test <- testIndZIP;
     
   } else if (test == "testIndTobit") { ## Tobit regression
     test <- testIndTobit;
     
   } else if (test == "censIndCR") {
     test <- censIndCR;
     
   } else if (test == "censIndWR") {
     test <- censIndWR;
     
   } else if (test == "censIndER") {
     test <- censIndER;
     
   } else if (test == "censIndLLR") {
     test <- censIndLLR;
     
   } else if (test == "testIndClogit") {
     test <- testIndClogit;
     
   } else if (test == "testIndBinom") {
     test <- testIndBinom;
     
   } else if (test == "testIndLogistic") {
     test <- testIndLogistic;
     
   } else if (test == "testIndMultinom") {
     test <- testIndMultinom;
     
   } else if (test == "testIndOrdinal") {
     test <- testIndOrdinal;
     
   } else if (test == "testIndQBinom") {
     test <- testIndQBinom;
     
   } else if (test == "testIndQPois") {
     test <- testIndQPois;
     
   } else if (test == "gSquare") {
     test <- gSquare;
     
   } else if (test == "testIndSPML") {
     test <- testIndSPML
     if ( !is.matrix(target) )   target <- cbind( cos(target), sin(target) )
   }
   #more tests here
 } else {
   stop('invalid test option');
 }
 
  dataset <- as.data.frame(dataset)  
  n.tests <- 0
  alpha <- log(threshold)
  kapa_pval <- NULL
  
  varsize <- dim(dataset)[2]
  if ( ( typeof(max_k) != "double" ) | max_k < 1 )   stop('invalid max_k option');
  if ( max_k > varsize )   max_k <- varsize;
  if ( (typeof(threshold) != "double") | threshold <= 0 | threshold > 1 )   stop('invalid threshold option');
  
  if ( is.null(ini) ) {
    mod <- univregs(target, dataset, test = test, wei = wei, ncores = ncores) 
    pval <- mod$pvalue
    n.tests <- dim(dataset)[2]
  } else pval <- ini$pvalue
  vars <- which(pval < alpha) 
  if ( length(vars) > 0 ) {
    sela <- which.min(pval) 
  } else  sela <- vars
  vars <- setdiff(vars, sela)
  
  ## 1 selected
  if ( length(vars) > 0  &  max_k >= 1 ) {
    a <- paste("kappa=", 1:max_k, sep = "")
    kapa_pval <- sapply(a, function(x) NULL)
    pval2 <- cond.regs(target, dataset, xIndex = vars, csIndex = sela, test = test, wei = wei, ncores = 1)$pvalue
    kapa_pval[[ 1 ]] <- rbind(pval2[vars], vars, sela)
    n.tests <- n.tests + length(vars)
    pval[vars] <- pmax(pval[vars], pval2[vars]) 
    ide <- which(pval[vars] < alpha)
    vars <- vars[ide]
    sel <- which.min(pval[vars])
    sela <- c(sela, vars[sel] )
    vars <- setdiff(vars, vars[sel])  
  }  ## end  if ( length(vars) > 0  &  max_k >= 1 ) {

  ## 2 selected
  if ( length(vars) > 0  &  max_k >= 1 ) {
    pval2 <- cond.regs(target, dataset, xIndex = vars, csIndex = sela[2], test = test, wei = wei, ncores = 1)$pvalue
    kapa_pval[[ 1 ]] <- cbind( kapa_pval[[ 1 ]], rbind(pval2[vars], vars, sela[2]) )
    n.tests <- n.tests + length(vars)
    pval[vars] <- pmax(pval[vars], pval2[vars]) 
    ide <- which(pval[vars] < alpha)
    vars <- vars[ide]
    if ( length(vars) > 0  &  max_k >= 2 ) {
      pval2 <- cond.regs(target, dataset, xIndex = vars, csIndex = sela, test = test, wei = wei, ncores = 1)$pvalue
      kapa_pval[[ 2 ]] <- rbind( pval2[vars], vars, matrix( rep(sela, length(vars)), ncol = length(vars) ) )
      pval[vars] <- pmax(pval[vars], pval2[vars]) 
      ide <- which(pval[vars] < alpha)
      vars <- vars[ide]
    }  ## end  if ( max_k >= 2 ) {
    sel <- which.min(pval[vars])
    sela <- c(sela, vars[sel] )
    vars <- setdiff(vars, vars[sel])    
    dm <- length(sela)
  }  ## end  if ( length(vars) > 0  &  max_k >= 1 ) {

  ## 3 selected  
  if ( length(vars) > 0  &  max_k >= 1 ) {
    if ( max_k >= 1 ) {
      pval2 <- cond.regs(target, dataset, xIndex = vars, csIndex = sela[3], test = test, wei = wei, ncores = 1)$pvalue
      kapa_pval[[ 1 ]] <- cbind( kapa_pval[[ 1 ]], rbind(pval2[vars], vars, sela[3]) )
      n.tests <- n.tests + length(vars)
      pval[vars] <- pmax(pval[vars], pval2[vars]) 
      ide <- which(pval[vars] < alpha)
      vars <- vars[ide]
    }  ## end  if ( max_k >= 1 ) {
    if ( length(vars) > 0  &  max_k >= 2 ) {
      cand <- Rfast::comb_n(sela, 2)[, - 1]
      pval2 <- cond.regs(target, dataset, xIndex = vars, csIndex = cand[, 1], test = test, wei = wei, ncores = 1)$pvalue
      kapa_pval[[ 2 ]] <- cbind( kapa_pval[[ 2 ]], rbind( pval2[vars], vars, matrix( rep(cand[, 1], length(vars)), ncol = length(vars) ) ) )
      n.tests <- n.tests + length(vars)
      pval[vars] <- pmax(pval[vars], pval2[vars]) 
      ide <- which(pval[vars] < alpha)
      vars <- vars[ide]
      pval2 <- cond.regs(target, dataset, xIndex = vars, csIndex = cand[, 2], test = test, wei = wei, ncores = 1)$pvalue
      kapa_pval[[ 2 ]] <- cbind( kapa_pval[[ 2 ]], rbind( pval2[vars], vars, matrix( rep(cand[, 2], length(vars)), ncol = length(vars) ) ) )
      n.tests <- n.tests + length(vars)
      pval[vars] <- pmax(pval[vars], pval2[vars]) 
      ide <- which(pval[vars] < alpha)
      vars <- vars[ide]
    }  ## end  if ( max_k >= 2 ) {
    if ( length(vars) > 0  &  max_k >= 3 ) {
      pval2 <- cond.regs(target, dataset, xIndex = vars, csIndex = sela, test = test, wei = wei, ncores = 1)$pvalue
      kapa_pval[[ 3 ]] <- rbind( pval2[vars], vars, matrix( rep(sela, length(vars)), ncol = length(vars) ) )
      n.tests <- n.tests + length(vars)
      pval[vars] <- pmax(pval[vars], pval2[vars]) 
      ide <- which(pval[vars] < alpha)
      vars <- vars[ide]
    }  ## end  if ( max_k >= 3 ) {
    sel <- which.min(pval[vars])
    sela <- c(sela, vars[sel] )
    vars <- setdiff(vars, vars[sel])  
  }  ## end  if ( length(vars) > 0  &  max_k >= 2 ) {  
  
  ## 4 selected
  while ( length(vars) > 0  &  max_k >= 1 ) {
    dm <- length(sela)
    pval2 <- cond.regs(target, dataset, xIndex = vars, csIndex = sela[dm], test = test, wei = wei, ncores = 1)$pvalue
    kapa_pval[[ 1 ]] <- cbind( kapa_pval[[ 1 ]], rbind(pval2[vars], vars, sela[dm]) )
    n.tests <- n.tests + length(vars)
    pval[vars] <- pmax(pval[vars], pval2[vars]) 
    ide <- which(pval[vars] < alpha)
    vars <- vars[ide]
    if ( length(vars) > 0  &  max_k >= 2 ) {
      for ( i in 2:max_k ) {  
        cand <- Rfast::comb_n(sort(sela), i)
        cand <- cand[, which(cand == sela[dm], arr.ind = TRUE)[, 2], drop = FALSE ]
        j <- 1
        #for ( j in 1:dim(cand)[2] ) {
        while ( length(vars) > 0  &  j <= dim(cand)[2] ) {
          pval2 <- cond.regs(target, dataset, xIndex = vars, csIndex = cand[, j], test = test, wei = wei, ncores = 1)$pvalue
          kapa_pval[[ i ]] <- cbind( kapa_pval[[ i ]], rbind( pval2[vars], vars, matrix( rep(cand[, j], length(vars)), ncol = length(vars) ) ) )
          n.tests <- n.tests + length(vars)
          pval[vars] <- pmax(pval[vars], pval2[vars]) 
          ide <- which(pval[vars] < alpha)
          vars <- vars[ide]
          j <- j + 1
        }  ## end  while ( length(vars) > 0  &  j <= dim(cand)[2] ) {
      }  ## end  for ( i in 2:max_k ) {  
    }  ## end  if ( max_k >= 2 ) {
    sel <- which.min(pval[vars])
    sela <- c(sela, vars[sel] )
    vars <- setdiff(vars, vars[sel])  
  } ## end  while ( length(vars) > 0  &  max_k >= 1 ) {
    
  runtime <- proc.time() - runtime
    
  if ( backward  & length( sela ) > 0  ) {
    tic <- proc.time()
    pv <- pval[sela]
    sela <- sela[ order(pv) ]
    bc <- mmpcbackphase(target, dataset[, sela, drop = FALSE], test = test, wei = wei, max_k = max_k, threshold = threshold)
    met <- bc$met
    sela <- sela[met]
    pval[sela] <- bc$pvalues
    n.tests <- n.tests + bc$counter
    runtime <- runtime + proc.time() - tic
  }  
    
 list(selectedVars = sela, pvalues = pval, univ = mod$pvalue, kapa_pval = kapa_pval, max_k = max_k, threshold = alpha, 
      n.tests = n.tests, runtime = runtime, test = ci_test)
}





