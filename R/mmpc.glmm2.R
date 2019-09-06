mmpc.glmm2 <- function(target, reps = NULL, group, dataset, max_k = 3, threshold = 0.05, test = NULL, ini = NULL, wei = NULL, slopes = FALSE, ncores = 1) {
  
  runtime <- proc.time()
  
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any(is.na(dataset)) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
  }
 
  la <- length( unique( as.numeric(target) ) )
  
    #auto detect independence test in case of not defined by the user and the test is null or auto
    if ( is.null(test) )  {
      if ( survival::is.Surv(target) ) {
        test <- "testIndGLMMCR"
      } else if ( la > 2 & sum(target - round(target) ) != 0  &  is.null(wei)  &  !slopes  &  is.null(reps) ) {
        test <- "testIndLMM"
      } else if (la == 2) {
        test <- "testIndGLMMLogistic"
      } else if ( la > 2  &  sum( target - round(target) ) == 0 ) {
        test <- "testIndGLMMPois"
      } else if ( la > 2  &  is.ordered(target) ) {
        test <- "testIndGLMMOrdinal"
      } else test <- "testIndGLMMReg"
    }  
    #available conditional independence tests
    av_tests = c("testIndGLMMReg", "testIndGLMMLogistic", "testIndGLMMPois", "testIndGLMMGamma", 
                 "testIndGLMMOrdinal", "testIndGLMMNormLog", "testIndLMM", "testIndGLMMCR", "auto",  NULL);
    
    ci_test <- test
    if ( length(test) == 1 ) {   #avoid vectors, matrices etc
      test <- match.arg(test, av_tests, TRUE);
      #convert to closure type
      if ( test == "testIndGLMMReg" ) {
        test <- testIndGLMMReg;
      } else if ( test == "testIndGLMMReg"  &  is.null(wei)  &  (!slopes) ) {
        test <- testIndLMM;
      } else if ( test == "testIndGLMMLogistic" ) {
        test <- testIndGLMMLogistic;
      } else if ( test == "testIndGLMMPois" ) {
        test <- testIndGLMMPois;
      } else if ( test == "testIndGLMMGamma" ) {
        test <- testIndGLMMGamma;
      } else if ( test == "testIndGLMMNormLog" ) {
        test <- testIndGLMMNormLog;
      } else if ( test == "testIndGLMMOrdinal" ) {
        test <- testIndGLMMOrdinal;
      } else if ( test == "testIndGLMMCR" ) {
        test <- testIndGLMMCR;
      } else if ( test == "testIndLMM" ) {
        test <- testIndLMM
      }
      
    } else   stop('invalid test option');

  n.tests <- 0
  alpha <- log(threshold)
  kapa_pval <- NULL
  
  varsize <- dim(dataset)[2]
  if ( ( typeof(max_k) != "double" ) | max_k < 1 )   stop('invalid max_k option');
  if ( max_k > varsize )   max_k <- varsize;
  if ( (typeof(threshold) != "double") | threshold <= 0 | threshold > 1 )   stop('invalid threshold option');
  
  if ( is.null(ini) ) {
    mod <- glmm.univregs(target = target, reps = reps, id = group, dataset = dataset, targetID = -1, test = test, wei = wei, 
                              slopes = slopes, ncores = ncores)
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
    pval2 <- glmm.condregs(target = target, reps = reps, id = group, dataset = dataset, xIndex = vars, csIndex = sela, test = test, 
                           wei = wei, slopes = slopes, ncores = 1)$pvalue 
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
    pval2 <- glmm.condregs(target = target, reps = reps, id = group, dataset = dataset, xIndex = vars, csIndex = sela[2], test = test, 
                           wei = wei, slopes = slopes, ncores = 1)$pvalue 
    kapa_pval[[ 1 ]] <- cbind( kapa_pval[[ 1 ]], rbind(pval2[vars], vars, sela[2]) )
    n.tests <- n.tests + length(vars)
    pval[vars] <- pmax(pval[vars], pval2[vars]) 
    ide <- which(pval[vars] < alpha)
    vars <- vars[ide]
    if ( max_k >= 2 ) {
      pval2 <- glmm.condregs(target = target, reps = reps, id = group, dataset = dataset, xIndex = vars, csIndex = sela, test = test, 
                             wei = wei, slopes = slopes, ncores = 1)$pvalue 
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
      pval2 <- glmm.condregs(target = target, reps = reps, id = group, dataset = dataset, xIndex = vars, csIndex = sela[3], test = test, 
                             wei = wei, slopes = slopes, ncores = 1)$pvalue 
      kapa_pval[[ 1 ]] <- cbind( kapa_pval[[ 1 ]], rbind(pval2[vars], vars, sela[3]) )
      n.tests <- n.tests + length(vars)
      pval[vars] <- pmax(pval[vars], pval2[vars]) 
      ide <- which(pval[vars] < alpha)
      vars <- vars[ide]
    }  ## end  if ( max_k >= 1 ) {
    if ( length(vars) > 0  &  max_k >= 2 ) {
      cand <- Rfast::comb_n(sela, 2)[, - 1]
      pval2 <- glmm.condregs(target = target, reps = reps, id = group, dataset = dataset, xIndex = vars, csIndex = cand[, 1], test = test, 
                             wei = wei, slopes = slopes, ncores = 1)$pvalue 
      kapa_pval[[ 2 ]] <- cbind( kapa_pval[[ 2 ]], rbind( pval2[vars], vars, matrix( rep(cand[, 1], length(vars)), ncol = length(vars) ) ) )
      n.tests <- n.tests + length(vars)
      pval[vars] <- pmax(pval[vars], pval2[vars]) 
      ide <- which(pval[vars] < alpha)
      vars <- vars[ide]
      pval2 <- glmm.condregs(target = target, reps = reps, id = group, dataset = dataset, xIndex = vars, csIndex = cand[, 2], test = test, 
                             wei = wei, slopes = slopes, ncores = 1)$pvalue 
      kapa_pval[[ 2 ]] <- cbind( kapa_pval[[ 2 ]], rbind( pval2[vars], vars, matrix( rep(cand[, 2], length(vars)), ncol = length(vars) ) ) )
      n.tests <- n.tests + length(vars)
      pval[vars] <- pmax(pval[vars], pval2[vars]) 
      ide <- which(pval[vars] < alpha)
      vars <- vars[ide]
    }  ## end  if ( max_k >= 2 ) {
    if ( length(vars) > 0  &  max_k >= 3 ) {
      pval2 <- glmm.condregs(target = target, reps = reps, id = group, dataset = dataset, xIndex = vars, csIndex = sela, test = test, 
                             wei = wei, slopes = slopes, ncores = 1)$pvalue 
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
    pval2 <- glmm.condregs(target = target, reps = reps, id = group, dataset = dataset, xIndex = vars, csIndex = sela[dm], test = test, 
                           wei = wei, slopes = slopes, ncores = 1)$pvalue 
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
        while ( length(vars) > 0  &  j < dim(cand)[2] ) {
          pval2 <- glmm.condregs(target = target, reps = reps, id = group, dataset = dataset, xIndex = vars, csIndex = cand[, j], test = test, 
                                 wei = wei, slopes = slopes, ncores = 1)$pvalue 
          kapa_pval[[ i ]] <- cbind( kapa_pval[[ i ]], rbind( pval2[vars], vars, matrix( rep(cand[, j], length(vars)), ncol = length(vars) ) ) )
          n.tests <- n.tests + length(vars)
          pval[vars] <- pmax(pval[vars], pval2[vars]) 
          ide <- which(pval[vars] < alpha)
          vars <- vars[ide]
          j <- j + 1
        }  ## end  for ( j in 1:dim(cand)[2] ) {
      }  ## end  for ( i in 2:max_k ) {  
    }  ## end  if ( max_k >= 2 ) {
    sel <- which.min(pval[vars])
    sela <- c(sela, vars[sel] )
    vars <- setdiff(vars, vars[sel])  
  } ## end  while ( length(vars) > 0  &  max_k >= 1 ) {
  
  runtime <- proc.time() - runtime
  
  # if ( backward  & length( sela ) > 0  ) {
  #   tic <- proc.time()
  #   pv <- pval[sela]
  #   sela <- sela[ order(pv) ]
  #   bc <- mmpcbackphase(target, dataset[, sela, drop = FALSE], test = test, wei = wei, max_k = max_k, threshold = threshold)
  #   met <- bc$met
  #   sela <- sela[met]
  #   pvalues[sela] <- bc$pvalue
  #   n.tests <- n.tests + bc$counter
  #   runtime <- runtime + proc.time() - tic
  # }  
  
  list(selectedVars = sela, pvalues = pval, univ = mod$pvalue, kapa_pval = kapa_pval, max_k = max_k, threshold = alpha, 
       n.tests = n.tests, runtime = runtime, test = ci_test, slope = slopes)
}





