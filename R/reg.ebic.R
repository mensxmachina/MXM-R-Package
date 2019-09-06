reg.ebic <- function(target, dataset, sela, test = "testIndLogistic", gam = NULL, wei = NULL) {
  
  dm <- dim(dataset)
  n <- dm[1]
  p <- dm[2]
  logn <- log(n)
  
  if ( is.null(gam) ) {
    con <- 2 - log(p) / logn
  } else con <- 2 * gam
  if ( con < 0 )  con <- 0
  M <- length(sela)
  
  if (test == "testIndReg") {
    fit2 <- lm( target ~., data = dataset[, sela], weights = wei )
    ebic <- BIC(fit2) + con * lchoose(p, M)
    
  } else if (test == "testIndLogistic") {
    fit2 <- glm( target ~ dataset[, sela], family = binomial, weights = wei, model = FALSE )
    ebic <- BIC(fit2) + con * lchoose(p, M)

  } else if (test == "testIndMultinom") {  
    fit2 <- nnet::multinom( target ~ dataset[, sela], weights = wei, trace = FALSE )
    ebic <- BIC(fit2) + con * lchoose(p, M)
    
  } else if (test == "testIndPois") {
    fit2 <- glm( target ~ dataset[, sela], family = poisson, weights = wei, model = FALSE )
    ebic <- BIC(fit2) + con * lchoose(p, M)

  } else if (test == "testIndNB") {
    fit2 <- MASS::glm.nb( target ~ dataset[, sela], weights = wei, model = FALSE )
    lik2 <- BIC(fit2) + con * lchoose(p, M)
    
  } else if (test == "testIndGamma") {
    fit2 <- glm( target ~ dataset[, sela], family = Gamma(log), weights = wei, model = FALSE )
    ebic <- BIC(fit2) + con * lchoose(p, M)
    
  } else if (test == "testIndNormLog") {
    fit2 <- glm( target ~ dataset[, sela], family = gaussian(log), weights = wei, model = FALSE )
    ebic <- BIC(fit2) + con * lchoose(p, M)
    
  } else if (test == "testIndBinom") {
    y <- target/wei
    fit2 <- glm( target ~ dataset[, sela], family = binomial, weights = wei, model = FALSE )
    ebic <- BIC(fit2) + con * lchoose(p, M)
  
  } else if (test == "testIndClogit") {  
    case <- target[, 1]
    id <- target[, 2]
    fit2 <- try( survival::clogit( case ~ . + strata(id), data = dataset[, sela] ), silent = TRUE)
    if ( identical( class(fit2), "try-error" ) ) {
      ebic <- NULL
    } else  ebic <- BIC(fit2) + con * lchoose(p, M)

  } else if (test == "censIndCR") {  
    fit2 <- try( survival::coxph( target ~., data = dataset[, sela], weights = wei ), silent = TRUE )
    if ( identical( class(fit2), "try-error" )  ) {
      ebic <- NULL
    } else  ebic <- BIC(fit2) + con * lchoose(p, M)

  } else if (test == "censIndWR") {  
    fit2 <- try( survival::survreg( target ~., data = dataset[, sela], weights = wei, control = list(iter.max = 10000) ), silent = TRUE )
    if ( identical( class(fit2), "try-error" ) ) {
      ebic <- NULL
    } else  ebic <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn +  con * lchoose(p, M)
    
  } else if (test == "testIndMMReg") {
    fit2 <- MASS::rlm( target ~ dataset[, sela], method = "MM", maxit = 2000 )
    ebic <- BIC(fit2) + con * lchoose(p, M)
   
  } else if (test == "testIndOrdinal") {
    fit2 <- ordinal::clm( target ~ dataset[, sela], weights = wei, model = FALSE ) 
    ebic <- BIC(fit2) + con * log(p)

  } else if (test == "testIndTobit") {
    fit2 <- try( survival::survreg( target ~ dataset[, sela], weights = wei, dist = "gaussian" ), silent = TRUE )
    if ( identical( class(fit2), "try-error" ) ) {
      ebic <- NULL
    } else  ebic <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn +  con * lchoose(p, M)
    
  }
  
  
  
  
}