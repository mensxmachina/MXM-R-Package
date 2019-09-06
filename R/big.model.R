big.model <- function(y, x, test) {
  
  if ( test == "testIndLogistic" ) {
    if ( is.null(x) ) {
      n <- length(y)
      p <- sum(y)/n
      devi <- -2 * ( n * p * log(p) + (n - n * p) * log(1 - p) )
      res <- list(be = log( p/(1 - p) ), devi = devi, phi = NA)
    } else {
      mod <- Rfast::glm_logistic(x, y, maxiters = 5000)
      res <- list(be = mod$be, devi = mod$devi, phi = NA)
    }  
    
  } else if ( test == "testIndMultinom" ) {
    if ( is.null(x) ) {
      mod <- Rfast::multinom.mle( Rfast::design_matrix(y, ones = FALSE) )
      p <- mod$prob
      res <- list(be = log(p[-1] / p[1] ), devi = - 2 * mod$loglik, phi = NA)
    } else {
      mod <- try( Rfast::multinom.reg(x, y, maxiters = 10000), silent = TRUE )
      if ( identical(class, "try-error") ) {
        res <- list(be = NULL, devi = NULL, phi = NA)
      } else   res <- list(be = mod$be, devi = mod$devi, phi = NA)
    }
    
  } else if ( test == "testIndPois" ) {
    if ( is.null(x) ) {
      n <- length(y)
      m <- sum(y)/n
      devi <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * log(m)
      res <- list(be = log(m), devi = devi, phi = NA)
    } else {
      mod <- Rfast::glm_poisson(x, y)
      res <- list(be = mod$be, devi = mod$devi, phi = NA)
    }
    
  } else if ( test == "testIndQPois" ) {
    if ( is.null(x) ) {
      n <- length(y)
      m <- sum(y)/n
      devi <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * log(m)
      res <- list(be = log(m), devi = devi, phi = NA)
    } else {
      mod <- Rfast::qpois.reg(x, y, maxiters = 5000)
      res <- list(be = mod$be, devi = mod$devi, phi = mod$phi)
    }
    
  } else if ( test == "censIndWR" ) {
    if ( is.null(x) ) {
      mod <- Rfast::weibull.mle(y, maxiters = 500)
      res <- list(be = log(mod$param[1]), devi = 2 * mod$loglik, phi = NA)
    } else {
      mod <- survival::survreg(y ~ x, control = list(iter.max = 10000) )
      res <- list(be = mod$coefficients, devi = 2 * mod$loglik[2], phi = NA)
    }
    
  } else if ( test == "testIndFisher" ) {
    mod <- try( Rfast::lmfit( cbind(1, x), y ), silent = TRUE )
    if ( identical(class, "try-error") ) {
      res <- list(be = NULL, devi = NULL, phi = NA)
    } else {
      dm <- dim(x)
      devi <- sum( (mod$residuals)^2 )
      res <- list(be = mod$be, devi = devi, phi = devi/(dm[1] - dm[2] - 1) )
    }  
  }
  
  res
  
}
