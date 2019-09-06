big.gomp <- function(target = NULL, dataset, tol = qchisq(0.95, 1) + log(dim(x)[1]), test = "testIndFisher", method = "ar2") {

  tic <- proc.time()
  oop <- options(warn = -1) 
  on.exit( options(oop) )  
  
  if ( is.null(target) ) {
    if (test == "censIndCR" | test == "censIndWR") {
      y <- Surv(dataset[, 1], dataset[, 2])
      x <- bigmemory::sub.big.matrix(dataset, firstCol = 3)
    } else {  
      y <- dataset[, 1]
      x <- bigmemory::sub.big.matrix(dataset, firstCol = 2)
    } 
  } else  y <- target 
  
  n <- dim(x)[1]
  phi <- NULL

  if ( test == "testIndFisher" ) {
  tic <- proc.time()
  ######### SSE
    if ( method == "sse" ) {
      rho <- Rfast::Var(y) * (n - 1)
      r <- cor(y, x[])
      sel <- which.max(abs(r))
      sela <- sel
      res <- .lm.fit(x[, sel, drop = FALSE], y)$residuals
      rho[2] <- sum(res^2)
      i <- 2
      while ( (rho[i - 1] - rho[i])/(rho[i - 1]) > tol  &  i < n ) {
        i <- i + 1
        r <- cor(res, x[])
        r[sela] <- NA
        sel <- which.max(abs(r))
        sela <- c(sela, sel)
        res <- .lm.fit(x[, sela], y)$residuals
        rho[i] <- sum(res^2)
        ind[sela] <- 0
      } ## end while ( (rho[i - 1] - rho[i])/(rho[i - 1]) > tol  &  i < n )
	len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Vars", "|sse|")
    ######### loglik
    } else if ( method == "loglik" ) {
      rho <- n * log( Rfast::Var(y) ) 
      r <- cor(y, x[])
      sel <- which.max( abs(r) )
      sela <- sel
      res <- .lm.fit(x[, sel, drop = FALSE], y)$residuals
      rho[2] <- n * log( sum(res^2)/(n - i) ) 
      i <- 2
      while ( rho[i - 1] - rho[i] > tol & i < n ) {
        i <- i + 1
        r <- cor(res, x[])
	  r[sela] <- NA
        sel <- which.max(abs(r))
        sela <- c(sela, sel)
        res <- .lm.fit(x[, sela], y)$residuals
        rho[i] <- n * log( sum(res^2)/(n - i) )
      } ## end while ( rho[i - 1] - rho[i] > tol & i < n )
	len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len] + n * log(2 * pi) + n )
      colnames(res) <- c("Vars", "BIC")
    ######### adjusted R-square
    } else if (method == "ar2") {
      down <- Rfast::Var(y) * (n - 1)
      rho <- 0
      r <- cor(y, x[])
      sel <- which.max( abs(r) )
      sela <- sel
      res <- .lm.fit(x[, sel, drop = FALSE], y)$residuals
      r2 <- 1 - sum(res^2)/down
      rho[2] <- 1 - (1 - r2) * (n - 1)/(n - 2)
	i <- 2
      while ( rho[i] - rho[i - 1] > tol & i < n ) {
        i <- i + 1
        r <- cor(res, x[])
        r[sela] <- NA
        sel <- which.max(abs(r))
        sela <- c(sela, sel)
        res <- .lm.fit(x[, sela], y)$residuals
        r2 <- 1 - sum(res^2)/down
        rho[i] <- 1 - (1 - r2) * (n - 1)/(n - i - 1)
      } ## end while ( rho[i] - rho[i - 1] > tol & i < n ) 
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Vars", "adjusted R2")
    }
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi, res = res)
      
  } else {
  
  if (test == "testIndGamma") {
    tic <- proc.time()
    mod <- glm( y ~ 1, family = Gamma(log) )
    rho <- mod$deviance
    phi <- summary(mod)[[ 14 ]]
    res <- y - fitted(mod)
    ela <- as.vector( cor(res, x[]) )
    sel <- which.max( abs(ela) )  
    sela <- sel
    names(sela) <- NULL
    mod <- glm(y ~ x[, sela], family = Gamma(log) )
    res <-  mod$residuals
    rho[2] <- mod$deviance
    phi[2] <- summary(mod)[[ 14 ]]
    i <- 2
    while ( (rho[i - 1] - rho[i]) / phi[i] > tol ) {
      i <- i + 1
	r <- cor(res, x[])
	r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- glm(y ~ x[, sela], family = Gamma(log) )
      res <- y - fitted(mod)
      rho[i] <- mod$deviance
      phi[i] <- summary(mod)[[ 14 ]]
    } ## end while ( (rho[i - 1] - rho[i]) > tol )
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi[1:len], res = res)
	
  } else if ( test == "testIndNormLog" ) {
    tic <- proc.time()
    ini <- Rfast::normlog.mle(y)
    m <- ini$param[1]
    rho <- sum( (y - ini$param[1])^2 )
    phi <- mod$devi/(n - 2)
    ela <- cor(y - m, x[])
    sel <- which.max( abs(ela) )
    sela <- sel
    names(sela) <- NULL 
    mod <- Rfast::normlog.reg(y, x[, sel])
    res <- y - exp( mod$be[1] + x[, sel] * mod$be[2] )
    rho[2] <- mod$deviance
    phi[2] <- mod$devi/(n - 2)
    i <- 2
    while ( (rho[i - 1] - rho[i]) / phi[i] > tol ) {
      i <- i + 1
	r <- cor(res, x[])
      r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- Rfast::normlog.reg(y, x[, sela])        
      res <- y - as.vector( exp( mod$be[1] + x[, sela] %*% mod$be[-1] ) ) 
      rho[i] <- mod$deviance
	phi[i] <- mod$deviance/(n - length(mod$be) )
    }
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi[1:len], res = res)

  } else if ( test == "testIndLogistic" ) {
    tic <- proc.time()
    n <- dim(x)[1]
    p <- sum(y)/n
    rho <-  -2 * (n * p * log(p) + (n - n * p) * log(1 - p)) 
    ela <- as.vector( cor(y - p, x[]) )
    sel <- which.max( abs(ela) )
    sela <- sel
    names(sela) <- NULL
    mod <- Rfast::glm_logistic(x[, sel], y)
    est <- exp(-mod$be[1] - x[, sel] * mod$be[2])
    res <- y - 1/(1 + est)
    rho[2] <- mod$devi
    i <- 2
    while ( (rho[i - 1] - rho[i]) > tol ) {
      i <- i + 1
	r <- cor(res, x[])
	r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- Rfast::glm_logistic(x[, sela], y)
      est <- as.vector(exp(-mod$be[1] - x[, sela] %*% mod$be[-1]))
      res <- y - 1/(1 + est)
      rho[i] <- mod$devi
    }
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi, res = res)

  } else if (test == "testIndQBinom") {
    tic <- proc.time()
    p <- sum(y)/n
    y0 <- 1 - y
    rho <- 2 * sum(y * log(y/p), na.rm = TRUE) + 2 * sum(y0 * log(y0/(1 - p)), na.rm = TRUE)    
    phi <- 1
    ela <- as.vector( cor(y - p, x[]) )
    sel <- which.max(abs(ela))
    sela <- sel
    names(sela) <- NULL
    mod <- Rfast::prop.reg(y, x[, sel], varb = "glm")
    est <- exp(-mod$info[1, 1] - x[, sel] * mod$info[2, 1])
    p <- 1/(1 + est)
    res <- y - p
    rho[2] <- 2 * sum(y * log(y/p), na.rm = TRUE) + 2 * sum(y0 * log(y0/(1 - p)), na.rm = TRUE)
    phi[2] <- mod$phi
    i <- 2
    while ((rho[i - 1] - rho[i])/phi[i] > tol) {
      i <- i + 1
     	r <- cor(res, x[])
	    r[sela] <- NA
      sel <- which.max(abs(r))
      sela <- c(sela, sel)
      mod <- Rfast::prop.reg(y, x[, sela], varb = "glm")
      est <- as.vector( exp(-mod$info[1, 1] - x[, sela] %*% mod$info[-1, 1]) )
      p <- 1/(1 + est)
      res <- y - p
      rho[i] <- 2 * sum(y * log(y/p), na.rm = TRUE) + 2 * sum(y0 * log(y0/(1 - p)), na.rm = TRUE)
      phi[i] <- mod$phi
    }
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi[1:len], res = res)
	
  } else if ( test == "testIndPois" ) {
    tic <- proc.time()
    m <- sum(y)/n
    rho <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * log(m)
    ela <- as.vector( cor(y - m, x[]) )
    sel <- which.max( abs(ela) )
    sela <- sel
    names(sela) <- NULL
    mod <- Rfast::glm_poisson(x[, sel], y)
    res <- y - exp( mod$be[1] + x[, sel] * mod$be[2] )
    rho[2] <- mod$devi
    i <- 2
    while ( (rho[i - 1] - rho[i]) > tol ) {
      i <- i + 1
	    r <- cor(res, x[])
	    r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- Rfast::glm_poisson(x[, sela], y)        
      res <- y - as.vector( exp( mod$be[1] + x[, sela] %*% mod$be[-1] ) ) 
      rho[i] <- mod$devi
    }
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi, res = res)
	
  } else if ( test == "testIndQPois" ) {
    tic <- proc.time()
    m <- sum(y)/n
    rho <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * log(m)
    phi <- 1
    ela <- as.vector( cor(y - m, x[]) )
    sel <- which.max( abs(ela) )
    sela <- sel
    names(sela) <- NULL
    mod <- Rfast::qpois.reg(x[, sel], y)
    phi[2] <- mod$phi
    res <- y - exp( mod$be[1, 1] + x[, sel] * mod$be[2, 1] )
    rho[2] <- mod$devi
    i <- 2
    while ( (rho[i - 1] - rho[i]) / phi[i] > tol ) {
      i <- i + 1
	    r <- cor(res, x[])
    	r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- Rfast::qpois.reg(x[, sela], y)        
      res <- y - as.vector( exp( mod$be[1, 1] + x[, sela] %*% mod$be[-1, 1] ) ) 
      rho[i] <- mod$devi
      phi[i] <- mod$phi 
    }
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi[1:len], res = res)
	
  } else if (test == "testIndNB") {
    tic <- proc.time()
    mod <- MASS::glm.nb(y ~ 1)
    rho <-  - 2 * as.numeric( logLik(mod) ) 
    res <-  y - fitted(mod)
    ela <- as.vector( cor(res, x[]) )
    sel <- which.max( abs(ela) )  
    sela <- sel
    names(sela) <- NULL
    mod <- MASS::glm.nb(y ~ x[, sela], control = list(epsilon = 1e-08, maxit = 100, trace = FALSE) )
    res <-  mod$residuals
    rho[2] <-  - 2 * as.numeric( logLik(mod) )
    i <- 2
    while ( (rho[i - 1] - rho[i]) > tol ) {
      i <- i + 1
	    r <- cor(res, x[])
	    r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- MASS::glm.nb( y ~ x[, sela], control = list(epsilon = 1e-08, maxit = 100, trace = FALSE) )        
      res <- y - fitted(mod)
      rho[i] <-  - 2 * as.numeric( logLik(mod) )
    } ## end while ( (rho[i - 1] - rho[i]) > tol )
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi, res = res)
           
  } else if ( test == "testIndMMReg") {
    tic <- proc.time()
    mod <- MASS::rlm(y ~ 1, method = "MM", maxit = 2000)
    rho <-  - 2 * as.numeric( logLik(mod) ) 
    res <- mod$residuals
    ela <- as.vector( cor(res, x[]) )
    sel <- which.max( abs(ela) )  
    sela <- sel
    names(sela) <- NULL
    mod <- MASS::rlm(y ~ x[, sela], method = "MM", maxit = 2000 )
    res <- mod$residuals
    rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
    i <- 2
    while ( (rho[i - 1] - rho[i]) > tol ) {
      i <- i + 1
    	r <- cor(res, x[])
    	r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- try( MASS::rlm(y ~ x[, sela], method = "MM", maxit = 2000 ), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        rho[i] <- rho[i - 1]
      } else {
        res <- mod$residuals
        rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
      }  
    } ## end while ( (rho[i - 1] - rho[i]) > tol )
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi, res = res)
      
  } else if ( test == "testIndRQ") {
    tic <- proc.time()
    mod <- quantreg::rq(y ~ 1)
    rho <-  - 2 * as.numeric( logLik(mod) ) 
    res <- mod$residuals
    ela <- as.vector( cor(res, x[]) )
    sel <- which.max( abs(ela) )  
    sela <- sel
    names(sela) <- NULL
    mod <- quantreg::rq(y ~ x[, sela])
    res <- mod$residuals
    rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
    i <- 2
    while ( (rho[i - 1] - rho[i]) > tol ) {
      i <- i + 1
      r <- cor(res, x[])
    	r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- try( quantreg::rq(y ~ x[, sela]), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        rho[i] <- rho[i - 1]
      } else {
        res <- mod$residuals
        rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
      }  
    } ## end while ( (rho[i - 1] - rho[i]) > tol )
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi, res = res)
      
  } else if (test == "testIndOrdinal") {
    tic <- proc.time()
    mod <- MASS::polr(y ~ 1)
    rho <-  - 2 * as.numeric( logLik(mod) )
    res <- ord.resid(y, mod$fitted.values)
    ela <- as.vector( cor(res, x[]) )
    sel <- which.max( abs(ela) )     
    sela <- sel
    names(sela) <- NULL
    mod <- MASS::polr(y ~ x[, sela])
    res <- ord.resid(y, mod$fitted.values) 
    rho[2] <-  - 2 * as.numeric( logLik(mod) )
    i <- 2
    while ( (rho[i - 1] - rho[i]) > tol ) {
      i <- i + 1
      r <- cor(res, x[])
	    r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- MASS::polr(y ~ x[, sela])        
      res <- ord.resid(y, mod$fitted.values)
      rho[i] <-  - 2 * as.numeric( logLik(mod) )
    } ## end while ( (rho[i - 1] - rho[i]) > tol )
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi, res = res)
      
  } else if ( test == "testIndTobit") {
    tic <- proc.time()
    mod <- survival::survreg(y ~ 1, dist = "gaussian")
    rho <-  - 2 * as.numeric( logLik(mod) ) 
    res <- resid(mod)
    ela <- as.vector( cor(res, x[]) )
    sel <- which.max( abs(ela) )  
    sela <- sel
    names(sela) <- NULL
    mod <- survival::survreg(y ~ x[, sela], dist = "gaussian" )
    res <- resid(mod)
    rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
    i <- 2
    while ( (rho[i - 1] - rho[i]) > tol ) {
      i <- i + 1
    	r <- cor(res, x[])
	    r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- try( survival::survreg(y ~ x[, sela], dist = "gaussian" ), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        rho[i] <- rho[i - 1]
      } else {
        res <- resid(mod)
        rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
        r[sela] <- NA
      }  
    } ## end while ( (rho[i - 1] - rho[i]) > tol )
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi, res = res)

  } else if ( test == "censIndCR") {
    tic <- proc.time()
    mod <- survival::coxph(y ~ 1)
    rho <-  - 2 * summary( mod)[[1]]
    res <- mod$residuals   ## martingale residuals
    ela <- as.vector( cor(res, x[]) )
    sel <- which.max( abs(ela) )  
    sela <- sel
    names(sela) <- NULL
    mod <- survival::coxph(y ~ x[, sela] )
    res <-  mod$residuals
    rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
    i <- 2
    while ( (rho[i - 1] - rho[i]) > tol ) {
      i <- i + 1
      r <- cor(res, x[])
	    r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- try( survival::coxph(y ~ x[, sela] ), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        rho[i] <- rho[i - 1]
      } else {
        res <- mod$residuals   ## martingale residuals
        rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
      }  
    } ## end while ( (rho[i - 1] - rho[i]) > tol )
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi, res = res)
    
  } else if ( test == "censIndWR") {
    tic <- proc.time()
    mod <- survival::survreg(y ~ 1)
    rho <-  - 2 * as.numeric( logLik(mod) ) 
    res <- resid(mod)
    ela <- as.vector( cor(res, x[]) )
    sel <- which.max( abs(ela) )  
    sela <- sel
    names(sela) <- NULL
    mod <- survival::survreg(y ~ x[, sela], control = list(iter.max = 10000) )
    res <- resid(mod)
    rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
    if ( is.na(rho[2]) )  rho[2] <- rho[1]
    i <- 2
    while ( (rho[i - 1] - rho[i]) > tol ) {
      i <- i + 1
      r <- cor(res, x[])
	    r[sela] <- NA
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- try( survival::survreg(y ~ x[, sela], control = list(iter.max = 10000) ), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        rho[i] <- rho[i - 1]
      } else {  
        res <- resid(mod)
        rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
      }  
    } ## end while ( (rho[i - 1] - rho[i]) > tol )
    len <- length(sela)
    res <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(res) <- c("Selected Vars", "Deviance") 
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = phi, res = res)

  } ##  end if (test == "testIndNB")
  
  }  ## end  if else (test == "testIndFisher") 
  
  result
}
