####################
#### Univariate ridge regression
####################
ridge.reg <- function(target, dataset, lambda, B = 1, newdata = NULL) {
  ## target is the dependent variable and can be a matrix as well
  ## However we only use it with a univariate target variable
  ## dataset contains the independent, CONTINUOUS ONLY, variables
  ## lambda is the ridge regularization parameter
  ## if lambda=0, the classical multivariate regression is implemented
  ## B is for bootstrap estimation of the standard errors of the betas
  ## newdata is the new independent variables values 
  ## whose values of y you want to estimate
  ## by default newdata is NULL
  target <- as.vector(target)
  dataset <- as.matrix(dataset)
  dm <- dim(dataset)
  n <- dm[1]  ## sample size
  p <- dm[2]  ## dimensionality of dataset
  my <- sum(target) / n
  yy <- target - my  ## center the dependent variables
  xtx <- crossprod(dataset)
  lamip <- lambda * diag(p)
  
  W <- solve( xtx + lamip )
  betas <- W %*% crossprod(dataset, yy) 
  est <- dataset %*% betas + my
  va <- Rfast::Var(target - est) * (n - 1) / (n - p - 1)
  # vab <- kronecker(va, W %*% xtx %*% W  ) 
  # seb <- as.vector( sqrt( diag(vab) ) )
  vab <- va * mahalanobis(W, numeric(p), xtx, inverted = TRUE)
  seb <- sqrt( vab )

  if (B > 1) { ## bootstrap estimation of the standard errors
    be <- matrix(nrow = B, ncol = p )
    for ( i in 1:B) {
       id <- sample(1:n, n, replace = TRUE)
       yb <- yy[id]     ;     xb <- dataset[id, ]
       be[i, ] <- solve( crossprod(xb) + lamip, crossprod(xb, yb) )
    }
    seb <- Rfast::colVars(be, std = TRUE) ## bootstrap standard errors of betas
  } 
  
  be <- as.vector(betas)
  ## seb contains the standard errors of the coefficients
  if ( is.null( colnames(dataset) ) ) {
    names(seb) <- paste("X", 1:p, sep = "")
    names(be) <- paste("X", 1:p, sep = "")
  } else  names(seb) <- names(be) <- colnames(dataset)
  
  est <- NULL  
  if ( !is.null(newdata) ) {  
    newdata <- as.matrix(newdata)
    newdata <- matrix(newdata, ncol = p)
    est <- as.vector( newdata %*% beta ) + my 
  } 
  
  list(beta = be, seb = seb, est = est)
}


