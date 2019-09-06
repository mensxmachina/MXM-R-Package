iamb <- function(target, dataset, threshold = 0.05, wei = NULL, test = NULL, user_test = NULL, stopping = "BIC", tol = 2, ncores = 1, back = "iambbs") {
  
  mod <- fs.reg(target, dataset, ini = NULL, threshold, wei, test, user_test, stopping, tol, ncores) 
  poies <- mod$info[, 1]
  
  if ( length(poies) <= 1 )  {
    res <- list(vars = poies, mod = mod, mess = paste("No backward regression performed") ) 

  } else {
    
    if ( back == "iambbs" ) {
      mod2 <- iamb.bs(target, dataset[, poies, drop = FALSE], threshold, wei, test, user_test)
      sel <- mod2$mat[, 1]
    } else {
      mod2 <- bs.reg(target, dataset[, poies, drop = FALSE], threshold, wei, test, user_test)
      sel <- mod2$mat[, 1]
    }  
    poies <- poies[sel]
    mod2$mat[, 1] <- poies
    res <- list(vars = poies, mod = mod2)
    
  }  
    
  res
}