kfbed.gee.reg <- function(y, x, id, reps = NULL, univ = NULL, alpha = 0.05, wei = NULL, K = 0:5, test = "testIndGEEReg", correl = "exchangeable", se = "jack") { 
    
  a <- fbed.gee.reg(target = y, dataset = x, id = id, reps = reps, ini = univ, threshold = alpha, wei = wei, K = max(K), test = test, correl = correl, se = se) 
    
  info <- a$info
  k <- dim(info)[1]
  mod <- list()
  if ( k > 0  & info[1, 1] > 0 ) {
    if (info[k, 1] == info[k - 1, 1])   k <- k - 1
    sel <- info[1:k, 1]
    for (i in 1:k)   mod[[ i ]] <- a$res[1:sel[i], ]   ## end if (backward) 
    names(mod) <- paste("K=", 0:(k-1), sep = "" )
  }  ## end if (k > 0)
  
  list(res = a, mod = mod)
}
