kfbed.glmm.reg <- function(y, x, id, univ = NULL, alpha = 0.05, wei = NULL, K = 0:5, method = "LR", gam = NULL, backward = TRUE, test = "testIndGLMMReg") {
    
  a <- fbed.glmm.reg(target = y, dataset = x, id = id, ini = univ, threshold = alpha, wei = wei, K = max(K), method = method, gam = gam, backward = FALSE, test = test)
    
  info <- a$info
  res <- a$res
  k <- dim(info)[1]
  mod <- list()
  if ( k > 0  & info[1, 1] > 0 ) {
    if (info[k, 1] == info[k - 1, 1])   k <- k - 1
    sel <- info[1:k, 1]
    if (backward) {
      if (method == "LR") {
        
        for (i in 1:k) {
          b <- glmm.bsreg(y, x[, a$res[1:sel[i], 1], drop = FALSE], id = id, threshold = alpha, wei = wei, test = test)
          
          if ( typeof(b) == "list" ) {
            mod[[ i ]] <- cbind(a$res[b$mat[, 1], 1], b$mat[, 3], b$mat[, 2])
            colnames(mod[[ i ]]) <- c("Vars", "stat", "log p-value")
          } else   mod[[ i ]] <- NULL
        }  ## end for
      } else {
        
        for (i in 1:k) {
          b <- ebic.bsreg(y, x[, res[1:sel[i], 1], drop = FALSE], test = test, wei = wei, gam = gam) 
          b <- ebic.glmm.bsreg(y, x[, res[1:sel[i], 1], drop = FALSE], id = id, wei = wei, gam = gam, test = test) 
          if ( typeof(b) == "list" ) {
            mod[[ i ]] <- b$mat
          } else  mod[[ i ]] <- NULL
        }      
        
      }  ## end if (method == "LR")
      
    }  else     for (i in 1:k)   mod[[ i ]] <- a$res[1:sel[i], ]   ## end if (backward) 
    names(mod) <- paste("K=", 0:(k-1), sep = "" )
  }  ## end if (k > 0)

  list(res = a, mod = mod)
}
