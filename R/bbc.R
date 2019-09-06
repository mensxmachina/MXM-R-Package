bbc <- function(predictions, target, metric = "auc.mxm", conf = 0.95, B = 1000) {
  
  dm <- dim(predictions)
  n <- dm[1]    ;    p <- dm[2]
  out.perf <- numeric(B)
  
  if (metric == "auc.mxm") {
    target <- as.numeric( as.factor(target) ) - 1
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      in.perf <- Rfast::colaucs(target[ind], predictions[ind, ])
      poio <- which.max(in.perf)
      out.perf[i] <- Rfast::auc(target[-ind], predictions[-ind, poio])
    } 
    
  } else if (metric == "fscore.mxm") {
    target <- as.numeric( as.factor(target) ) - 1
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      for (j in 1:p) {
        tab <- table(target[ind], predictions[ind, j])
        prec <- tab[2, 2]/(tab[2, 2] + tab[1, 2])
        rec <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
        in.perf[j] <- prec * rec / (prec + rec)    
      }
      poio <- which.max(in.perf)
      tab <- table(target[-ind], predictions[-ind, poio])
      prec <- tab[2, 2]/(tab[2, 2] + tab[1, 2])
      rec <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
      out.perf[i] <- 2 * prec * rec / (prec + rec)
    }

  } else if (metric == "euclid_sens.spec.mxm") {
    target <- as.numeric( as.factor(target) ) - 1
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      for (j in 1:p) {
        tab <- table(target[ind], predictions[ind, j])
        spec <- tab[1, 1]/(tab[1, 1] + tab[1, 2])
        sens <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
        in.perf[j] <- sqrt( (1 - sens)^2 + (1 - spec)^2 )
      }
      poio <- which.max(in.perf)
      tab <- table(target[-ind], predictions[-ind, poio])
      spec <- tab[1, 1]/(tab[1, 1] + tab[1, 2])
      sens <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
      out.perf[i] <- sqrt( (1 - sens)^2 + (1 - spec)^2 )
    }
    
  } else if (metric == "acc.mxm") {
    target <- as.numeric( as.factor(target) ) - 1
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      in.perf <- 1 - Rfast::colmeans( predictions[ind, ] - target[ind] )
      poio <- which.max(in.perf)
      out.perf[i] <- 1 - mean( predictions[-ind, poio] - target[-ind] )
    } 
  
  } else if (metric == "acc_multinom.mxm") {
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      in.perf <- 1 - Rfast::colmeans( predictions[ind, ] - target[ind] )
      poio <- which.max(in.perf)
      out.perf[i] <- 1 - mean( predictions[-ind, poio] - target[-ind] )
    } 
  
  } else if (metric == "mse.mxm") {
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      in.perf <-  - Rfast::colmeans( (predictions[ind, ] - target[ind])^2 ) 
      poio <- which.max(in.perf)
      out.perf[i] <-  - mean( (predictions[-ind, poio] - target[-ind])^2 ) 
    } 
    
  } else if (metric == "pve.mxm") {
    co <- length(target) - 1
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      in.perf <- 1 - Rfast::colsums( (predictions[ind, ] - target[ind])^2 ) / ( co * Rfast::Var(target[ind]) )
      poio <- which.max(in.perf)
      out.perf[i] <-  1 - sum( (predictions[-ind, poio] - target[-ind])^2 ) / ( co * Rfast::Var(target[-ind]) )
    } 
    
  } else if (metric == "ord_mae.mxm") {
    target <- as.numeric(target)
  	for (i in 1:dim(predictions)[2] ) 	predictions[, i] <- as.numeric(predictions[, i])
	  predictions <- as.matrix(predictions)
	  
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      in.perf <-  - Rfast:: colmeans( abs(predictions[ind, ] - target[ind]) ) 
      poio <- which.max(in.perf)
      out.perf[i] <-  - mean( abs(predictions[-ind, ] - target[-ind]) ) 
    } 
    
  } else if (metric == "mae.mxm") {
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      in.perf <-  - Rfast:: colmeans( abs(predictions[ind, ] - target[ind]) ) 
      poio <- which.max(in.perf)
      out.perf[i] <-  - mean( abs(predictions[-ind, ] - target[-ind]) ) 
    } 
    
  } else if (metric == "ci.mxm") {
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      for (j in 1:p)   in.perf[j] <- survival::survConcordance(target[ind] ~ predictions[ind, j])$concordance
      poio <- which.max(in.perf)
      out.perf[i] <- survival::survConcordance(target[-ind] ~ predictions[-ind, poio])$concordance
    }
    
  } else if (metric == "ciwr.mxm") {
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      for (j in 1:p)   in.perf[j] <- 1 - survival::survConcordance(target[ind] ~ predictions[ind, j])$concordance
      poio <- which.max(in.perf)
      out.perf[i] <- 1 - survival::survConcordance(target[-ind] ~ predictions[-ind, poio])$concordance
    }
    
  } else if (metric == "poisdev.mxm") {
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      in.perf <-  - colSums( target[ind] * log(target[ind] / predictions[ind, ]), na.rm = TRUE ) 
      poio <- which.max(in.perf)
      out.perf[i] <-  - 2 * sum( target[-ind] * log(target[-ind] / predictions[-ind, poio]), na.rm = TRUE ) 
    }
    
  } else if (metric == "binomdev.mxm") {
    in.perf <- numeric(p)
    ya <- target[, 1]     ;    N <- target[, 2]
    yb <- N - ya
    esta <- predictions     ;    estb <- N - esta
    for (i in 1:B) {
      ind <- sample.int(n, n, replace = TRUE)
      in.perf <-  - colSums( ya[ind] * log(ya[ind] / esta[ind, ]), na.rm = TRUE ) - colSums( yb[ind] * log(yb[ind] / estb[ind, ]), na.rm = TRUE ) 
      poio <- which.max(in.perf)
      out.perf[i] <-  - 2 * sum( ya[-ind] * log(ya[-ind] / esta[-ind, poio]), na.rm = TRUE ) - 2 * sum( yb[-ind] * log(yb[-ind] / estb[-ind, poio]), na.rm = TRUE ) 
    }
    
  } ## end  if (metric == " ") {
  a <- 0.5 * (1 - conf )
  ci <- quantile(out.perf, probs = c(a, 1 - a) )  
  list(bbc.perf = mean(out.perf), out.perf = out.perf, ci = ci )
  
}  

