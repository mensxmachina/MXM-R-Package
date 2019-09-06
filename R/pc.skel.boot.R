pc.skel.boot <- function(dataset, method = "pearson", alpha = 0.01, R = 199, ncores = 1) {
  ## dataset contains the data, it must be a matrix 
  ## type can be either "pearson" or "spearman" for continuous variables OR
  ## "cat" for categorical variables
  ## alpha is the level of significance, set to 0.01 by default
  ## ncores is for parallel computations
  G <- Rfast::pc.skel(dataset = dataset, method = method, alpha = alpha, R = 1)$G
  title <- deparse( substitute(dataset) )
  dm <- dim(dataset)
  n <- dm[1]
  p <- dm[2]
  
  if (ncores <= 1) {  ## one core
    gboot <- matrix(0, nrow = R, ncol = p^2)
    for (i in 1:R) {
      id <- sample(n, n, replace = TRUE)
      gb <- Rfast::pc.skel(dataset = dataset[id, ], method = method, alpha = alpha, R = 1)$G
      gboot[i, ] <- as.vector(gb)
    }  ## end for (i in 1:R)

  } else {  ## parallel computations
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    gboot <- foreach( i = 1:R, .combine = rbind, .export = "pc.skel", .packages = "Rfast" ) %dopar% {
      id <- sample(n, n, replace = TRUE)
      gb <- Rfast::pc.skel(dataset = dataset[id, ], method = method, alpha = alpha, R = 1)$G
      return( as.vector(gb) )
    }
    stopCluster(cl)
  }  ## end if (ncores <= 1)
  Gboot <- Rfast::colmeans(gboot)
  Gboot <- matrix(Gboot, nrow = p, ncol = p)
  list(G = G, Gboot = Gboot, title = title)
}    