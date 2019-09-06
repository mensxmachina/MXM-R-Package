rint.regs <- function(target, dataset, targetID = -1, id, reps = NULL, tol = 1e-08) {
  
  univariateModels <- list();
  dm <- dim(dataset)
  n <- dm[1]
  D <- dm[2]
  ind <- 1:D
  
  if (targetID != -1 ) {
    target <- dataset[, targetID]
    dataset[, targetID] <- rnorm(n)
  }   
  poia <- NULL
  poia <- Rfast::check_data(dataset)
  if ( sum(poia > 0) )  ind[poia] <- 0
  
  if ( is.null(reps) ) {
    mod <- Rfast::rint.regs(target, dataset, id, logged = TRUE)
    mod[ is.na(mod) ] <- 0
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]  ## pf(stat, 1, n - 4, lower.tail = FALSE, log.p = TRUE)
    
  } else {
    
    sy <- as.vector( rowsum(target, id) )
    ni <- tabulate(id)
    my <- sy / ni
    ni2 <- ni^2
    Xi <- cbind(1, reps, dataset[, 1])
    sxall <- rowsum(cbind(1, reps, dataset), id)
    
    funa <- function(d, n, ni, ni2, S, hi2)   sum( log1p(ni * d) ) + n * log(S - d * sum(ni2 * hi2/ (1 + ni * d) ) )    
    stat <- numeric(D)
    for (i in ind) {
      Xi[, 3] <- dataset[, i]
      xx <- crossprod(Xi)
      # sx <- rowsum(Xi, id)
      sx <- sxall[, c(1, 2, i + 2)]
      sxy <- crossprod(Xi, target)
      mx <- sx / ni
      #############
      # b1 <- solve(xx, sxy)  ## not much of a difference!!
      # S <- sum( (target - Xi %*% b1)^2 )  ## not much of a difference!!
      mod <- .lm.fit(Xi, target)
      b1 <- mod$coefficients
      S <- sum( mod$residuals^2 )
      hi2 <- ( my - mx %*% b1 )^2
      mod <- optimise(funa, c(0, 70), n = n, ni = ni, ni2 = ni2, S = S, hi2 = hi2, tol = tol)
      d <- mod$minimum 
      b2 <- solve( xx - d * crossprod(sx/(1+ ni * d), sx), sxy - d * crossprod( sx, sy/(1 + ni * d) ) )   
      k <- 2
      while ( sum( abs(b2 - b1) ) > tol  &  k < 100) {
        k <- k + 1
        b1 <- b2
        S <- sum( (target - Xi %*% b1)^2 )
        hi2 <- ( my - mx %*% b1 )^2
        mod <- optimise(funa, c(0, 70), n = n, ni = ni, ni2 = ni2, S = S, hi2 = hi2, tol = tol)
        d <- mod$minimum 
        b2 <- solve( xx - d * crossprod(sx/(1+ ni * d), sx), sxy - d * crossprod( sx, sy/(1 + ni * d) ) ) 
      }
      se <- (S - d * sum(ni^2 * hi2/ (1 + ni * d) ) )/n
      seb <- solve( xx - d * crossprod(sx/(1+ ni * d), sx) )[3, 3] * se 
      stat[i] <- b2[3]^2 / seb
    }
    
    univariateModels$stat <- stat
    univariateModels$pvalue <- pf(stat, 1, n - 5, lower.tail = FALSE, log.p = TRUE)
    
  }  ## else if ( is.null(reps) )

  if (targetID != - 1) {
    univariateModels$stat[targetID] <- 0
    univariateModels$pvalue[targetID] <- log(1)
  }
  
  univariateModels
}   



















# rint.regs_old <- function(target, dataset, targetID = -1, id, reps = NULL, tol = 1e-08) {
#   
#   univariateModels <- list();
#   dm <- dim(dataset)
#   n <- dm[1]
#   D <- dm[2]
#   
#   if (targetID != -1 ) {
#     target <- dataset[, targetID]
#     dataset[, targetID] <- rnorm(n)
#   }   
#   poia <- NULL
#   poia <- Rfast::check_data(dataset)
#   if ( sum(poia > 0) )  dataset[, poia] <- rnorm(n * length(poia) )
# 
#   if ( is.null(reps) ) {
#     # ni <- tabulate(id)
#     # ni2 <- ni^2
#     # funa <- function(d, n, ni, ni2, S, hi2)   sum( log1p(ni * d) ) + n * log(S - d * sum(ni2 * hi2/ (1 + ni * d) ) )    
#     # sy <- as.vector( rowsum(target, id) )
#     # my <- sy / ni
#     # Sy <- sum(sy)
#     # r <- as.vector( cov(target, dataset) )
#     # mesi <- Sy/n
#     # xs <- colsums(dataset)
#     # xs2 <- colsums(dataset^2)
#     # vx <- (xs2 - xs^2/n)/(n - 1)
#     # b <- r/vx
#     # a <- mesi - b * xs/n
#     # be <- cbind(a, b)   
#     # stat <- numeric(D)
#     # sx <- cbind(ni, ni)
#     # xx <- matrix( c(n, 0, 0, 0), 2, 2)
#     # sxy <- c(Sy, Sy)
#     # ##############
#     # for (i in 1:D) {
#       # Xi <- dataset[, i]
#       # xx[2:3] <- xs[i]
#       # xx[4] <- xs2[i] 
#       # sx[, 2] <- rowsum(Xi, id)
#       # sxy[2] <- sum(Xi * target)
#       # mx <- sx[, 2]/ni
#       # #############
#       # b1 <- be[i, ] 
#       # S <- sum( (target - b1[1] - b1[2] * Xi )^2 )
#       # hi2 <- ( my - b1[1] - b1[2] * mx )^2
#       # mod <- optimise(funa, c(0, 50), n = n, ni = ni, ni2 = ni2, S = S, hi2 = hi2, tol = tol)
#       # d <- mod$minimum 
#       # tcom <- t( sx/(1+ ni * d) )
#       # A <- xx - d * tcom %*% sx
#       # B <- sxy - d * tcom %*% sy
#       # down <- A[1,1] * A[2, 2] - A[1, 2]^2
#       # b2 <- c(A[2, 2] * B[1] - A[1, 2] * B[2], - A[1, 2] * B[1] + A[1, 1] * B[2])/down
#       # ###########
#       # while ( sum( abs(b2 - b1) ) > tol ) {
#         # b1 <- b2
#         # S <- sum( (target - b1[1] - b1[2] * Xi )^2 )
#         # hi2 <- ( my - b1[1] - b1[2] * mx )^2
#         # mod <- optimise(funa, c(0, 50), n = n, ni = ni, ni2 = ni2, S = S, hi2 = hi2, tol = tol)
#         # d <- mod$minimum 
#         # tcom <- t( sx/(1+ ni * d) )
#         # A <- xx - d * tcom %*% sx
#         # B <- sxy - d * tcom %*% sy
#         # down <- A[1,1] * A[2, 2] - A[1, 2]^2
#         # b2 <- c(A[2, 2] * B[1] - A[1, 2] * B[2], - A[1, 2] * B[1] + A[1, 1] * B[2])/down
#       # }
#       # se <- (S - d * sum(ni2 * hi2/ (1 + ni * d) ) )/n
#       # seb <- A[1, 1] / down * se
#       # stat[i] <- b2[2]^2 / seb
#     # }
# 	mod <- Rfast::rint.regs(target, dataset, id, logged = TRUE) 
#     univariateModels$stat <- mod[, 1]
#     univariateModels$pvalue <- mod[, 2]  ## pf(stat, 1, n - 4, lower.tail = FALSE, log.p = TRUE)
# 
#   } else {
#     
#     sy <- as.vector( rowsum(target, id) )
#     ni <- tabulate(id)
#     my <- sy / ni
#     ni2 <- ni^2
#     Xi <- cbind(1, reps, dataset[, 1])
# 
#     funa <- function(d, n, ni, ni2, S, hi2)   sum( log1p(ni * d) ) + n * log(S - d * sum(ni2 * hi2/ (1 + ni * d) ) )    
#     stat <- numeric(D)
#     for (i in 1:D) {
#       Xi[, 3] <- dataset[, i]
#       xx <- crossprod(Xi)
#       sx <- rowsum(Xi, id)
#       sxy <- crossprod(Xi, target)
#       mx <- sx / ni
#       #############
#       mod <- .lm.fit(Xi, target)
#       b1 <- mod$coefficients
#       S <- sum( mod$residuals^2 )
#       hi2 <- ( my - mx %*% b1 )^2
#       mod <- optimise(funa, c(0, 70), n = n, ni = ni, ni2 = ni2, S = S, hi2 = hi2, tol = tol)
#       d <- mod$minimum 
#       b2 <- solve( xx - d * crossprod(sx/(1+ ni * d), sx), sxy - d * crossprod( sx, sy/(1 + ni * d) ) )   
#       k <- 2
#       while ( sum( abs(b2 - b1) ) > tol  &  k < 100) {
#         k <- k + 1
#         b1 <- b2
#         S <- sum( (target - Xi %*% b1)^2 )
#         hi2 <- ( my - mx %*% b1 )^2
#         mod <- optimise(funa, c(0, 70), n = n, ni = ni, ni2 = ni2, S = S, hi2 = hi2, tol = tol)
#         d <- mod$minimum 
#         b2 <- solve( xx - d * crossprod(sx/(1+ ni * d), sx), sxy - d * crossprod( sx, sy/(1 + ni * d) ) ) 
#       }
#       se <- (S - d * sum(ni^2 * hi2/ (1 + ni * d) ) )/n
#       seb <- solve( xx - d * crossprod(sx/(1+ ni * d), sx) )[3, 3] * se 
#       stat[i] <- b2[3]^2 / seb
#     }  
#     univariateModels$stat = stat
#     univariateModels$pvalue = pf(stat, 1, n - 4, lower.tail = FALSE, log.p = TRUE)
#   
#   }  ## else if ( is.null(reps) )
#   if ( sum(poia>0) > 0 ) {
#     univariateModels$stat[poia] = 0
#     univariateModels$pvalue[poia] = log(1)
#   }
# 
#   univariateModels
# }   


