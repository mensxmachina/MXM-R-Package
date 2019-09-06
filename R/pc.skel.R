pc.skel <- function(dataset, method = "pearson", alpha = 0.01, rob = FALSE, R = 1, stat = NULL, ini.pvalue = NULL) {
  ## dataset contains the data, it must be a matrix 
  ## type can be either "pearson" or "spearman" for continuous variables OR
  ## "cat" for categorical variables
  ## alpha is the level of significance, set to 0.05 by default
  ## rob is TRUE or FALSE and is supported by type = "pearson" only, i.e. it is to be used
  ## for robust estimation of Pearson correlation coefficient only
  title <- deparse( substitute(dataset) )
  nam <- colnames(dataset)
  n <- dim(dataset)[2]
  if ( method != "distcor"  &  method != "comb.fast" &  method != "comb.mm"  &  rob == FALSE )  {
    res <- Rfast::pc.skel(dataset = dataset, method = method, alpha = alpha, R = R)
    if ( is.null( nam ) ) {
      colnames(res$G) = rownames(res$G) = paste("X", 1:n, sep = "")
      colnames(res$stat) = rownames(res$stat) = paste("X", 1:n, sep = "")
      colnames(res$pvalue) = rownames(res$pvalue) = paste("X", 1:n, sep = "")
    } else {
      colnames(res$G) = rownames(res$G) = nam
      colnames(res$stat) = rownames(res$stat) = nam
      colnames(res$pvalue) = rownames(res$pvalue) = nam
    }  
    info <- summary( Rfast::rowsums(res$G) )
    density <- sum(res$G) / n / ( n - 1 ) 
    res$density = density
    res$info = info
    res$title = title 
    
  } else {
    
    alpha <- log(alpha)
    if ( any( is.na(dataset) ) ) {
      #dataset = as.matrix(dataset);
      warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
      if ( is.matrix(dataset) )  {
        dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
      } else {
        poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
        for ( i in poia )  {
          xi <- dataset[, i]
          if ( is.numeric(xi) ) {                    
            xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
          } else if ( is.factor( xi ) )   xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
          dataset[, i] <- xi
        }
      }
    }  ## end if ( any(is.na(dataset)) ) 
    if ( method == "pearson" ) {
      ci.test <- condi 
      type <- method
      rob <- TRUE
    } else if ( method == "distcor" ) {
      ci.test <- dist.condi
      type <- NULL
      rob <- FALSE
    } else if (method == "comb.fast" ) {
      ci.test <- ci.fast
      type <- NULL
      rob <- FALSE
    } else if (method == "comb.mm" ) {
      ci.test <- ci.mm
      type <- NULL
      rob <- FALSE
    }  
    dm <- dim(dataset)
    n <- dm[2]
    m <- dm[1]
    k <- 0  ## initial size of the conditioning set
    G <- matrix(2, n, n)  # 3 sep-set indicates a subset of variables which eliminate given edge  
    ## If an element has the number 2 it means there is connection, otherwiser it will have 0
    diag(G) = -100
    durat = proc.time()
    
    if ( method == "pearson") {
      
      if ( is.null(stat) | is.null(ini.pvalue) ) {
        stat = pv = matrix(0, n, n)
        for ( i in 1:c(n - 1) ) {
          for ( j in c(i + 1):n ) {
            ro <- condi(i, j, 0, dataset, type = "pearson", rob = TRUE) 
            stat[i, j] = ro[1]
            stat[j, i] = ro[1]
            pv[i, j] = ro[2]
            pv[j, i] = ro[2]
          }
        }  ## end for ( i in 1:c(n - 1) )
      } else    pv <- ini.pvalue
      pvalue = pv  ## p-values
      dof = matrix(m - 3, n, n)
      stadf = stat / dof
      
    } else if ( method == "distcor" )  {
      if ( is.null(stat) | is.null(ini.pvalue) ) {
        stat = pv = matrix(0, n, n)
        for ( i in 1:c(n - 1) ) {
          for ( j in c(i + 1):n ) {
            ro <- dist.condi(i, j, 0, dataset, R = R) 
            stat[i, j] <- ro[1]
            stat[j, i] <- ro[1]
            pv[i, j] <- ro[2]
            pv[j, i] <- ro[2]
          }
        }  ## end for ( i in 1:c(n - 1) )
        pvalue <- pv  ## p-values
        stadf <- stat 
      } else  pv <- ini.pvalue
    } else if ( method == "comb.fast" )  {
      
      if ( is.null(stat) | is.null(ini.pvalue) ) {
        stat = pv = matrix(0, n, n)
        for ( i in 1:c(n - 1) ) {
          for ( j in c(i + 1):n ) {
            ro <- ci.fast(i, j, NULL, dataset) 
            stat[i, j] <- ro[1]
            stat[j, i] <- ro[1]
            pv[i, j] <- ro[2]
            pv[j, i] <- ro[2]
          }
        }  ## end for ( i in 1:c(n - 1) )
        pvalue <- pv  ## p-values
        stadf <- stat 
      } else  pv <- ini.pvalue
    } else if ( method == "comb.mm" )  {
      if ( is.null(stat) | is.null(ini.pvalue) ) {
        stat = pv = matrix(0, n, n)
        for ( i in 1:c(n - 1) ) {
          for ( j in c(i + 1):n ) {
            ro <- ci.mm(i, j, NULL, dataset) 
            stat[i, j] <- ro[1]
            stat[j, i] <- ro[1]
            pv[i, j] <- ro[2]
            pv[j, i] <- ro[2]
          }
        }  ## end for ( i in 1:c(n - 1) )
        pvalue <- pv  ## p-values
        stadf <- stat 
      } else  pv <- ini.pvalue
    }
    ini.pvalue <- pv

    pv <- pvalue
    pv[ lower.tri(pv) ] = 2 
    G[pvalue > alpha] <- 0   ## removes edges from non significantly related pairs
    
    diag(pvalue) = diag(pv) = 0
    ina <- 1:n 
    sep <- list()
    n.tests <- NULL
    #### some more initial stuff 
    dial <- which( pv <= alpha, arr.ind = TRUE )
    zeu <- cbind( dial, stadf[ dial ], pv[ dial ] )  ## all significant pairs of variables
    zeu <- zeu[ order( - zeu[, 4], zeu[, 3] ), , drop = FALSE] ## order of the pairs based on their strength
    duo = dim(zeu)[1]  ## number of pairs to be checked for conditional independence
    n.tests[1] = n * (n - 1) / 2
    #### main search
    if (duo == 0) {
      diag(G) <- 0
      res <- list(kappa = k, G = G) 
    } else {
      
      ell <- 2
      ## Execute PC algorithm: main loop
      while ( k < ell & duo > 0 )  {
        k <- k + 1   ## size of the seperating set will change now
        tes <- 0
        met <- matrix(0, duo, k + 2)
        
        for ( i in 1:nrow(zeu) ) {
          candpair <- zeu[i, 1:2]
          adjx <- ina[ G[ candpair[1], ] == 2 ]   ;   lx <- length(adjx)  ## adjacents to x
          adjy <- ina[ G[ candpair[2], ] == 2 ]   ;   ly <- length(adjy)  ## adjacents to y
          
          if ( lx >= k )  {
            pvalx <- pvalue[ candpair[1], adjx ]
            infox <- cbind( adjx, pvalx)
            infox <- infox[ order( - pvalx ), ]
            if ( !is.matrix(infox) ) {
              samx <- cbind( infox[1], infox[2] )
            } else  samx <- cbind( t( combn(infox[, 1], k) ), t( combn(infox[, 2], k) ) )  ## factorial, all possible unordered pairs
          }  ## end if ( lx >= k )
          
          if ( ly >= k ) {
            pvaly <- pvalue[ candpair[2], adjy ]
            infoy <- cbind(adjy, pvaly)
            infoy <- infoy[ order( - pvaly ), ]
            if ( !is.matrix(infoy) ) {
              samy <- cbind( infoy[1], infoy[2] )
            } else  samy = cbind( t( combn(infoy[, 1], k) ), t( combn(infoy[, 2], k) ) )  ## factorial, all possible unordered pairs
          }  ## end if ( ly >= k )
          
          if ( !is.null(samx) ) sx = 1  else sx = 0
          if ( !is.null(samy) ) sy = 1  else sy = 0 
          sam <- rbind( samx * sx, samy * sy ) 
          sam <- unique(sam)
          ## sam contains either the sets of k neighbours of X, or of Y or of both
          ## if the X and Y have common k neighbours, they are removed below
          rem = intersect( zeu[i, 1:2], sam )
          if ( length(rem) > 0 ) {
            pam <- list()
            for ( j in 1:length(rem) ) {
              pam[[ j ]] <- as.vector( which(sam == rem[j], arr.ind = TRUE)[, 1] ) 
            }
          }  ## end if ( length(rem) > 0 )
          
          pam <- unlist(pam)
          sam <- sam[ - pam, ] 
          
          if ( !is.matrix(sam) ) {
            sam = matrix( sam, nrow = 1 ) 
          } else if ( nrow(sam) == 0 ) {
            G = G 
            
          } else { 
            if ( k == 1 ) {
              sam = sam[ order( sam[, 2 ] ), ]
            } else {
              an <- t( apply(sam[, -c(1:2)], 1, sort, decreasing = TRUE) )
              sam <- cbind(sam[, 1:2], an)
              nc <- ncol(sam)
              sam2 <- as.data.frame( sam[, nc:1] )     
              sam2 <- sam2[ do.call( order, as.list( sam2 ) ), ] 
              sam <- as.matrix( sam2[, nc:1] )
            }  ##  end if (k == 1) 
          }  ## end if if ( !is.matrix(sam) )
          
          if ( dim(sam)[1] == 0 ) {
            G <- G  
          } else {
            a <- ci.test( candpair[1], candpair[2], sam[1, 1:k], dataset, type = type, rob = rob, R = R )
            b <- a[2]
            if ( a[2] > alpha ) {
              G[ candpair[1], candpair[2] ] = 0  ## remove the edge between two variables
              G[ candpair[2], candpair[1] ] = 0  ## remove the edge between two variables 
              met[i, ] <- c( sam[1, 1:k], a[1:2] )
              tes <- tes + 1 
            } else {
              m <- 1
              while ( a[2] < alpha  &  m < nrow(sam) ) {
                m <- m + 1
                a <- ci.test( candpair[1], candpair[2], sam[m, 1:k], dataset, type = type, rob = rob, R = R )
                b <- c(b, a[2])
                tes <- tes + 1
              }  ## end while ( a[2] < alpha  &  m < nrow(sam) )
              if (a[2] > alpha) {
                G[ candpair[1], candpair[2] ] = 0  ## remove the edge between two variables
                G[ candpair[2], candpair[1] ] = 0  ## remove the edge between two variables
                met[i, ] = c( sam[m, 1:k], a[1:2] ) 
              }  ## end if (a[2] > alpha)
            }  ## end if ( a[2] > alpha )
            pvalue[ candpair[1], candpair[2] ] = max(b, pvalue[ candpair[1], candpair[2] ] )
            pvalue[ candpair[2], candpair[1] ] = max(b, pvalue[ candpair[1], candpair[2] ] )
          }  ## end if ( dim(sam)[1] == 0 ) 
          sam = samx = samy = NULL
        }  ## end for ( i in 1:nrow(zeu) )
        
        ax = ay = list()
        lx = ly = numeric( duo )
        for ( i in 1:duo ) {
          ax[[ i ]] = ina[ G[ zeu[i, 1], ] == 2 ]  ;  lx[i] = length( ax[[ i ]] )
          ay[[ i ]] = ina[ G[ zeu[i, 2], ] == 2 ]  ;  ly[i] = length( ay[[ i ]] ) 
        }
        ell <- max(lx, ly)
        id <- which( rowSums(met) > 0 )
        if (length(id) == 1) {
          sep[[ k ]] = c( zeu[id, 1:2], met[id, ] )
        } else  sep[[ k ]] = cbind( zeu[id, 1:2], met[id, ] )
        zeu <- zeu[-id, , drop = FALSE]  
        duo <- dim(zeu)[1]
        n.tests[ k + 1 ] = tes
      }  ## end while ( k <= ell & duo > 0 )
      
      G <- G/2
      diag(G) <- 0
      durat <- proc.time() - durat
      ###### end of the algorithm
      for ( l in 1:k ) { 
        if ( is.matrix(sep[[ l ]]) ) {
          colnames( sep[[ l ]] ) <- c("X", "Y", paste("SepVar", 1:l), "stat", "logged.p-value")
          #sepa =  sepa[ order(sepa[, 1], sepa[, 2] ), ]
        } else {
          if ( length(sep[[ l ]]) > 0)   names( sep[[ l ]] ) <- c("X", "Y", paste("SepVar", 1:l), "stat", "logged.p-value")
        }
      }   ## end for ( l in 1:k ) 
      #######################
    }  
    n.tests <- n.tests[ n.tests>0 ]
    k <- length(n.tests) - 1
    sepset <- list()
    
    if (k == 0) {
      sepset = NULL
    } else {
      for ( l in 1:k ) {
        if ( is.matrix( sep[[ l ]] ) )  {
          nu <- nrow( sep[[ l ]] )
          if ( nu > 0 ) sepset[[ l ]] <- sep[[ l ]][1:nu, ]
        } else sepset[[ l ]] = sep[[ l ]]    
      }
    }  ## end if (k == 0) 
    names(n.tests) <- paste("k=", 0:k, sep ="")
    info <- summary( Rfast::rowsums(G) )
    density <- sum(G) / n / ( n - 1 ) 
    if ( is.null( colnames(dataset) ) ) {
      colnames(G) = rownames(G) = paste("X", 1:n, sep = "")
      colnames(stat) = rownames(stat) = paste("X", 1:n, sep = "")
      colnames(pvalue) = rownames(pvalue) = paste("X", 1:n, sep = "")
    } else {
      colnames(G) = rownames(G) = nam
      colnames(stat) = rownames(stat) = nam
      colnames(pvalue) = rownames(pvalue) = nam
    }  
    res <- list(stat = stat, ini.pvalue = ini.pvalue, pvalue = pvalue, runtime = durat, kappa = k, n.tests = n.tests, density = density, info = info, G = G, sepset = sepset, title = title )
  }  ## end if ( method != distcor & rob = FALSE )
  ##################
  res
}