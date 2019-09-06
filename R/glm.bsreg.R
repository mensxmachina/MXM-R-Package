glm.bsreg <- function(target, dataset, threshold = 0.05, wei = NULL, test = NULL) {
  
  threshold <- log(threshold)
  dm <- dim(dataset)
  if ( is.null(dm) ) {
    n <- length(target)
    p <- 1
  } else {
    n <- dm[1]  ## sample size 
    p <- dm[2]  ## number of variables
  }  
  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. 
    No backward procedure was attempted")
  } else {
  #check for NA values in the dataset and replace them with the variable median or the mode
    if ( any(is.na(dataset)) ) {
      #dataset = as.matrix(dataset);
      warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
      if ( is.matrix(dataset) )  {
        dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
      } else {
        poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
        for( i in poia )  {
          xi <- dataset[, i]
          if( is.numeric(xi) )   {                    
            xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
          } else if ( is.factor( xi ) ) {
            xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
          }
          dataset[, i] <- xi
        }
      }
    }
    dataset <- as.data.frame(dataset)
    
    if ( is.null(test) ) {
      if ( is.matrix(target) )  {
        test =="testIndBinom" 
      } else if ( is.factor(target)  |  length( unique(target) ) == 2 ) {
        ci_test <- test <- "testIndLogistic"
      } else if ( length( unique(target) ) > 2  ) {
        if ( sum( round(target) - target ) == 0 ) {
          ci_test <- test <- "testIndPois"
        } else  ci_test <- test <- "testIndReg"  
      }
    }
    ##################################
    # target checking and initialize #
    ################################## 
    if ( test == "testIndBinom" )  {
      
      ci_test <- "testIndBinom"
      wei <- target[, 2]
      ywei <- target[, 1] / wei
      
      tic <- proc.time()
      ini <- glm( ywei ~., data = dataset, family = binomial(logit), weights = wei, y = FALSE, model = FALSE )
      tab <- drop1( ini, test = "Chisq" )
      dof <- tab[-1, 1]
      stat <- tab[-1, 4]

      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 
      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 
      
      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = ci_test, final = ini ) 
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE]

      i <- 1  
      
          while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
            
            ini <- glm( ywei ~., data = dat, family = binomial(logit), weights = wei, y = FALSE, model = FALSE )
            i <- i + 1
            k <- p - i + 1
            if ( k == 1 ) {
              mod <- glm(ywei ~ 1, family = binomial(logit), weights = wei)
              stat <- 2 * ( logLik(ini) - logLik(mod) )
              dof <- length( coef(ini) ) - length( coef(mod) ) 
              pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
              
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], pval, stat) )
                dat <- dataset[, -info[, 1], drop = FALSE ] 
                mat <- mat[-sel, , drop = FALSE] 
              } else {
                info <- rbind(info, c(0, -10, -10)) 
                final <- ini
                mat[, 2:3] <- c(pval, stat)
              }  
            } else {
              tab <- drop1( ini, test = "Chisq" )
              dof <- tab[-1, 1]
              stat <- tab[-1, 4]
              mat[, 2:3] <- cbind( pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
              sel <- which.max( mat[, 2] )
            
             if ( mat[sel, 2] < threshold ) {
               final <- ini
               info <- rbind(info, c(0, -10, -10) )
              
             } else {
               info <- rbind(info, mat[sel, ] )
               mat <- mat[-sel, , drop = FALSE] 
               dat <- dataset[, -info[, 1], drop = FALSE ] 
            }
           }  
          }
          
        runtime <- proc.time() - tic		
        info <- info[ info[, 1] > 0, , drop = FALSE]
        res <- list(runtime = runtime, info = info, mat = mat, ci_test = ci_test, final = final ) 

      }   
      
    } else {

   ############ 
   ###  Poisson or logistic regression
   ###########
	
    if ( test == "testIndLogistic" | test == "testIndPois" ) {

	    tic <- proc.time()

      if ( length( unique(target) ) == 2 ) {
        oiko <- binomial(logit)
        ci_test <- "testIndLogistic"
        
      } else if ( length( unique(target) )  > 2  &  sum( round(target) - target ) == 0 ) {
        oiko <- poisson(log)
        ci_test <- "testIndPois"
      }
      ini <- glm( target ~.,  data = dataset, family = oiko, weights = wei, y = FALSE, model = FALSE )
      tab <- drop1( ini, test = "Chisq" )
      dof <- tab[-1, 1]
      stat <- tab[-1, 4]
      
      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 
      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 

      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = ci_test, final = ini ) 
        
      } else {
       
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE] 

      i <- 1  

         while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
           i <- i + 1
           k <- p - i + 1
           ini <- glm( target ~., data = dat, family = oiko, weights = wei, y = FALSE, model = FALSE )
           
           if ( k == 1 ) {
             mod <- glm(target ~ 1, data = dat, family = oiko, weights = wei)
             stat <- 2 * ( logLik(ini) - logLik(mod) )
             dof <- length( coef(ini) ) - length( coef(mod) ) 
             pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
             
             if (pval > threshold ) {
               final <- "No variables were selected"
               info <- rbind(info, c(mat[, 1], pval, stat) )
               dat <- dataset[, -info[, 1], drop = FALSE ]
               mat <- mat[-sel, , drop = FALSE] 
             } else {
               info <- rbind(info, c(0, -10, -10)) 
               final <- ini
               mat[, 2:3] <- c(pval, stat)
             }
             
           } else {
           
             tab <- drop1( ini, test = "Chisq" )
             dof <- tab[-1, 1]
             stat <- tab[-1, 4]
             mat[, 2:3] <- cbind( pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
             sel <- which.max( mat[, 2] )
      
             if ( mat[sel, 2] < threshold ) {
               final <- ini
               info <- rbind(info, c(0, -10, -10) )

             } else {
               info <- rbind(info, mat[sel, ] )
               mat <- mat[-sel, , drop = FALSE] 
               dat <- dataset[, -info[, 1], drop = FALSE ]
             }
           }  
         }

      runtime <- proc.time() - tic		
      info <- info[ info[, 1] > 0, , drop = FALSE]
      res <- list(runtime = runtime, info = info, mat = mat, ci_test = ci_test, final = final ) 
      
    }
    ############ 
    ###  Linear regression
    ############
    } else { 
      
	    tic <- proc.time()
      if (test == "testIndReg")  {
        ci_test <- "testIndReg"
	    
        ini <- lm( target ~., data = dataset, weights = wei, y = FALSE, model = FALSE )
        df2 <- n - length( coef(ini) )
        tab <- drop1( ini, test = "F" )
        df1 <- tab[-1, 1]
        stat <- tab[-1, 5]
        mat <- cbind(1:p, pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE), stat )

      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 

      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 

      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = ci_test, final = ini ) 
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE] 

      i <- 1  

          while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
            i <- i + 1
            k <- p - i + 1
            
            ini <- lm( target ~.,  data = dat, weights = wei, y = FALSE, model = FALSE )
            df2 <- n - length( coef(ini) )
            tab <- drop1( ini, test = "F" )
            df1 <- tab[-1, 1]
            stat <- tab[-1, 5]
            pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE)
            
            if ( k == 1 ) {
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], pval, stat) )
                dat <- dataset[, -info[, 1], drop = FALSE ]
                mat <- mat[-sel, , drop = FALSE] 
              } else {
                info <- rbind(info, c(0, -10, -10)) 
                final <- ini
                mat[, 2:3] <- c(pval, stat)
              }  
            } else {
              
              mat[, 2:3] <- cbind( pval, stat )
              
              sel <- which.max( mat[, 2] )
              if ( mat[sel, 2] < threshold ) {
                final <- ini
                info <- rbind(info, c(0, -10, -10) )
              } else {
                info <- rbind(info, mat[sel, ] )
                mat <- mat[-sel, , drop = FALSE] 
                dat <- dataset[, -info[, 1], drop = FALSE] 
              }  ## end if ( mat[sel, 2] < threshold )
            }  ## end if ( k == 1 )
            
          }  ## end while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )
        
        runtime <- proc.time() - tic		
        info <- info[ info[, 1] > 0, , drop = FALSE]
        res <- list(runtime = runtime, info = info, mat = mat, ci_test = ci_test, final = final ) 

      }  ## end if ( mat[sel, 2] < threshold )
      
     } else if (test == "testIndMMReg")  {
      
      tic <- proc.time()
      ci_test <- "testIndMMReg"

        mat <- matrix(ncol = 3, nrow = p)
        mat[, 1] <- 1:p
        ini <- MASS::rlm( target ~., data = dataset, method = "MM", maxit = 2000)
        lik1 <- as.numeric( logLik(ini) )
        dofini <- length( coef(ini) ) 
        if (p == 1) {
          fit2 <- MASS::rlm( target ~ 1, data = dataset, method = "MM", maxit = 2000 )
          mat[1, 3] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
          dof <- dofini - length( coef(fit2) )
          mat[1, 2] = pchisq( mat[1, 3], dof, lower.tail = FALSE, log.p = TRUE )  
        } else {
          for (j in 1:p) {
            fit2 <- MASS::rlm( target ~., data = dataset[, -j, drop = FALSE], method = "MM", maxit = 2000 )
            mat[j, 3] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
            dof <- dofini - length( coef(fit2) )
            mat[j, 2] = pchisq( mat[j, 3], dof, lower.tail = FALSE, log.p = TRUE )  
          }
        }  ## end if (p == 1) 			
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 
      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 
      
      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = ci_test, final = ini ) 
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE] 
        
        i <- 1  
        
          while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
            i <- i + 1
            k <- p - i + 1

              ini <- MASS::rlm( target ~., data = dat, method = "MM", maxit = 2000 )
              lik1 <- as.numeric( logLik(ini) )
              dofini <- length( coef(ini) )		
              
              if ( k == 1 ) {
                mod <- MASS::rlm(target ~ 1, data = dat, method = "MM", maxit = 2000)
                stat <- 2 * abs( lik1 - logLik(mod) )
                dof <- dofini - length( coef(mod) ) 
                pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
                
                if (pval > threshold ) {
                  final <- "No variables were selected"
                  info <- rbind(info, c(mat[, 1], pval, stat) )
                  dat <- dataset[, -info[, 1], drop = FALSE ]
                  mat <- mat[-sel, , drop = FALSE] 
                } else {
                  info <- rbind(info, c(0, -10, -10)) 
                  final <- ini
                  mat[, 2:3] <- c(pval, stat)
                }  
              } else { 
                
                for (j in 1:k) {
                  fit2 <- MASS::rlm( target ~., data = dat[, -j, drop = FALSE], method = "MM", maxit = 2000 )
                  mat[j, 3] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
                  dof <- dofini - length( coef(fit2) )
                  mat[j, 2] <- pchisq( mat[j, 3], dof, lower.tail = FALSE, log.p = TRUE )  
                }  
                
                sel <- which.max( mat[, 2] )
                if ( mat[sel, 2] < threshold ) {
                  final <- ini
                  info <- rbind(info, c(0, -10, -10) )
                } else {
                  info <- rbind(info, mat[sel, ] )
                  mat <- mat[-sel, , drop = FALSE] 
                  dat <- dataset[, -info[, 1] ,drop = FALSE]
                }  ## end if ( mat[sel, 2] < threshold )
              }
          }  ## end while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 ) 
          
          runtime <- proc.time() - tic		
          info <- info[ info[, 1] > 0, , drop = FALSE]
          res <- list(runtime = runtime, info = info, mat = mat, ci_test = ci_test, final = final ) 

      }  ## end if ( mat[sel, 2] < threshold )
      
      }
      
     }   ## end else 
    
      
    }

    
  }  
  res
}    

     


 
  
    