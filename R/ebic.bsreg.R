ebic.bsreg <- function(target, dataset, test = NULL, wei = NULL, gam = NULL) {
  
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any( is.na(dataset) ) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(za){ za[which(is.na(za))] = median(za, na.rm = TRUE) ; return(za) } ) 
    } else {
      poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
      for ( i in poia )  {
        xi <- dataset[, i]
        if ( is.numeric(xi) ) {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) ) {
          xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        }
        dataset[, i] <- xi
      }
    }
  }
  
  zevar <- Rfast::check_data(dataset)
  if ( sum( zevar > 0 ) > 0 )  dataset[, zevar] <- rnorm( dim(dataset)[1] * length(zevar) )
  dataset <- as.data.frame(dataset)
  
  if ( test == "testIndReg" ) {
    result <- ebic.lm.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if ( test == "testIndPois" ) {
    result <- ebic.glm.bsreg(target, dataset, gam = gam, wei = wei, type = "poisson")
  
  } else if ( test == "testIndLogistic" ) {
    result <- ebic.glm.bsreg(target, dataset, wei = wei, gam = gam, type = "logistic")
    
  } else if ( test == "testIndBinom" ) {
    result <- ebic.glm.bsreg(target, dataset, gam = gam, wei = wei, type = "binomial")
    
  } else if ( test == "testIndNB" ) {
    result <- ebic.nb.bsreg(target, dataset, gam = gam, wei = wei)

  } else if ( test == "testIndMultinom" ) {
    result <- ebic.multinom.bsreg(target, dataset, gam = gam, wei = wei)  

  } else if ( test == "testIndOrdinal" ) {
    result <- ebic.ordinal.bsreg(target, dataset, gam = gam, wei = wei)
    
  } else if ( test == "testIndMMReg" ) {
    result <- ebic.mm.bsreg(target, dataset, gam = gam, wei = wei)
    
  } else if ( test == "censIndCR" ) {
    result <- ebic.cr.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if ( test == "censIndWR" ) {
    result <- ebic.wr.bsreg(target, dataset, gam = gam, wei = wei)
    
  } else if ( test == "censIndLLR" ) {
    result <- ebic.llr.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if ( test == "testIndBeta" ) {
    result <- ebic.beta.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if ( test == "testIndZIP" ) {
    result <- ebic.zip.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if ( test == "testIndGamma" ) {
    result <- ebic.glm.bsreg(target, dataset, gam = gam, wei = wei, type = "gamma")
  
  } else if ( test == "testIndNormLog" ) {
    result <- ebic.glm.bsreg(target, dataset, gam = gam, wei = wei, type = "gaussian")
  
  } else if ( test == "testIndTobit" ) {
    result <- ebic.tobit.bsreg(target, dataset, gam = gam, wei = wei)
    
  } else if ( test == "testIndClogit" ) {
    result <- ebic.clogit.bsreg(target, dataset, gam = gam, wei = wei)
    
  } else if ( test == "testIndSPML" )  {
    result <- ebic.spml.bsreg(target, dataset, gam = gam) 
  }
  
  back.rem <- result$info[, 1]
  back.n.tests <- dim(dataset)[2]:dim(result$mat)[1] 
  result$mat <- result$mat
  result$back.rem <- back.rem
  result$back.n.tests <- sum(back.n.tests)
  result$runtime <- result$runtime 
  result
}   
