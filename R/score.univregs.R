score.univregs <- function(target, dataset, test) {
  
  univariateModels <- list();
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]

  ## Beta regression 
  if ( identical(test, testIndBeta) ) {
    mod <- Rfast::score.betaregs(target, dataset, logged = TRUE )
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]
	
  ## Negative Binomial 
  } else if ( identical(test, testIndNB) ) {
    mod <- Rfast::score.negbinregs(target, dataset, logged = TRUE )
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]
	
  ## Poisson
  } else if ( identical(test, testIndPois) ) {
    mod <- Rfast::score.glms(target, dataset, oiko = "poisson", logged = TRUE )
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]
	
  ## logistic regression
  } else if ( identical(test, testIndLogistic) ) { 
    mod <- Rfast::score.glms(target, dataset, oiko = "binomial", logged = TRUE)  
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]
	
    ## multinomial regression
  } else if ( identical(test, testIndMultinom) ) { 
    mod <- Rfast::score.multinomregs(target, dataset, logged = TRUE)  
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]
	
	## Gamma regression
  } else if ( identical(test, testIndGamma) ) { 
    mod <- Rfast::score.gammaregs(target, dataset, logged = TRUE)  
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]
	
  # Weibull
  } else if ( identical(test, censIndWR) ) {
    mod <- Rfast::score.weibregs(target, dataset, logged = TRUE )
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]
  } else univariateModels <- NULL
  
  if ( !is.null(univariateModels) )  {
    id <- which( is.na(univariateModels$stat) )
    univariateModels$stat[id] <- 0
    univariateModels$pvalue[id] <- 1
  }
  univariateModels
  
}