ebicScore <- function(target, dataset, test, wei = NULL, targetID = -1, ncores = 1, gam = NULL) {
  #how many tests
  dm <- dim(dataset)
  cols <- dm[2]
  rows <- dm[1]
  bic <- numeric(cols) 
  
  if (targetID != -1 ) {
    target <- dataset[, targetID]
    dataset[, targetID] <- rbinom(rows, 1, 0.5)
  }   
  id <- NULL
  id <- Rfast::check_data(dataset)
  if ( sum(id > 0) )  dataset[, id] <- rnorm(rows * length(id) )
  
  if ( is.null(gam) ) {
    con <- 2 - log(cols) / log(rows)
  } else con <- 2 * gam
  if ( (con) < 0 )  con <- 0
  
  if ( is.null(gam) ) {
    con <- 2 - log(cols) / log(rows)
    if ( (con) < 0 )  con <- 0
  } else con <- 2 * gam
  
  for (i in 1:cols)  bic[i] <- test(target, dataset, xIndex = i, csIndex = 0, wei = wei)
  
  if ( gam != 0 ) {
    bic <- bic + con * log(cols)
  } else bic <- bic
  if ( targetID != - 1 )   bic[targetID] <- Inf
  if ( sum(id>0) > 0 )   bic[id] <- Inf
  
  list(ebic = bic)
}
