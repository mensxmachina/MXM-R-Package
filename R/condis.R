condis <- function(ind1, ind2, cs1, cs2, Var, dat, type = "pearson", rob = FALSE, max_k = 2, R = 1 ) {
  
  if ( type == "pearson" | type == "spearman" ) {
    if (R == 1) {
      funa <- pearson_condis
    } else funa <- pearson_condis.rob    
  } else if (type == "cat") {
    funa <- cat_condis
    type <- Rfast::colrange(dat, cont = FALSE)
  } else if (type == "distcor") {
    funa <- distcor_condis
  } else if (type == "ci.mm" | type == "ci.fast") {
    funa <- comb_condis
  } 
  
  funa(ind1 = ind1, ind2 = ind2, cs1 = cs1, cs2 = cs2, Var = Var, dat = dat, type = type, rob = rob, max_k = max_k, R = R)
}

