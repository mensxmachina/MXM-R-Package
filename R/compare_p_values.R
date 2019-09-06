compare_p_values = function(pval, pval2, stat, stat2) {
  if ( length(pval) == 0 | length(pval2) == 0 | length(stat) == 0 | length(stat2) == 0 ) {
    return(FALSE);
  } else {
    if ( is.na(pval2)  | is.na(stat2)  | is.na(pval) | is.na(stat) ) {
      pval2 = 0.0;
      return(FALSE);       #(pval < pval2);
    } else {
      if ( pval == pval2 ) {
        return (stat > stat2);
      } else {
       return ( pval < pval2);
      }
    }
  }
}