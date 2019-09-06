shd <- function(est, true) {
  mat <- table(est, true)
  a1 <- sum(est == 0  &  true == 1) + sum(est == 1 & true == 0)
  a2 <- sum(est == 0 & true == 2) + sum(est == 2 & true == 0) 
  a3 <- sum(est == 1 & true == 2) + sum(est == 2 & true == 1)
  a4 <- sum(est == 3 & true == 2) 
  list(mat = mat, dis = a1 + a2 + a3 + a4)
}

   ## SHD as defined in Tsamardinos et al. (2006)
   ##  True   Estimated   Penalty
   ##  -                  1   (a1)
   ##         -           1   (a1)
   ##  ->                 1   (a2)
   ##         <-          1   (a2)
   ##  ->     -           1   (a3)
   ##  -      <-          1   (a3)
   ##  ->     <-          1   (a4)               
 



   

  