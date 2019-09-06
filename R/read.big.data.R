read.big.data <- function(path,  sep = ",", header = FALSE) {
  bigmemory::read.big.matrix(path, type = "double", header = header)
}
 