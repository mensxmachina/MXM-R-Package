\name{Cross-validation for ridge regression}
\alias{ridgereg.cv}
\title{
Cross validation for the ridge regression
}

\description{
Cross validation for the ridge regression is performed using the TT estimate of bias (Tibshirani and Tibshirani, 2009).
There is an option for the GCV criterion which is automatic. 
}

\usage{
ridgereg.cv( target, dataset, K = 10, lambda = seq(0, 2, by = 0.1), auto = FALSE, 
seed = FALSE, ncores = 1, mat = NULL )
}

\arguments{
  \item{target}{
  A numeric vector containing the values of the target variable. If the values are proportions or percentages, 
  i.e. strictly within 0 and 1 they are mapped into R using log( target/(1 - target) ).
}
  \item{dataset}{
  A numeric matrix containing the variables. Rows are samples and columns are features.
}
  \item{K}{
  The number of folds. Set to 10 by default. 
}
  \item{lambda}{
  A vector with the a grid of values of \eqn{\lambda} to be used. 
}
  \item{auto}{
  A boolean variable. If it is TRUE the GCV criterion will provide an automatic answer for the best $lambda$. Otherwise k-fold cross validation is performed.
}
  \item{seed}{
  A boolean variable. If it is TRUE the results will always be the same. 
}
 \item{ncores}{
  The number of cores to use. If it is more than 1 parallel computing is performed.  
}
 \item{mat}{
  If the user has its own matrix with the folds, he can put it here. It must be a matrix with K columns, each column is a fold and it contains the positions of the data,
 i.e. numbers, not the data. For example the first column is c(1,10,4,25,30), the second is c(21, 23,2, 19, 9) and so on.
}
}

\details{
 The lm.ridge command in MASS library is a wrapper for this function. If you want a fast choice of \eqn{\lambda}, then specify auto = TRUE and the \eqn{\lambda} which minimizes the generalised cross-validation criterion will be returned. Otherise a k-fold cross validation is performed and the estimated performance is bias corrected as suggested by Tibshirani and Tibshirani (2009). 
}

\value{
A list including:
\item{mspe}{
If auto is FALSE the values of the mean prediction error for each value of \eqn{\lambda}. 
}
\item{lambda}{
If auto is FALSE the \eqn{\lambda} which minimizes the MSPE. 
}
\item{performance}{
If auto is FALSE the minimum bias corrected MSPE along with the estimate of bias. 
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

\references{
Hoerl A.E. and R.W. Kennard (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1):55-67.

Brown P. J. (1994). Measurement, Regression and Calibration. Oxford Science Publications.

Tibshirani R.J., and Tibshirani R. (2009). A bias correction for the minimum error rate in cross-validation. The Annals of Applied Statistics 3(2): 822-829.

}
\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\note{
The values can be extracted with the $ symbol, i.e. this is an S3 class output. 
}

\seealso{
\code{ \link{ridge.reg} }
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(runif(200 * 40, 1, 100), nrow = 200 ) 
#the target feature is the last column of the dataset as a vector
target <- dataset[, 40]
a1 <- ridgereg.cv(target, dataset, auto = TRUE)
a2 <- ridgereg.cv( target, dataset, K = 10, lambda = seq(0, 1, by = 0.1) )
}

\keyword{ Ridge regression }