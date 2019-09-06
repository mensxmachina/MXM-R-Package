\name{Ridge regression}
\alias{ridge.reg}
\title{
Ridge regression
}

\description{
Regularisation via ridge regression is performed. 

}
\usage{
ridge.reg(target, dataset, lambda, B = 1, newdata = NULL)
}

\arguments{
  \item{target}{
  A numeric vector containing the values of the target variable. If the values are proportions or percentages, 
  i.e. strictly within 0 and 1 they are mapped into R using log( target/(1 - target) ).
}
  \item{dataset}{
  A numeric matrix containing the variables. Rows are samples and columns are features.
}
  \item{lambda}{
  The value of the regularisation parameter \eqn{\lambda}.
}
  \item{B}{
  Number of bootstraps. If B = 1 no bootstrap is performed and no standard error for the regression coefficients is returned. 
}
  \item{newdata}{
  If you have new data and want to predict the value of the target put them here, otherwise, leave it NULL.
}
}

\details{
There is also the lm.ridge command in MASS library if you are interested in ridge regression. 
}

\value{
A list including:
\item{beta}{
The regression coefficients if no bootstrap is performed. If bootstrap is performed their standard error appears as well.
}
\item{seb}{
The standard erorr of the regression coefficients. If bootstrap is performed their bootstrap estimated standard error appears.
}
\item{est}{
The fitted values if no new data are available. If you have used new data these will be the predicted target values. 
}
}

\references{
Hoerl A.E. and R.W. Kennard (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1): 55-67.

Brown P. J. (1994). Measurement, Regression and Calibration. Oxford Science Publications.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ridgereg.cv} }
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(runif(100 * 30, 1, 100), nrow = 100 ) 
#the target feature is the last column of the dataset as a vector
target <- dataset[, 10]
dataset <- dataset[, -10]
a1 <- ridge.reg(target, dataset, lambda = 0.5, B = 1, newdata = NULL)
a2 <- ridge.reg(target, dataset, lambda = 0.5, B = 100, newdata = NULL) 
}

\keyword{ Ridge regression }