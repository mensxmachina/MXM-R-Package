\name{Ridge regression coefficients plot}
\alias{ridge.plot}
\title{
Ridge regression
}
\description{
A plot of the regularised parameters is shown. 

}
\usage{
ridge.plot(target, dataset, lambda = seq(0, 5, by = 0.1) ) 
}

\arguments{
\item{target}{
A numeric vector containing the values of the target variable. If the values are proportions or percentages, 
i.e. strictly within 0 and 1 they are mapped into R using log( target/(1 - target) ). In any case, they must be continuous only.
}
\item{dataset}{
A numeric matrix containing the continuous variables. Rows are samples and columns are features.
}
\item{lambda}{
A grid of values of the regularisation parameter \eqn{\lambda}.
}
}

\details{
For every value of \eqn{\lambda} the coefficients are obtained. They are plotted versus the \eqn{\lambda} values. 
}

\value{
A plot with the values of the coefficients as a function of \eqn{\lambda}.
}

\references{
Hoerl A.E. and R.W. Kennard (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1): 55-67.

Brown P. J. (1994). Measurement, Regression and Calibration. Oxford Science Publications.
}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr>, Vincenzo Lagani <vlagani@csd.uoc.gr> and Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ridge.reg}, \link{ridgereg.cv}
}
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix( runif(300 * 20, 1, 20), nrow = 300 ) 
#the target feature is the last column of the dataset as a vector
target <- dataset[, 20]
dataset <- dataset[, -20]
ridge.plot(target, dataset)
}

\keyword{ Ridge regression }