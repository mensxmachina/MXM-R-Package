\name{Regression modeler}
\alias{modeler}
\title{
Generic regression modelling function
}

\description{
Generic regression modelling function.
}

\usage{
modeler(target, dataset = NULL, test = "testIndFisher")
}

\arguments{
\item{target}{
The target (dependent) variable. It can be a numerical variable, factor, ordinal factor, percentages, or time to event. 
}
\item{dataset}{
The predictor variable(s). It can be a vector, a matrix with continuous only variables. If there are no predictor variables leave this NULL.
}
\item{test}{
Unlike \code{\link{reg.fit}} this accepts the test. The test argument is exactly like in all feature selection methods. 
This function accepts the following:  "testIndReg", "testIndPois", "testIndNB", "testIndLogistic", 
"testIndMMReg", "testIndRQ", "testIndBinom", "censIndCR", "censIndWR", "censIndLLR", "testIndBeta", "testIndGamma, "testIndNormLog", 
"testIndTobit", "testIndQPois", "testIndQBinom", "testIndFisher", "testIndMultinom" and "testIndOrdinal". 
}
}

\details{
This is a generic regression function designed for continuous predictor variables only. It was useful for me so I decided to epxort it.
}

\value{
A list including:
\item{mod}{
The fitted model.
}
\item{dev}{
The deviance. For some models though ("testIndMMReg", "testIndRQ", "censIndCR", "censIndWR", "testIndTobit", "testIndBeta", "testIndNB",
""testIndQPois", "testIndQBinom") this contains twice the log-likelihood.
}
\item{bic}{
The BIC of the model. This is NA for the "testIndQPois" and "testIndQBinom" because they are quasi likhelidood models and hence have no BIC. 
}
\item{res}{
The residuals of the fitted model.
}
}

\references{
Almost the same as in \code{\link{CondIndTests}}.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{reg.fit}, \link{fbedreg.bic}, \link{mmpc.model}, \link{ridge.reg} }
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(runif(100 * 5, 1, 100), nrow = 100 ) 
#the target feature is the last column of the dataset as a vector
target <- dataset[, 1]
dataset <- dataset[, -1]
a <- modeler(target, dataset)
}

\keyword{ Regression model }