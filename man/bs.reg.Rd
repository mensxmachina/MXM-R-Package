\name{Backward selection regression}
\alias{bs.reg}
\title{
Variable selection in regression models with backward selection
}

\description{
Variable selection in regression models with backward selection
}

\usage{
bs.reg(target, dataset, threshold = 0.05, wei = NULL, test = NULL, user_test = NULL)
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details. For the gamma regression this must a vector with strictly positive numbers (no zeros allowed). 
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). In either case, only two cases are avaialble, either all data are continuous, or categorical. 
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
The regression model to use. Available options are most of the tests for SES and MMPC. 
The ones NOT available are "censIndER", "testIndMVreg", "testIndSpearman". 
If you want to use multinomial or ordinal logistic regression,  make sure your target is factor. 
See also \code{\link{SES}} and \code{\link{CondIndTests}} for the tests.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, the "test" argument is ignored.
}
  %\item{maxit}{
    %The maximum number of iterations required by the MM estimators in the robust binary or poisson regrression. This argument is %taken into account only if the argument robust is set to TRUE.
    %}
}
\details{
If the sample size is less than the number of variables a meesage will appear and no backward regression is performed. 
}

\value{
The output of the algorithm is S3 object including:
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{info}{
A matrix with the non selected variables and their latest test statistics and p-values.
}
\item{mat}{
A matrix with the selected variables and their latest statistics and logged p-values.
}

\item{ci_test}{
The conditional independence test used. 
}
\item{final}{
The final regression model.
}
}

\author{
Michail Tsagris
 
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
  \code{\link{glm.fsreg}, \link{lm.fsreg}, \link{bic.fsreg}, \link{bic.glm.fsreg}, \link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)
dataset <- matrix( runif(500 * 10, 1, 100), ncol = 10 )
target <- rnorm(500)
a <- bs.reg(target, dataset, threshold = 0.05, test = "testIndRQ") 
b <- bs.reg(target, dataset, threshold = 0.05, test = "testIndReg") 
b2 <- bs.reg(target, dataset, threshold = 0.05, test = "testIndFisher") 
}

\keyword{ Backward regression }
\keyword{ Markov Blanket }
\keyword{ Variable Selection }
