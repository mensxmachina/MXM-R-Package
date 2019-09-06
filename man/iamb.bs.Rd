\name{IAMB backward selection phase}
\alias{iamb.bs}
\title{
IAMB backward selection phase
}

\description{
IAMB backward selection phase.
}

\usage{
iamb.bs(target, dataset, threshold = 0.05, wei = NULL, test = NULL, user_test = NULL)
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. 
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = observations). In either case, only two cases are avaialble, either all data are continuous, or categorical. 
}
\item{threshold}{
Threshold (suitable values in (0,1)) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
The regression model to use. Available options are most of the tests for SES and MMPC. The ones NOT available are "gSquare", "censIndER", "testIndMVreg", "testIndClogit", "testIndSpearman" and "testIndFisher".
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when robust is set to TRUE.
}
\item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, the "test" argument is ignored.
}
}
\details{
IAMB stands for Incremental Association Markov Blanket. The algorithm comprises of a forward selection and a modified backward selection process. This functions does the modified backward selection process. In the usual backward selection, among the non singificant variabels, the one with the maximum p-value is dropped. So, one variable is removed at every step. In the IAMB backward phase, at aevery step, all non significant variables are removed. This makes it a lot faster. 
}

\value{
The output of the algorithm is a list of an S3 object including:
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{ci_test}{
The conditional independence test used. 
}
\item{vars}{
The selected variables.
}
\item{mat}{
A matrix with the selected variables and their latest test statistic and p-value. If no variable is selected this is NULL.
}
\item{final}{
The final regression model.
}
}

\references{
Tsamardinos, I., Aliferis, C. F., Statnikov, A. R., & Statnikov, E. (2003). Algorithms for Large Scale Markov Blanket Discovery. In FLAIRS conference, pp. 376-380.
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

a1 <- iamb.bs(target, dataset, threshold = 0.05, test = "testIndRQ") 
a2 <- bs.reg(target, dataset, threshold = 0.05, test = "testIndRQ") 
}

\keyword{ Backward regression }
\keyword{ Markov Blanket }
\keyword{ Variable Selection }
