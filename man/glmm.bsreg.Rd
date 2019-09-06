\name{Backward selection regression for GLMM}
\alias{glmm.bsreg}
\title{
Backward selection regression for GLMM
}

\description{
Backward selection regression for GLMM
}

\usage{
glmm.bsreg(target, dataset, id, threshold = 0.05, wei = NULL, test = "testIndGLMMReg") 
}

\arguments{
\item{target}{
The class variable. This can be a numerical vector with continuous data, binary or discrete valued data. It can also be a factor variable with two levels only.
}
\item{dataset}{
The dataset; provide a numerical a matrix (columns = variables, rows = samples).
}
\item{id}{
This is a numerical vector of the same size as target denoting the groups or the subjects. 
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL.
}
\item{test}{
This is for the type of regression to be used, "testIndGLMMReg", for Gaussian regression, "testIndGLMMLogistic for logistic regression or "testIndGLMMPois" for Poisson regression.
}
%\item{ncores}{
%How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
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
A matrix with the non selected variables and their latest test statistics and logged p-values.
}
\item{mat}{
A matrix with the selected variables and their latest statistics and p-values.
}
\item{final}{
The final regression model.
}
}

\references{
Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.
}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{ \link{fbed.glmm.reg}, \link{ebic.glmm.bsreg}, \link{MMPC.glmm} }
}

\examples{
\dontrun{
require(lme4)
data(sleepstudy)
reaction <- sleepstudy$Reaction
days <- sleepstudy$Days
subject <- sleepstudy$Subject
x <- matrix(rnorm(180 * 10),ncol = 10) ## unrelated preidctor variables
m1 <- glmm.bsreg(Reaction, x, subject) 
m2 <- MMPC.glmm(target = reaction, group = subject, dataset = x)
}
}

