\name{Backward selection regression for GLMM using the eBIC}
\alias{ebic.glmm.bsreg}
\title{
Backward selection regression for GLMM using the eBIC
}

\description{
Backward selection regression for GLMM using the eBIC
}

\usage{
ebic.glmm.bsreg(target, dataset, id, wei = NULL, gam = NULL, test = "testIndGLMMReg") 
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
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{gam}{
In case the method is chosen to be "eBIC" one can also specify the \eqn{gamma} parameter. The default value is "NULL", so that the value is automatically calculated.
}
\item{test}{
This is for the type of regression to be used, "testIndGLMMReg", for Gaussian regression, "testIndGLMMLogistic for logistic regression or "testIndGLMMPois" for Poisson regression.
}
%\item{ncores}{
%How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
%}
}

\details{
The algorithm is a variation of the usual forward selection. At every step, the most significant variable enters the selected variables set. 
In addition, only the significant variables stay and are further examined. The non signifcant ones are dropped. This goes until no variable 
can enter the set. The user has the option to redo this step 1 or more times (the argument K). In the end, a backward selection is performed
to remove falsely selected variables.
}

\value{
A list including:
\item{runtime}{
The runtime required.
}
\item{info}{
A matrix with the number of variables and the number of tests performed (or models fitted) at each round (value of K). 
}
\item{mat}{
A matrix with the selected variables and their eBIC.
}
}

\references{
Borboudakis G. and Tsamardinos I. (2019). Forward-backward selection with early dropping. Journal of Machine Learning Research, 20(8): 1-39.

Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.
}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{ \link{fbed.glmm.reg}, \link{glmm.bsreg}, \link{MMPC.glmm} }
}

\examples{
\dontrun{
require(lme4)
data(sleepstudy)
reaction <- sleepstudy$Reaction
days <- sleepstudy$Days
subject <- sleepstudy$Subject
x <- matrix(rnorm(180 * 20),ncol = 20) ## unrelated preidctor variables
m1 <- ebic.glmm.bsreg(reaction, x, id = subject) 
m2 <- MMPC.glmm(reaction, group = subject, dataset = x)
}
}

