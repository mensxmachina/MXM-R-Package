\name{eBIC for many regression models}
\alias{ebic.regs}
\title{
eBIC for many regression models
}

\description{
eBIC for many regression models.
}

\usage{
ebic.regs(target, dataset, xIndex, csIndex,  gam = NULL, test = NULL, wei = NULL, 
ncores = 1)
}

\arguments{
\item{target}{
The target (dependent) variable. It must be a numerical vector, a factor or a Surv object. 
}
\item{dataset}{
The indendent variable(s). This can be a matrix or a dataframe with continuous only variables, a data frame with mixed or only categorical variables. 
}
\item{xIndex}{
The indices of the variables whose association with the target you want to test.
}
\item{csIndex}{
The index or indices of the variable(s) to condition on. If this is 0, the the function \code{\link{univregs}} will be called.
}
\item{gam}{
In case the method is chosen to be "eBIC" one can also specify the \eqn{gamma} parameter. The default value is "NULL", so that the value is automatically calculated.
}
\item{test}{
One of the following: \bold{testIndBeta}, \bold{testIndReg}, \bold{testIndLogistic}, \bold{testIndOrdinal}, 
\bold{testIndPois}, \bold{testIndZIP}, \bold{testIndNB}, \bold{testIndClogit}, \bold{testIndBinom}, \bold{testIndIGreg}, 
\bold{censIndCR}, \bold{censIndWR}, \bold{censIndER}, \bold{censIndLLR} \bold{testIndMultinom}, \bold{testIndTobit}, \bold{testIndSPML}, \bold{testIndGamma} or \bold{testIndNormLog}. 

Note that in all cases you must give the name of the test, without " ". 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
}

\details{
This function is more as a help function for MMPC, but it can also be called directly by the user. In some, one should specify the regression model to use and the function will perform all simple regressions, i.e. all regression models between the target and each of the variables in the dataset. The function does not check for zero variance columns, only the "univregs" and related functions do. 
}

\value{
A list including:
\item{stat}{
The value of the test statistic.
}
\item{pvalue}{
The logarithm of the p-value of the test. 
}
}

\references{
Chen J. and Chen Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. 
Biometrika, 95(3): 759-771.

Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.

McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.

}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{univregs}, \link{SES}, \link{MMPC}, \link{CondIndTests}, \link{reg.fit}, \link{ridge.reg} }
}

\examples{
y <- rpois(100, 15)
x <- matrix( rnorm(100 * 10), ncol = 10)
a1 <- univregs(y, x, test = testIndPois)
a2 <- perm.univregs(y, x, test = permPois)
a3 <- wald.univregs(y, x, test = waldPois)
a4 <- cond.regs(y, as.data.frame(x), xIndex = 1:4, csIndex = 5, test = testIndPois)
}

\keyword{ Regression models }
