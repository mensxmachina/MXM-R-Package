\name{Conditional independence regression based tests}
\alias{cond.regs}
\alias{glmm.condregs}
\alias{gee.condregs}
\title{
Conditional independence regression based tests
}

\description{
Conditional independence regression based tests.
}

\usage{
cond.regs(target, dataset, xIndex, csIndex, test = NULL, wei = NULL, ncores = 1)

glmm.condregs(target, reps = NULL, id, dataset, xIndex, csIndex, test, wei = NULL, 
slopes = FALSE, ncores = 1)

gee.condregs(target, reps = NULL, id, dataset, xIndex, csIndex, test, wei = NULL, 
correl = "echangeable", se = "jack", ncores = 1)     
}

\arguments{
\item{target}{
The target (dependent) variable. It must be a numerical vector. 
}
\item{dataset}{
The indendent variable(s). For the "univregs" this can be a matrix or a dataframe with continuous only variables, a data frame with mixed or only categorical variables. For the "wald.univregs", "perm.univregs", "score.univregs" and "rint.regs" this can only by a numerical matrix.
}
\item{xIndex}{
The indices of the variables whose association with the target you want to test.
}
\item{csIndex}{
The index or indices of the variable(s) to condition on. If this is 0, the the function \code{\link{univregs}} will be called.
}
\item{test}{
For the "cond.regs" one of the following: \bold{testIndBeta}, \bold{testIndReg}, \bold{testIndLogistic}, \bold{testIndOrdinal}, \bold{testIndPois}, \bold{testIndQPois}, \bold{testIndZIP}, \bold{testIndNB}, \bold{testIndClogit}, \bold{testIndBinom}, \bold{testIndQBinom}, \bold{testIndIGreg}, \bold{censIndCR}, \bold{censIndWR}, \bold{censIndER}, \bold{censIndLLR}, \bold{testIndMMReg}, \bold{testIndMVreg}, \bold{testIndMultinom}, \bold{testIndTobit}, \bold{testIndGamma}, \bold{testIndNormLog} or \bold{testIndSPML} for a circular target. 

For the "glmm.condregs" one of the following: \bold{testIndLMReg}, \bold{testIndGLMMLogistic}, 
\bold{testIndGLMMPois}, \bold{testIndGLMMOrdinal},  \bold{testIndGLMMGamma} or \bold{testIndGLMMNormLog}. 

For the "gee.condregs" one of the following: \bold{testIndGEEReg}, \bold{testIndGEELogistic}, 
\bold{testIndGEEPois}, \bold{testIndGEEGamma} or \bold{testIndGEENormLog}. 

Note that in all cases you must give the name of the test, without " ". 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample 
sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other 
cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the 
first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction 
in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
\item{id}{
A numerical vector of the same length as target with integer valued numbers, such as 1, 2, 3,... (zeros, negative values and 
factors are not allowed) specifying the clusters or subjects. This argument is for the rint.regs (see details for more information).
}
\item{reps}{
If you have measurements over time (lognitudinal data) you can put the time here (the length must be equal to the length of the target) 
or set it equal to NULL. (see details for more information).
}
\item{slopes}{
Should random slopes for the ime effect be fitted as well? Default value is FALSE. 
}
\item{correl}{
The correlation structure. For the Gaussian, Logistic, Poisson and Gamma regression this can be either "exchangeable" 
(compound symmetry, suitable for clustered data) or "ar1" (AR(1) model, suitable for longitudinal data). 
}
\item{se}{
The method for estimating standard errors. This is very important and crucial. The available options for Gaussian, Logistic, 
Poisson and Gamma regression are: a)  'san.se', the usual robust estimate. b) 'jack': if approximate jackknife variance estimate 
should be computed. c) 'j1s': if 1-step jackknife variance estimate should be computed and d) 'fij': logical indicating if fully 
iterated jackknife variance estimate should be computed. If you have many clusters (sets of repeated measurements) "san.se" is 
fine as it is astmpotically correct, plus jacknife estimates will take longer. If you have a few clusters, then maybe 
it's better to use jacknife estimates. 

The jackknife variance estimator was suggested by Paik (1988), which is quite suitable for cases when the number of subjects is 
small (K < 30), as in many biological studies. The simulation studies conducted by Ziegler et al. (2000) and Yan and Fine (2004) 
showed that the approximate jackknife estimates are in many cases in good agreement with the fully iterated ones. 
}
}

\details{
This function is more as a help function for MMPC, but it can also be called directly by the user. In some, one should specify 
the regression model to use and the function will perform all simple regressions, i.e. all regression models between the target 
and each of the variables in the dataset. The function does not check for zero variance columns, only the "univregs" and 
related functions do. 

If you want to use the GEE methodology, make sure you load the library geepack first. 
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

Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data. Journal of the American Statistical Association, 93(443): 1068-1077.
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
a4 <- cond.regs(y, as.data.frame(x), xIndex = 1:9, csIndex = 10, test = testIndPois)
}

\keyword{ Regression models }
