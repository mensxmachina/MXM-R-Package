\name{Univariate regression based tests}
\alias{univregs}
\alias{ebic.univregs} 
\alias{perm.univregs}
\alias{wald.univregs}
\alias{score.univregs}
\alias{big.score.univregs}
\alias{rint.regs}
\alias{glmm.univregs}
\alias{gee.univregs}

\title{
Univariate regression based tests
}

\description{
Univariate regression based tests.
}

\usage{
univregs(target, dataset, targetID = -1, test = NULL, user_test = NULL, wei = NULL, 
ncores = 1) 

ebic.univregs(target, dataset, targetID = -1, test = NULL, user_test = NULL, 
wei = NULL, ncores = 1, gam = NULL) 

wald.univregs(target, dataset, targetID = - 1, test = NULL, user_test = NULL, 
wei = NULL, ncores = 1)

perm.univregs(target, dataset, targetID = -1, test = NULL, user_test = NULL, 
wei = NULL, threshold = 0.05, R = 999, ncores = 1) 

score.univregs(target, dataset, test) 

big.score.univregs(target = NULL, dataset, test) 

rint.regs(target, dataset, targetID = -1, id, reps = NULL, tol = 1e-08)

glmm.univregs(target, reps = NULL, id, dataset, targetID = -1, test, wei = NULL, 
slopes = FALSE, ncores = 1)

gee.univregs(target, reps = NULL, id, dataset, targetID = -1, test, wei = NULL, 
correl = "echangeable", se = "jack", ncores = 1)
}

\arguments{
\item{target}{
The target (dependent) variable. It must be a numerical vector. In the case of "big.score.univregs" this can also be NULL if it is included in the first column of the dataset.
}
\item{dataset}{
The indendent variable(s). For the "univregs" this can be a matrix or a dataframe with continuous only variables, a data frame with mixed or only categorical variables. For the "wald.univregs", "perm.univregs", "score.univregs" and "rint.regs" this can only by a numerical matrix. For the "big.score.univregs" this is a big.matrix object with continuous data only.
}
\item{targetID}{
This is by default -1. If the target is not a variable but is in included in the dataset, you can specify the column of dataset playing the role of the target. 
}
\item{test}{
For the "univregs" this can only be one of the following: \bold{testIndFisher}, \bold{testIndSpearman}, \bold{gSquare}, \bold{testIndBeta}, \bold{testIndReg}, \bold{testIndLogistic}, \bold{testIndMultinom}, \bold{testIndOrdinal}, 
\bold{testIndPois}, \bold{testIndZIP}, \bold{testIndNB}, \bold{testIndClogit}, \bold{testIndBinom}, \bold{testIndIGreg}, \bold{censIndCR}, \bold{censIndWR}, \bold{censIndER}, \bold{censIndLLR}, \bold{testIndTobit}, \bold{testIndGamma}, \bold{testIndNormLog} or \bold{testIndSPML} for a circular target. For the \bold{testIndSpearman} the user must supply the ranked target and 
ranked dataset. The reason for this, is because this function is called internally by SES and MMPC, the ranks have 
already been applied and there is no reason to re-rank the data. 

For the "ebic.univregs" this can only be one of the following: \bold{testIndFisher}, \bold{testIndBeta}, \bold{testIndReg}, \bold{testIndLogistic}, \bold{testIndMultinom}, \bold{testIndOrdinal}, \bold{testIndPois}, \bold{testIndZIP}, \bold{testIndNB}, \bold{testIndClogit}, \bold{testIndBinom}, \bold{censIndCR}, \bold{censIndWR}, \bold{censIndER}, \bold{censIndLLR}, \bold{testIndTobit}, \bold{testIndGamma}, \bold{testIndNormLog} or \bold{testIndSPML} for a circular target.

For the "wald.univregs" this can only be one of the following: \bold{waldBeta}, \bold{waldCR}, \bold{waldWR}, 
\bold{waldER}, \bold{waldLLR}, \bold{waldTobit}, \bold{waldClogit}, \bold{waldLogistic}, \bold{waldPois}, \bold{waldNB}, \bold{waldBinom}, \bold{waldZIP}, \bold{waldMMReg}, \bold{waldIGreg}, \bold{waldOrdinal}, \bold{waldGamma} or \bold{waldNormLog}.

For the "perm.univregs" this can only be one of the following: \bold{permFisher}, \bold{permReg}, \bold{permRQ}, 
\bold{permBeta}, \bold{permCR}, \bold{permWR}, \bold{permER}, \bold{permLLR}, \bold{permTobit}, \bold{permClogit}, \bold{permLogistic}, \bold{permPois}, \bold{permNB}, \bold{permBinom}, \bold{permgSquare}, \bold{permZIP}, \bold{permMVreg}, \bold{permIGreg}, \bold{permGamma} or \bold{permNormLog}.

For the "score.univregs" and "big.score.univregs" this can only be one of the following: \bold{testIndBeta}, 
\bold{testIndLogistic}, \bold{testIndMultinom}, \bold{testIndPois}, \bold{testIndNB}, \bold{testIndGamma} or 
\bold{censIndWR} but with no censored values and simply the vector, 
not a Surv object. Ordinal regression is not supported. 

For the "glmm.univregs" this can only be one of the following: \bold{testIndGLMMReg}, \bold{testIndLMM}, 
\bold{testIndGLMMLogistic}, \bold{testIndGLMMOrdinal}, \bold{testIndGLMMPois}, \bold{testIndGLMMGamma}, 
\bold{testIndGLMMNormLog} or \bold{testIndGLMMCR}. 

For the "gee.univregs" this can only be one of the following: \bold{testIndGLMMReg}, \bold{testIndGLMMLogistic}, 
\bold{testIndGLMMOrdinal}, \bold{testIndGLMMPois}, \bold{testIndGLMMGamma} or \bold{testIndGLMMNormLog}. 

\bold{Note that in all cases you must give the name of the test, without " "}. 
}
\item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, the "test" 
argument is ignored.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{gam}{
This is for ebic.univregs only and it refers to the \eqn{gamma} value of the eBIC. If it is NULL, the default value is calculated and if gam is zero, the usual BIC is returned. see Chen and Chen (2008) for more information.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05. The reason for this is to speed up the computations. 
The permutation p-value is calculated as the proportion of times the permuted test statistic is higher than the observed test statistic. When running the
permutations, if the proportion is more than 0.05, the process stops. A decision must be made fast, and if the non rejection decision has been made, there 
is no need to perform the rest permutations; the decision cannot change.
}
\item{R}{
This is for the permutations based regression models. It is the number of permutations to apply. Note, that not all the number of permutations need be performed. If the number of times the test statistic is greater than the observed test statistic is more than threshold * R, the iterations stop, as a decision has already been made. 
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
\item{id}{
A numerical vector of the same length as target with integer valued numbers, such as 1, 2, 3,... (zeros, negative values and factors are not allowed) specifying the clusters or subjects. This argument is for the rint.regs (see details for more information).
}
\item{reps}{
If you have measurements over time (lognitudinal data) you can put the time here (the length must be equal to the length of the target) or set it equal to NULL. (see details for more information).
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm in the random effects models.
}
\item{slopes}{
Should random slopes for the ime effect be fitted as well? Default value is FALSE. 
}
\item{correl}{
The correlation structure. For the Gaussian, Logistic, Poisson and Gamma regression this can be either "exchangeable" (compound symmetry, suitable for clustered data) or "ar1" (AR(1) model, suitable for longitudinal data). For the ordinal logistic regression its only the "exchangeable" correlation sturcture.  
}
\item{se}{
The method for estimating standard errors. This is very important and crucial. The available options for Gaussian, Logistic, Poisson and Gamma regression are: a)  'san.se', the usual robust estimate. b) 'jack': if approximate jackknife variance estimate should be computed. 
c) 'j1s': if 1-step jackknife variance estimate should be computed and d) 'fij': logical indicating if fully iterated jackknife variance estimate should be computed. If you have many clusters (sets of repeated measurements) "san.se" is fine as it is astmpotically correct, plus jacknife estimates will take longer. If you have a few clusters, then maybe it's better to use jacknife estimates. 

The jackknife variance estimator was suggested by Paik (1988), which is quite suitable for cases when the number of subjects is small (K < 30), as in many biological studies. The simulation studies conducted by Ziegler et al. (2000) and Yan and Fine (2004) showed that the approximate jackknife estimates are in many cases in good agreement with the fully iterated ones. This is obsolete for "testIndGEEOrdinal", but is here for compatibility reasons.
}
}

\details{
This function is more as a help function for SES and MMPC, but it can also be called directly by the user. In some, one should specify the regression model to use and the function will perform all simple regressions, i.e. all regression models between the target and each of the variables in the dataset. 

For the score.univregs, score based tests are used which are extremely fast. 

For the rint.regs, we perform linear mixed models (weights are not allowed) with random intercepts only (no ranodm slopes). This function works for clustered or longitudinal data. The covariance structure we impose is compound symmetry, hence for longitudinal data, this may not be the best option, yet it will work. 

If you want to use the GEE methodology, make sure you load the library geepack first. 
}

\value{
In the case of "ebic.univregs" a list with one element 
\item{ebic}{ 
The eBIc of every model. If the i-th value is "Inf" it means that the i-th variable had zero variance. It had the same value everywhere.
}
For all other cases a list including:
\item{stat}{
The value of the test statistic. If the i-th value is zero (0) it means that the i-th variable had zero variance. It had the same value everywhere.
}
\item{pvalue}{
The logarithm of the p-value of the test. If the i-th value is zero (0) it means that the i-th variable had zero variance. It had the same value everywhere.
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
\code{ \link{cond.regs}, \link{SES}, \link{MMPC}, \link{CondIndTests}, \link{reg.fit}, \link{ridge.reg} }
}

\examples{
y <- rpois(50, 15)
x <- matrix( rnorm(50 * 7), ncol = 7)
a1 <- univregs(y, x, test = testIndPois)
a2 <- perm.univregs(y, x, test = permPois)
a3 <- wald.univregs(y, x, test = waldPois)
}

\keyword{ Regression models }
