\name{Conditional independence tests for left censored data }
\alias{testIndTobit}
\alias{waldTobit}
\alias{permTobit}
\title{
Conditional independence test for survival data 
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS. This test can based on the Cox (semi-parametric) regression or on the Weibull (parametric) regression.
}

\usage{
testIndTobit(target, dataset, xIndex, csIndex, wei = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL)

waldTobit(target, dataset, xIndex, csIndex, wei = NULL,  
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL)

permTobit(target, dataset, xIndex, csIndex, wei = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
threshold = 0.05, R = 999)
}

\arguments{
\item{target}{
A Survival object (class Surv from package survival) with left censored data. This test works with responses that are left censored. See the examples below (the command \link{Surv} or the final example in the \link{survreg} documentation of the "survival" package) for more information on how to create the target variable.  
}
\item{dataset}{
A numeric matrix or data frame, in case of categorical predictors (factors), containing the variables for performing the test. Rows as samples and columns as features. In the cases of "waldTobit" and "permTobit" this is strictly a matrix. 
}
\item{xIndex}{
The index of the variable whose association with the target we want to test.
}
\item{csIndex}{
The indices of the variables to condition on. If you have no variables set this equal to 0.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{univariateModels}{
Fast alternative to the hash object for univariate test. List with vectors "pvalues" (p-values), "stats" (statistics) and "flags" (flag = TRUE if the test was succesful) representing the univariate association of each variable with the target. Default value is NULL.
}
\item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use the hash-based implementation of the statistics of SES. Default value is FALSE. If TRUE you have to specify the stat_hash argument and the pvalue_hash argument.
}
\item{stat_hash}{
A hash object which contains the cached generated statistics of a SES run in the current dataset, using the current test.
}
\item{pvalue_hash}{
A hash object which contains the cached generated p-values of a SES run in the current dataset, using the current test.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. 
}
\item{R}{
The number of permutations, set to 999 by default. There is a trick to avoind doing all permutations. As soon as the number of times the permuted test statistic is more than the observed test statistic is more than 50 (in this example case), the p-value has exceeded the signifiance level (threshold value) and hence the predictor variable is not significant. There is no need to continue do the extra permutations, as a decision has already been made. 
}
}

\details{
Tobit regression is performed. The implied model is Gaussian with left censored data.

If hash = TRUE, censIndCR, censIndWR and censIndER require the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. 
These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order to speed up next runs of the current statistic test. 
If "SESoutput" is the output of a SES run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.

Important: Use these arguments only with the same dataset that was used at initialization.

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".

The log-likelihood ratio test used in "testIndTobit" requires the fitting of two models. The Wald test used in "waldTobit", requires fitting of only one model, the full one. The significance of the variable is examined only. Only continuous (or binary) predictor variables are currently accepted in this test. 
}

\value{
A list including:
\item{pvalue}{
A numeric value that represents the logarithm of the generated p-value.
}
\item{stat}{
A numeric value that represents the generated statistic.
}
\item{stat_hash}{
The current hash object used for the statistics. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
\item{pvalue_hash}{
The current hash object used for the p-values. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
}

\references{
Tobin James (1958). Estimation of relationships for limited dependent variables. Econometrica. 26(1): 24-36. 
}

\author{
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\note{
This test uses the functions "survreg" and Surv of the package survival and the function anova (analysis of variance) of the package stats.
}

\seealso{
\code{\link{SES}, \link{censIndWR}, \link{testIndFisher}, \link{gSquare}, \link{testIndLogistic}, \link{Surv}, \link{anova}, \link{CondIndTests}}
}

\examples{
require(survival, quietly = TRUE)

x <- matrix( rnorm(100 * 30), ncol = 30)
y <- x[, 1] - x[, 2] + rnorm(100, 5)
y[y < 0 ] <- 0
y <- survival::Surv(y, y>0, type = 'left') 
  
#run the censIndCR   conditional independence test
testIndTobit(y, x, xIndex = 12, csIndex = c(5, 7, 4) )
waldTobit(y, x, xIndex = 12, csIndex = c(5, 7, 4) )
permTobit(y, x, xIndex = 12, csIndex = c(5, 7, 4), R = 499 )
  
#run the SES algorithm using the censIndCR conditional independence
#test for the survival class variable
a <- MMPC(y, x, max_k = 2, threshold = 0.05, test = "testIndTobit");
}

\keyword{ Survival }
\keyword{ Tobit regression }
\keyword{ Log likelihood ratio }
\keyword{ Conditional independence test }
