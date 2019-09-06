\name{Correlation based conditonal independence tests}
\alias{testIndFisher}
\alias{testIndSpearman}
\alias{testIndMMFisher}
\alias{permFisher}
\alias{permMMFisher}
\alias{permDcor}
\title{
Fisher and Spearman conditional independence test for continuous class variables
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS.
}

\usage{
testIndFisher(target, dataset, xIndex, csIndex, wei = NULL, statistic = FALSE, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, 
pvalue_hash = NULL)

testIndMMFisher(target, dataset, xIndex, csIndex, wei = NULL, statistic = FALSE, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, 
pvalue_hash = NULL)

testIndSpearman(target, dataset, xIndex, csIndex, wei = NULL, statistic = FALSE, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, 
pvalue_hash = NULL)

permFisher(target, dataset, xIndex, csIndex, wei = NULL, statistic = FALSE, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, 
pvalue_hash = NULL, threshold = 0.05, R = 999)

permMMFisher(target, dataset, xIndex, csIndex, wei = NULL, statistic = FALSE, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, 
pvalue_hash = NULL, threshold = 0.05, R = 999)

permDcor(target, dataset, xIndex, csIndex, wei = NULL, statistic = FALSE, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, 
pvalue_hash = NULL, threshold = 0.05, R = 499)
}

\arguments{
\item{target}{
A numeric vector containing the values of the target variable. If the values are proportions or percentages, i.e. strictly within 0 and 1 they are mapped into R using log( target/(1 - target) ). 
This can also be a list of vectors as well. In this case, the metanalytic approach is used. 
}
\item{dataset}{
A numeric matrix containing the variables for performing the test. Rows as samples and columns as features.
}
\item{xIndex}{
The index of the variable whose association with the target we want to test.
}
\item{csIndex}{
The indices of the variables to condition on. If you have no variables set this equal to 0.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. This is not used with "permDcor".
}
\item{statistic}{ 
A boolean variable indicating whether the test statistics (TRUE) or the p-values should be combined (FALSE). See the details about this.
For the permFisher test this is not taken into account.
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
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05. This is actually obsolete here, but has to be in order tyo have a concise list of input arguments across the same family of functions.
}
\item{R}{
The number of permutations to use. The default value is 999. For the "permDcor" this is set to 499. 
}
}

\details{
If hash = TRUE, testIndFisher requires the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.

Important: Use these arguments only with the same dataset that was used at initialization.

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".

Note that if the \code{\link{testIndReg}} is used instead the results will not be be the same, unless the sample size is very large. This is because the Fisher test uses the t distribution stemming from the Fisher's z transform and not the t distribution of the correlation coefficient.

BE CAREFUL with testIndSpearman. The Pearson's correlation coefficient is actually calculated. So, you must have transformed the data into their ranks before plugging them here. The reason for this is to speed up the computation time, as this test can be used in SES, MMPC and mmhc.skel. The variance of the Fisher transformed Spearman's correlation is \eqn{\frac{1.06}{n-3}} and the variance of the Fisher transformed Pearson's correlation coefficient is \eqn{\frac{1}{n-3}}.

When performing the above tests with multiple datasets, the test statistic and the p-values are combined in a meta-analytic way. Is up to the user to decide whether to use the fixed effects model approach and combine the test statistics (statistic = TRUE), or combine the p-values as Fisher suggested (statistic = FALSE). 

The argument R is useful only for the permFisher and permDcor tests. The permDcor test uses the distance correlation instead of the usual Pearson or Spearman correlations.

TestIndMMFisher does a robust estimation of the correlation via MM regression.
}

\value{
A list including:
\item{pvalue}{
A numeric value that represents the logarithm of the generated p-value due to Fisher's method (see reference below).
}
\item{stat}{
A numeric value that represents the generated statistic due to Fisher's method (see reference below).
}
\item{stat_hash}{
The current hash object used for the statistics. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
\item{pvalue_hash}{
The current hash object used for the p-values. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
}

\references{
Fisher R. A. (1925). Statistical methods for research workers. Genesis Publishing Pvt Ltd.

Fisher R. A. (1948). Combining independent tests of significance. American Statistician, 2(5), 30--31

Fisher R. A. (1915). Frequency distribution of the values of the correlation coefficient in samples from an indefinitely large population. Biometrika, 10(4): 507--521.

Fieller E. C., Hartley H. O. and Pearson E. S. (1957). Tests for rank correlation coefficients. I. Biometrika, 44(3/4): 470--481.

Fieller E. C. and Pearson E. S. (1961). Tests for rank correlation coefficients. II. Biometrika, 48(1/2): 29--40.

Hampel F. R., Ronchetti E. M., Rousseeuw P. J., and Stahel W. A. (1986). Robust statistics: the approach based on influence functions. John Wiley & Sons. 

Pearson, K. (1895). Note on regression and inheritance in the case of two parents. Proceedings of the Royal Society of London, 58, 240--242.

Peter Spirtes, Clark Glymour, and Richard Scheines. Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, second edition, January 2001.

Lee Rodgers J., and Nicewander W.A. (1988). "Thirteen ways to look at the correlation coefficient." The American Statistician 42(1): 59--66.

Shevlyakov G. and Smirnov P. (2011). Robust Estimation of the Correlation Coefficient: An Attempt of Survey. Austrian Journal of Statistics, 40(1 & 2): 147--156.

Szekely G.J. and Rizzo, M.L. (2014). Partial distance correlation with methods for dissimilarities. The Annals of Statistics, 42(6): 2382--2412.

Szekely G.J. and Rizzo M.L. (2013).  Energy statistics: A class of statistics based on distances. Journal of Statistical Planning and Inference 143(8): 1249--1272.
}

\author{
Vincenzo Lagani and Ioannis Tsamardinos

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Vincenzo Lagani <vlagani@csd.uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{testIndSpearman}, \link{testIndReg}, \link{SES}, \link{testIndLogistic}, \link{gSquare}, \link{CondIndTests}}
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(runif(300 * 50, 1, 1000), nrow = 50 )
#the target feature is the last column of the dataset as a vector
target <- dataset[, 50]
res1 <- testIndFisher(target, dataset, xIndex = 44, csIndex = 10)
res2 <- testIndSpearman(target, dataset, xIndex = 44, csIndex = 10)
res3 <- permFisher(target, dataset, xIndex = 44, csIndex = 10, R = 999)
res4 <- permDcor(target, dataset, xIndex = 44, csIndex = 10, R = 99)

#define class variable (here tha last column of the dataset)
dataset <- dataset[, -50]
#run the MMPC algorithm using the testIndFisher conditional independence test
mmpcObject <- MMPC(target, dataset, max_k = 3, threshold = 0.05, test = "testIndFisher")
}

\keyword{ Fisher's Test}
\keyword{ Conditional Independence Test }
