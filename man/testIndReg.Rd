\name{Conditional independence tests for continous univariate and multivariate data }
\alias{testIndReg}
\alias{testIndRQ}
\alias{testIndMVreg}
\alias{testIndMMReg}
\alias{permReg}
\alias{permMMReg}
\alias{permRQ}
\alias{permMVreg}
\alias{waldMMReg}

\title{
Linear (and non-linear) regression conditional independence test for continous univariate and multivariate response variables
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS. The pvalue is calculated by comparing a linear regression model based on the conditioning set CS against a model whose regressor are both X and CS. The comparison is performed through an F test the appropriate degrees of freedom on the difference between the deviances of the two models.
}

\usage{
testIndReg(target, dataset, xIndex, csIndex, wei = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL)

testIndRQ(target, dataset, xIndex, csIndex, wei = NULL,
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL)

testIndMVreg(target, dataset, xIndex, csIndex, wei = NULL,
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL)

testIndMMReg(target, dataset, xIndex, csIndex, wei = NULL,  
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL)

waldMMReg(target, dataset, xIndex, csIndex, wei = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL) 

permReg(target, dataset, xIndex, csIndex, wei = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
threshold = 0.05, R = 999)

permMMReg(target, dataset, xIndex, csIndex, wei = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
threshold = 0.05, R = 999)

permRQ(target, dataset, xIndex, csIndex, wei = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
threshold = 0.05, R = 999)

permMVreg(target, dataset, xIndex, csIndex, wei = NULL,
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
threshold = 0.05, R = 999)
}

\arguments{
\item{target}{
A numeric vector containing the values of the target variable. If the values are proportions or percentages, i.e. strictly within 0 and 1 they are mapped into R using log( target/(1 - target) ). In the case of testIndMVreg, the same takes place true. See details for more information.
}
\item{dataset}{
A numeric matrix or data frame, in case of categorical predictors (factors), containing the variables for performing the test. Rows as samples and columns as features. In the cases of "waldIGreg" and "waldMMreg" this is strictly a matrix. 
}
\item{xIndex}{
The index of the variable whose association with the target we want to test.
}
\item{csIndex}{
The indices of the variables to condition on. If you have no variables set this equal to 0.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. They are not take into account in the robust regression via MM estimation (testIndMMReg and permMMreg).
}
\item{univariateModels}{
Fast alternative to the hash object for univariate test. List with vectors "pvalues" (p-values) and "stats" (statistics) representing the univariate association of each variable with the target. Default value is NULL.
}
\item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use tha hash-based implementation of the statistics of SES. Default value is FALSE. If TRUE you have to specify the stat_hash argument and the pvalue_hash argument.
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
The number of permutations, set to 999 by default. There is a trick to avoind doing all permutations. As soon as the number of times the permuted test statistic is more than the observed test statistic is more than 50 (if threshold = 0.05 and R = 999), the p-value has exceeded the signifiance level (threshold value) and hence the predictor variable is not significant. There is no need to continue do the extra permutations, as a decision has already been made. 
}
}

\details{
If hash = TRUE, all three tests require the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.

Important: Use these arguments only with the same dataset that was used at initialization.

TestIndReg offers linear regression. 

testIndMMReg offers robust linear MM estimation regression. 

TestIndRQ offers quantile (median) regression as a robust alternative to linear regression. 

In both cases, if the dependent variable consists of proportions (values between 0 and 1) the logit transformation is applied and the tests are applied then. 

testIndMVreg is for multivariate continuous response variables. Compositional data are positive multivariate data and each vector (observation) sums to the same constant, usually taken 1 for convenience. A check is performed and if such data are found, the additive log-ratio (multivariate logit) transformation (Aitchison, 1986) is applied beforehand. Zeros are not allowed. 

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".

The Wald test used in "waldMMReg" requires fitting of only one model, the full one. The significance of the variable is examined only. Only continuous (or binary) predictor variables are currently accepted in this test. 
}

\value{
A list including:
\item{pvalue}{
A numeric value that represents the logarithm of the generated p-value due to linear regression (see reference below).
}
\item{stat}{
A numeric value that represents the generated statistic due to linear regression (see reference below).
}
\item{stat_hash}{
The current hash object used for the statistics. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
\item{pvalue_hash}{
The current hash object used for the p-values. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
}

\references{
Draper, N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition.

Hampel F. R., Ronchetti E. M., Rousseeuw P. J., and Stahel W. A. (1986). Robust statistics: the approach based on influence functions. John Wiley & Sons. 

Koenker R.W. (2005). Quantile regression. New York, Cambridge University Press.

Yohai, V. J. (1987). High breakdown-point and high efficiency robust estimates for regression. The Annals of Statistics, 15(2): 642-656.

Mardia, Kanti, John T. Kent and John M. Bibby. Multivariate analysis. Academic press, 1979.

John Aitchison. The Statistical Analysis of Compositional Data, Chapman & Hall; reprinted in 2003, with additional material, by The Blackburn Press.
}

\author{
Vincenzo Lagani, Ioannis Tsamardinos, Michail Tsagris and Giorgos Athineou

R implementation and documentation: Giorgos Athineou <athineou@csd.gr>, Vincenzo Lagani <vlagani@csd.uoc.gr> and Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{testIndRQ}, \link{testIndFisher}, \link{testIndSpearman}, \link{CondIndTests}}
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(runif(100 * 50, 1, 100), ncol = 50 )
#the target feature is the last column of the dataset as a vector
target <- dataset[, 50]
dataset <- dataset[, -50]

testIndReg(target, dataset, xIndex = 44, csIndex = 10)
testIndMMReg(target, dataset, xIndex = 44, csIndex = 10)
testIndRQ(target, dataset, xIndex = 44, csIndex = 10)
testIndIGreg(target, dataset, xIndex = 44, csIndex = 10)

#run the MMPC algorithm using the testIndReg conditional independence test
m1 <- MMPC(target, dataset, max_k = 3, threshold = 0.05, test = "testIndReg")
m2 <- MMPC(target, dataset, max_k = 3, threshold = 0.05, test = "testIndRQ")
}

\keyword{ Linear regression test}
\keyword{ Conditional independence test }
