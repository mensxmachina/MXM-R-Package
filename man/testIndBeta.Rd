\name{Conditional independence test for proportions/percentages}
\alias{testIndBeta}
\alias{permBeta}
\alias{waldBeta}
\title{
Beta regression conditional independence test for proportions/percentage class dependent variables and mixed predictors
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS. The pvalue is calculated by comparing a Beta regression model based on the conditioning set CS against a model whose regressor are both X and CS. The comparison is performed through a chi-square test with the appropriate degrees of freedom on the difference between the deviances of the two models.
}

\usage{
testIndBeta(target, dataset, xIndex, csIndex, wei = NULL,  
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL)

permBeta(target, dataset, xIndex, csIndex, wei = NULL,  
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
threshold = 0.05, R = 999)

waldBeta(target, dataset, xIndex, csIndex, wei = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL)
}

\arguments{
\item{target}{
A numeric vector containing the values of the target variable. They must be percentages or proportions, i.e. within the (0, 1) interval. Currently 0 and/or 1 values are not allowed.
}
\item{dataset}{
A numeric matrix or data frame, in case of categorical predictors (factors), containing the variables for performing the test. Rows as samples and columns as features. In the case of "waldBeta" this is strictly a matrix. 
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
If hash = TRUE, testIndBeta requires the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.

Important: Use these arguments only with the same dataset that was used at initialization. For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests". An alternative regression to this is "testIndReg" and "testIndRQ" .In these two latter cases, the logit transformation is first applied to the target variable. 

The log-likelihood ratio test used in "testIndBeta" requires the fitting of two models. The Wald test used in "waldBeta" requires fitting of only one model, the full one. The significance of the variable is examined only. Only continuous (or binary) predictor variables are currently accepted in this test. 

The \code{\link{testIndQBinom}} is another alternative to use.
}

\value{
A list including:
\item{pvalue}{
A numeric value that represents the logarithm of the generated p-value due to Beta regression (see reference below).
}
\item{stat}{
A numeric value that represents the generated statistic due to Beta regression (see reference below).
}
\item{stat_hash}{
The current hash object used for the statistics. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
\item{pvalue_hash}{
The current hash object used for the logged p-values. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
}

\references{
Ferrari S.L.P. and Cribari-Neto F. (2004). Beta Regression for Modelling Rates and Proportions. Journal of Applied Statistics, 31(7): 799-815.
}

\author{
Vincenzo Lagani, Ioannis Tsamardinos, Michail Tsagris and Giorgos Athineou

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr>, Vincenzo Lagani <vlagani@csd.uoc.gr> and Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{testIndQBinom}, \link{SES}, \link{testIndReg}, \link{testIndRQ}, \link{testIndFisher}, \link{CondIndTests}}
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(runif(100 * 10, 1, 1000), ncol = 10 ) 
#the target feature is the last column of the dataset as a vector
target <- dataset[, 10]
dataset <- dataset[, -10]
target <- target / (max(target) + 2 )

testIndBeta(target, dataset, xIndex = 4, csIndex = 9)

#run the MMPC algorithm using the testIndBeta conditional independence test
mmpcObject <- MMPC(target, dataset, max_k = 3, threshold = 0.05, test = "testIndBeta")
}

\keyword{Beta regression}
\keyword{ Conditional independence test }
