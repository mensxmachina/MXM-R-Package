\name{G-square conditional independence test for discrete data}
\alias{gSquare}
\alias{permgSquare}
\title{
G-square conditional independence test for discrete data
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent 
from 'TARGET' given a conditioning set CS. This test is based on the log likelihood ratio test.
}

\usage{
gSquare(target, dataset, xIndex, csIndex, wei = NULL,  
univariateModels = NULL, hash = FALSE, stat_hash = NULL, 
pvalue_hash = NULL)

permgSquare(target, dataset, xIndex, csIndex, wei = NULL,
univariateModels = NULL, hash = FALSE, stat_hash = NULL, 
pvalue_hash = NULL, threshold = 0.05, R = 999)
}

\arguments{
\item{target}{
A numeric vector containing the values of the target variable. \bold{The minimum value must be 0}.
}
\item{dataset}{
A numeric matrix containing the variables for performing the test. Rows as samples and columns as features. \bold{The minimum value must be 0}.
}
\item{xIndex}{
The index of the variable whose association with the target we want to test.
}
\item{csIndex}{
The indices of the variables to condition on.
}
\item{wei}{
This argument is not used in this test. 
}
\item{univariateModels}{
Fast alternative to the hash object for univariate tests. List with vectors "pvalues" (p-values), "stats" (statistics) and 
"flags" (flag = TRUE if the test was succesful) representing the univariate association of each variable with the target. 
Default value is NULL.
}
\item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use the hash-based implementation of the statistics of SES. 
Default value is FALSE. If TRUE you have to specify the stat_hash argument and the pvalue_hash argument.
}
\item{stat_hash}{
A hash object which contains the cached generated statistics of a SES run in the current dataset, 
using the current test.
}
\item{pvalue_hash}{
A hash object which contains the cached generated p-values of a SES run in the current dataset, 
using the current test.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05. This is actually obsolete here, but has to be in order to have a concise list of input arguments across the same family of functions.
}
\item{R}{
The number of permutations to use. The default value is 999.
}
}

\details{
If the number of samples is at least 5 times the number of the parameters to be estimated, the test is performed, otherwise, 
independence is not rejected (see Tsmardinos et al., 2006, pg. 43) 

If hash = TRUE, testIndLogistic requires the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the 
statistical test. These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order 
to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES run, then these objects can be 
retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.

Important: Use these arguments only with the same dataset that was used at initialization.

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".
}

\value{
A list including:
\item{pvalue}{
A numeric value that represents the logarithm of the generated p-value of the \eqn{G^2} test (see reference below).
}
\item{stat}{
A numeric value that represents the generated statistic of the \eqn{G^2} test (see reference below).
}
\item{stat_hash}{
The current hash object used for the statistics. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
\item{pvalue_hash}{
The current hash object used for the p-values. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
}

\references{
Tsamardinos, Ioannis, Laura E. Brown, and Constantin F. Aliferis. The max-min hill-climbing Bayesian network structure learning algorithm. 
Machine learning, 2006 65(1): 31--78. 
}

\author{
R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%}

\seealso{
\code{\link{SES}, \link{testIndFisher}, \link{testIndLogistic}, \link{censIndCR}, \link{CondIndTests}}
}

\examples{
#simulate a dataset with binary data
dataset <- matrix(rbinom(500 * 51, 1, 0.6), ncol = 51)
#initialize binary target
target <- dataset[, 51]
#remove target from the dataset
dataset <- dataset[, -51]

#run the gSquare conditional independence test for the binary class variable
results <- gSquare(target, dataset, xIndex = 44, csIndex = c(10,20) )
results

#run SES algorithm using the gSquare conditional independence test for the binary class variable
sesObject <- SES(target, dataset, max_k = 3, threshold = 0.05, test = "gSquare");
target <- as.factor(target)
sesObject2 <- SES(target, dataset, max_k = 3, threshold = 0.05, test = "testIndLogistic");
}

\keyword{ G-square test }
\keyword{ Log likelihood ratio }
\keyword{ Conditional independence test }
