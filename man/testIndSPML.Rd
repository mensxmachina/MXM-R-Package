\name{Conditional independence test for circular data}
\alias{testIndSPML}
\title{
Circular regression conditional independence test for circular class dependent variables and continuous predictors.
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS. The pvalue is calculated by comparing a Beta regression model based on the conditioning set CS against a model whose regressor are both X and CS. The comparison is performed through a log-likelihood ratio test using circular regression with 2 degrees of freedom.
}

\usage{
testIndSPML(target, dataset, xIndex, csIndex, wei = NULL, univariateModels = NULL, 
hash = FALSE, stat_hash = NULL, pvalue_hash = NULL)
}

\arguments{
\item{target}{
A numerical vector with the data expressed in radians, or a 2 column matrix with the cosine and sine of the response (i.e. unit vector). See examples for more information on this. 
}
\item{dataset}{
A numeric matrix with the variables for performing the test. Rows as samples and columns as features. currently this test accepts only continuous variables.
}
\item{xIndex}{
The index of the variable whose association with the target we want to test.
}
\item{csIndex}{
The indices of the variables to condition on. If you have no variables set this equal to 0.
}
\item{wei}{
This is not used in this test.
}
\item{univariateModels}{
Fast alternative to the hash object for univariate test. List with vectors "pvalues" (p-values), "stats" (statistics) representing the univariate association of each variable with the target. Default value is NULL.
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
%\item{threshold}{
%Threshold (suitable values in (0, 1)) for assessing p-values significance. 
%}
%\item{R}{
%The number of permutations, set to 999 by default. There is a trick to avoind doing all permutations. As soon as the number of times the permuted test statistic is more than the observed test statistic is more than 50 (if threshold = 0.05 and R = 999), the p-value has exceeded the signifiance level (threshold value) and hence the predictor variable is not significant. There is no need to continue do the extra permutations, as a decision has already been made. 
%}
}

\details{
If hash = TRUE, testIndSPML requires the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.
  
Important: Use these arguments only with the same dataset that was used at initialization. For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".
  
The log-likelihood ratio test used in "testIndSPML" requires the fitting of two models.  The significance of the variable is examined and Only continuous (or binary) predictor variables are currently accepted in this test. 
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
  \code{\link{univregs}, \link{SES}, \link{testIndReg}, \link{CondIndTests}}
}

\examples{
y <- runif(100, - pi, pi)  ## suppose these are radians
x <- matrix( rnorm(100 * 3), ncol = 3)
testIndSPML(y1, x, csIndex = 1, xIndex = 2)
## alternatively
y1 <- cbind(cos(y), sin(y)) ## a matrix with two columns
testIndSPML(y1, x, csIndex = 1, xIndex = 2)
}

\keyword{Beta regression}
\keyword{ Conditional independence test }
