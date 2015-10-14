\name{gSquare}
\alias{gSquare}
\title{
G square conditional independence test for discrete data.
}
\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS. This test is based on the log likelihood ratio test.
}
\usage{
gSquare(target, dataset, xIndex, csIndex, dataInfo = NULL, univariateModels = NULL,
hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, robust = FALSE)
}
\arguments{
  \item{target}{
  A numeric vector containing the values of the target variable.
}
  \item{dataset}{
  A numeric matrix containing the variables for performing the test. Rows as samples and columns as features.
}
  \item{xIndex}{
  The index of the variable whose association with the target we want to test.
}
  \item{csIndex}{
  The indices of the variables to condition on.
}
  \item{dataInfo}{
  list object with information on the structure of the data. Default value is NULL.
}
  \item{univariateModels}{
  Fast alternative to the hash object for univariate test. List with vectors "pvalues" (p-values), "stats" (statistics) and "flags" (flag = TRUE if the test was succesful) representing the univariate association of each variable with the target. Default value is NULL.
}
  \item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use the hash-based implementation of the statistics of SES. Default value is FALSE. If TRUE you have to specify the stat_hash argument and the pvalue_hash argument.
}
  \item{stat_hash}{
A hash object (hash package required) which contains the cached generated statistics of a SES run in the current dataset, using the current test.
}
  \item{pvalue_hash}{
A hash object (hash package required) which contains the cached generated p-values of a SES run in the current dataset, using the current test.
}
  \item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) 
to use a robust test. Currently it is not available for this test.
}
}
\details{
If hash = TRUE, testIndLogistic requires the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.

Important: Use these arguments only with the same dataset that was used at initialization.

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".
}
\value{
A list including:
\item{pvalue}{
A numeric value that represents the generated p-value due to G^2  test (see reference below).
}
\item{stat}{
A numeric value that represents the generated statistic due to G^2 test (see reference below).
}
\item{flag}{
A numeric value (control flag) which indicates whether the test was succesful (0) or not (1).
}
\item{stat_hash}{
The current hash object used for the statistics. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
\item{pvalue_hash}{
The current hash object used for the p-values. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
}
\references{
Tsamardinos, Ioannis, Laura E. Brown, and Constantin F. Aliferis. The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 2006 65(1): 31-78. 
}
\author{
R implementation and documentation: Giorgos Athineou <athineou@ics.forth.gr>
}
\note{
This test uses the functions gSquareBin and gSquareDis of the package pcalg in order to generate the pvalue for the G^2 test.
}
\seealso{
\code{\link{SES}, \link{testIndFisher}, \link{testIndLogistic}, \link{censIndLR}, \link{CondIndTests}}
}
\examples{
#simulate a dataset with binary data
dataset <- matrix(nrow = 50 , ncol = 101)
dataset <- apply(dataset, 2, function(i) sample(c(0,1),50, replace=TRUE))
#initialize binary target
target <- dataset[,101]
#remove target from the dataset
dataset <- dataset[,-101]

require(pcalg)

if(require("pcalg", quietly = TRUE))
{
  #run the gSquare conditional independence test for the binary class variable
  results <- gSquare(target, dataset, xIndex = 44, csIndex = c(10,20))
  results
  
  #require(gRbase) #for faster computations in the internal functions
  #run SES algorithm using the gSquare conditional independence test for the binary class variable
  sesObject <- SES(target , dataset , max_k=3 , threshold=0.05 , test="gSquare");
  #print summary of the SES output
  summary(sesObject);
  #plot the SES output
  plot(sesObject, mode="all");
}
}
\keyword{ G^2 }
\keyword{ Log Likelihood Ratio }
\keyword{ Conditional Independence Test }
