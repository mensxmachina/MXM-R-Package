\name{testIndFisher}
\alias{testIndFisher}
\title{
Fisher's conditional independence test for continous class variables.
}
\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS. 
}
\usage{
testIndFisher(target, dataset, xIndex, csIndex, dataInfo = NULL, univariateModels = NULL, 
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
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robustified version of Fisher's correlation coefficient via M-estimation available from the "rlm" in MASS. It takes more time than non robust version but it is suggested in case of outliers. Default value is FALSE. 
}
}
\details{
If hash = TRUE, testIndFisher requires the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.

Important: Use these arguments only with the same dataset that was used at initialization.

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".
}
\value{
A list including:
\item{pvalue}{
A numeric value that represents the generated p-value due to Fisher's method (see reference below).
}
\item{stat}{
A numeric value that represents the generated statistic due to Fisher's method (see reference below).
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
Peter Spirtes, Clark Glymour, and Richard Scheines. Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, second edition, January 2001.
}
\author{
Vincenzo Lagani and Ioannis Tsamardinos (Copyright 2012)

R implementation and documentation: Giorgos Athineou <athineou@ics.forth.gr> Vincenzo Lagani <vlagani@ics.forth.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{SES}, \link{testIndLogistic}, \link{gSquare}, \link{censIndLR}, \link{CondIndTests}}
}
\examples{

#simulate a dataset with continuous data
dataset <- matrix(runif(1000 * 200, 1, 1000), nrow = 1000 ) 
#the target feature is the last column of the dataset as a vector
target <- dataset[,200]
results <- testIndFisher(target, dataset, xIndex = 44, csIndex = 100)
#>results
#   $pvalue
#   [1] 0.5553586
#   
#   $stat
#   [1] 0.01869107
#   
#   $flag
#   [1] 1
#   
#   $stat_hash
#   NULL
#   
#   $pvalue_hash
#   NULL

#require(gRbase)  #for faster computations in the internal functions
#define class variable (here tha last column of the dataset)
target = 200;
#run the SES algorithm using the testIndFisher conditional independence test
sesObject <- SES(target , dataset , max_k=3 , threshold=0.05 , test="testIndFisher");
#print summary of the SES output
summary(sesObject);
#plot the SES output
plot(sesObject, mode="all");

}
\keyword{ Fisher's Test}
\keyword{ Conditional Independence Test }
