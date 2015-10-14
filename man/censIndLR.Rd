\name{censIndLR}
\alias{censIndLR}
\title{
Conditional independence test for survival data
}
\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS. This test is based on the widely used Cox regression model (Cox, 1972).
}
\usage{
censIndLR(target, dataset, xIndex, csIndex, dataInfo = NULL, univariateModels = NULL, 
hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, robust = FALSE)
}
\arguments{
  \item{target}{
  A Survival object (class Surv from package survival) containing the time to event data (time) and the status indicator vector (event). View \link{Surv} documentation for more information.
}
  \item{dataset}{
  A numeric matrix or data frame, in case of categorical predictors (factors), containing the variables for performing the test. Rows as samples and columns as features.
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
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robustified version of Cox regression. Currently the robust version is not available for this test.
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
A numeric value that represents the generated p-value.
}
\item{stat}{
A numeric value that represents the generated statistic.
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
V. Lagani and I. Tsamardinos (2010). Structure-based variable selection for survival data. Bioinformatics Journal 16(15): 1887-1894.

Cox,D.R. (1972) Regression models and life-tables. J. R. Stat. Soc., 34, 187-220.
}
\author{
R implementation and documentation: Vincenzo Lagani <vlagani@ics.forth.gr>, Giorgos Athineou <athineou@ics.forth.gr>
}
\note{
This test uses the functions coxph and Surv of the package survival and the function anova (analysis of variance) of the package stats.
}
\seealso{
\code{\link{SES}, \link{testIndFisher}, \link{gSquare}, \link{testIndLogistic}, \link{Surv}, \link{coxph}, \link{anova}, \link{CondIndTests}}
}
\examples{
#create a survival simulated dataset
dataset <- matrix(nrow = 1000 , ncol = 100)
dataset <- apply(dataset, 1:2, function(i) runif(1, 1, 100))
dataset <- as.data.frame(dataset);
timeToEvent = rep(0,1000)
event = rep(0,1000)
c = rep(0,1000)
for(i in 1:1000)
{
  timeToEvent[i] = dataset[i,1] + 0.5*dataset[i,30] + 2*dataset[i,65] + runif(1, 0, 1);
  event[i] = sample(c(0,1),1)
  c[i] = runif(1, 0, timeToEvent[i]-0.5)
  if(event[i] == 0)
  {
    timeToEvent[i] = timeToEvent[i] - c[i]
  }
}

require(survival)

#init the Surv object class feature
if(require(survival, quietly = TRUE))
{
  target <- Surv(time=timeToEvent, event=event)
  
  #run the censIndLR   conditional independence test
  require(stats)
  res = censIndLR(target, dataset, xIndex=12, csIndex=c(35,7,4))
  res
  
  #run the SES algorithm using the censIndLR conditional independence
  #test for the survival class variable
  
  #require(gRbase) #for faster computations in the internal functions
  sesObject <- SES(target , dataset , max_k=1 , threshold=0.05 , test="censIndLR");
  #print summary of the SES output
  summary(sesObject);
  #plot the SES output
  plot(sesObject, mode="all");
}
}
\keyword{ Survival }
\keyword{ Log Likelihood Ratio }
\keyword{ Conditional Independence Test }
