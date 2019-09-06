\name{Backward phase of MMPC}
\alias{mmpcbackphase}
\title{
Backward phase of MMPC
}

\description{
Backward phase of MMPC.
}

\usage{
mmpcbackphase(target, dataset, max_k = 3, threshold = 0.05, test = NULL,
wei = NULL, R = 1) 
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor 
or a Surv object. See also Details.
}
\item{dataset}{
The data-set; provide either a data frame or a matrix (columns = variables , rows = samples).
Alternatively, provide an ExpressionSet (in which case rows are samples and columns are features, 
see bioconductor for details).
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). Integer, default value is 3.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
The conditional independence test to use. Type the test without " ", e.g. type testIndFisher, Not "testIndFisher". 
Default value is NULL. See also \code{\link{CondIndTests}}.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL.
}
%\item{hash}{
%A boolean variable which indicates whether (TRUE) or not (FALSE) to store the statistics calculated during SES 
%execution in a hash-type object. Default value is %FALSE. If TRUE a hashObject is produced.
%}
%\item{hashObject}{
%A List with the hash objects generated in a previous run of SES or MMPC. 
%Each time SES runs with "hash=TRUE" it produces a list of hashObjects that can be re-used in order to speed up 
%next runs of SES or MMPC. 
%Important: the generated hashObjects should be used only when the same dataset is re-analyzed, possibly with 
%different values of max_k and threshold.
\item{R}{
The number of permutations, set to 1 by default (no permutations based test). There is a trick to avoind doing 
all permutations. As soon as the number of times the permuted test statistic is more than the observed test 
statistic is more than 50 (if threshold = 0.05 and R = 999), the p-value has exceeded the signifiance level 
(threshold value) and hence the predictor variable is not significant. There is no need to continue do the extra 
permutations, as a decision has already been made. 
}
}

\details{
For each of the selected variables (dataset) the function performs conditional independence tests where the 
conditioning sets are formed from the other variables. All possible combinations are tried until the variable 
becomes non significant. The maximum size of the conditioning set is equal to max_k. This is called in the 
\code{\link{MMPC}} when the backward phase is requested.
}

\value{
A list including:
\item{met}{ 
A numerical vector of size equal to the number of columns of the dataset.
}
\item{counter}{
The number of tests performed.
}
\item{pvalues}{
The maximum p-value for the association of each predictor variable.
}
}

\references{
Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. 
Machine learning, 65(1), 31-78.
}

\author{
Ioannis Tsamardinos, Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>  
}

\seealso{
\code{\link{MMPC}, \link{mmhc.skel}, \link{CondIndTests}, \link{cv.mmpc}}
}

\examples{
set.seed(123)

#simulate a dataset with continuous data
dataset <- matrix(runif(500 * 100, 1, 100), ncol = 100)

#define a simulated class variable 
target <- 3 * dataset[, 10] + 2 * dataset[, 100] + 3 * dataset[, 20] + rnorm(500, 0, 5)

# MMPC algorithm 
m1 <- MMPC(target, dataset, max_k = 3, threshold = 0.05, test="testIndFisher");
m2 <- MMPC(target, dataset, max_k = 3, threshold = 0.05, test="testIndFisher", backward = TRUE);
x <- dataset[, m1@selectedVars]
mmpcbackphase(target, x, test = testIndFisher)
}

\keyword{ SES }
\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }
\keyword{ Variable Selection }
