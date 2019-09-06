\name{Fast MMPC}
\alias{mmpc2}
\title{
A fast version of MMPC
}

\description{
A fast version of MMPC
}

\usage{
mmpc2(target, dataset, max_k = 3, threshold = 0.05, test = "testIndLogistic", 
ini = NULL, wei = NULL, ncores = 1, backward = FALSE)
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
\item{dataset}{
The data-set; provide either a data frame or a matrix (columns = variables , rows = samples).
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). Integer, default value is 3.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
One of the following: "testIndBeta", "testIndReg", "testIndRQ", "testIndLogistic", "testIndMultinom", 
"testIndOrdinal", "testIndPois", "testIndQPois", "testIndZIP", "testIndNB", "testIndClogit", 
"testIndBinom", "testIndQBinom", "testIndIGreg", "censIndCR", "censIndWR", "censIndER", "testIndMMReg",
"testIndMVreg", "testIndMultinom", "testIndOrdinal", "testIndTobit", "testIndGamma", "testIndNormLog" or "testIndSPML". 
}
\item{ini}{
This is a supposed to be a list. To avoid calculating the univariate associations (first step of SES, MMPC and of FBED) again, you can extract them from the first run of omne of these algorithms and plug them here. 
This can speed up the second run (and subequent runs of course) by 50\%. See the details and the argument "univ" in the output values. 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}

%\item{R}{
%The number of permutations, set to 999 by default. There is a trick to avoind doing all permutations. As soon as the number of times the permuted test statistic is more than the observed test statistic is more than 50 (if threshold = 0.05 and R = 999), the p-value has exceeded the signifiance level (threshold value) and hence the predictor variable is not significant. There is no need to continue do the extra permutations, as a decision has already been made. 
%}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of 
variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). 
The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
\item{backward}{
If TRUE, the backward (or symmetry correction) phase will be implemented. This removes any falsely included variables in the parents and children set of the target variable. It calls the \code{link{mmpcbackphase}} for this purpose.
}
}

\details{
MMPC tests each feature for inclusion (selection). For each featurer it performsa conditional independence tets. Each test requires fitting two regression models, one without the feature and one with the feature included. In this version, we have changed the order of the tests. We find all possible subsets of the already selected features and for each of them we test each feature. This way, only half of the regression models the usual MMPC fits, are fitted. Also, less tests will be performed. It is the same algorithm, with a change in the sequence. 

We have seen a 50\% in the computational time, but the drawback is that if you want to run MMPC with different value of "max$_$k" and "alpha", this is not possible from here. This function is for oa signle pair of "max$_$k" and "alpha" values. It saves no test statistics, only p-values, no hashing and hence is memory efficient, but contains less information than \code{\link{MMPC}}.
}

\value{
The output of the algorithm is an S3 object including:
\item{selectedVars}{
The selected variables, i.e., the signature of the target variable.
}
\item{pvalues}{
For each feature included in the dataset, this vector reports the strength of its association with the target in the context of all other variable. Particularly, this vector reports the max p-values found when the association of each variable with the target is tested against different conditional sets. Lower values indicate higher association.
}
\item{univ}{
A vector with the logged p-values of the univariate associations. This vector is very important for subsequent runs of MMPC with different hyper-parameters. After running SES with some hyper-parameters you might want to run MMPCagain with different hyper-parameters. To avoid calculating the univariate associations (first step) again, you can take this list from the first run of SES and plug it in the argument "ini" in the next run(s) of MMPC. This can speed up the second run (and subequent runs of course) by 50\%. See the argument "univ" in the output values. 
}
\item{kapa_pval}{
A list with the same number of elements as the max$_k$. Every element in the list is a matrix. The first row is the logged p-values, the second row is the variable whose conditional association with the target variable was tests and the other rows are the conditioning variables. 
}
\item{max_k}{
The max_k option used in the current run.
}
\item{threshold}{
The threshold option used in the current run.
}
\item{n.tests}{
The number of tests performed by MMPC will be returned. 
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{test}{
The character name of the statistic test used.
}
}

\references{
Feature Selection with the R Package MXM: Discovering Statistically Equivalent Feature Subsets, Lagani, V. and Athineou, 
G. and Farcomeni, A. and Tsagris, M. and Tsamardinos, I. (2017). Journal of Statistical Software, 80(7).

Tsamardinos, I., Aliferis, C. F., & Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal 
relations. In Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). ACM.

Brown, L. E., Tsamardinos, I., & Aliferis, C. F. (2004). A novel algorithm for scalable and accurate Bayesian network learning. 
Medinfo, 711-715.

Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 
65(1), 31-78.

Tsagis M. (2018). Guide on performing feature selection with the R package MXM.
http://mensxmachina.org/site/wp-content/uploads/2018/04/Guide-on-performing-feature-selection-with-the-R-package-MXM.pdf
}

\author{
Ioannis Tsamardinos, Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>. 
}

\seealso{
\code{\link{MMPC}, \link{certificate.of.exclusion2}}
}

\examples{
set.seed(123)

#simulate a dataset with continuous data
dataset <- matrix(runif(100 * 100, 1, 100), ncol = 100)

#define a simulated class variable 
target <- 3 * dataset[, 10] + 2 * dataset[, 100] + 3 * dataset[, 20] + rnorm(100, 0, 5)
m <- median(target)
target[target < m] <- 0
target[abs(target) > 0 ] <- 1

m1 <- mmpc2(target, dataset, max_k = 3, threshold = 0.05, test="testIndLogistic")
m1$selectedVars  ## S3 class, $, not @
m1$runtime

m2 <- MMPC(target, dataset, max_k = 3, threshold = 0.05, test="testIndLogistic")
m2@selectedVars  ## S3 class, @, not $
m2@runtime
}

\keyword{ SES }
\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }
\keyword{ Variable Selection }
