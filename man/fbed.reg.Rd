\name{Forward Backward Early Dropping selection regression}
\alias{fbed.reg}
\title{
Forward Backward Early Dropping selection regression
}

\description{
Forward Backward Early Dropping selection regression.
}

\usage{
fbed.reg(target, dataset, ini = NULL, test = NULL, threshold = 0.05, wei = NULL, 
K = 0, method = "LR", gam = NULL, backward = TRUE) 
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. 
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). 
}
\item{ini}{
If you already have the test statistics and the p-values of the univariate associations (the first step of FBED) supply them as a list with the names "stat" and "pvalue" respectively. If you have used the EBIc this list contains the eBIC of the univariate associations. Note, that the "gam" argument must be the same though.
}
\item{test}{
The available tests: "testIndReg", "testIndPois", "testIndNB", "testIndLogistic", 
"testIndMMReg", "testIndBinom", "censIndCR", "censIndWR", "testIndBeta", "testIndZIP", 
"testIndGamma, "testIndNormLog", "testIndTobit", "testIndQPois",  "testIndQBinom", 
"testIndFisher" "testIndMultinom", "testIndOrdinal", "testIndClogit" and "testIndSPML". Note that not all 
of them work with eBIC. The only test that does not have "eBIC" and no weights (input argument "wei") 
is the test "gSquare". The tests "testIndFisher" and "testIndSPML" have no weights.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when robust is set 
to TRUE. If you want to use the "testIndBinom", then supply the successes in the y and the trials here.
}
\item{K}{
How many times should the process be repeated? The default value is 0. You can also specify a range of values of K, say 0:4 for example.
}
\item{method}{
Do you want the likelihood ratio test to be performed ("LR" is the default value) or perform the selection using the "eBIC" criterion (BIC is a special case). 
}
\item{gam}{
In case the method is chosen to be "eBIC" one can also specify the \eqn{gamma} parameter. The default value is "NULL", so that the value is automatically calculated.
}
\item{backward}{
After the Forward Early Dropping phase, the algorithm proceeds witha the usual Backward Selection phase. The default value is set to TRUE. It is advised to perform this step as maybe some variables are false positives, they were wrongly selected. 

The backward phase using likelihood ratio test and eBIc are two different functions and can be called directly by the user. SO, if you want for example to perform a backard regression with a different threshold value, just use these two functions separately. 
}
%\item{ncores}{
%How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
%}
}

\details{
The algorithm is a variation of the usual forward selection. At every step, the most significant variable enters the selected variables set. In addition, only the significant variables stay and are further examined. The non signifcant ones are dropped. This goes until no variable can enter the set. The user has the option to re-do this step 1 or more times (the argument K). In the end, a backward selection is performed to remove falsely selected variables. Note that you may have specified, for example, K=10, but the maximum value FBED used can be 4 for example. 

The "testIndQPois" and "testIndQBinom" do not work for method = "eBIC" as there is no BIC associated with quasi Poisson and quasi binomial
regression models.

If you specify a range of values of K it returns the results of fbed.reg for this range of values of K. For example, instead of runnning fbed.reg with K=0, K=1, K=2 and so on, you can run fbed.reg with K=0:2 and the selected variables found at K=2, K=1 and K=0 are returned. Note also, that you may specify maximum value of K to be 10, but the maximum value FBED used was 4 (for example). 
}

\value{
If K is a single number a list including:
\item{univ}{
If you have used the log-likelihood ratio test this list contains the test statistics and the associated p-values of the univariate associations tests. If you have used the EBIc this list contains the eBIC of the univariate associations. Note, that the "gam" argument must be the same though.
}
\item{res}{
A matrix with the selected variables, their test statistic and the associated logged p-value.
}
\item{info}{
A matrix with the number of variables and the number of tests performed (or models fitted) at each round (value of K). This refers to the
forward phase only.
}
\item{runtime}{
The runtime required.
}
\item{back.rem}{
The variables removed in the backward phase.
}
\item{back.n.tests}{
The number of models fitted in the backward phase.
}

If K is a sequence of numbers a list tincluding:
\item{univ}{
If you have used the log-likelihood ratio test this list contains the test statistics and the associated p-values of the univariate associations tests. If you have used the EBIc this list contains the eBIC of the univariate associations.
}
\item{res}{
The results of fbed.reg with the maximum value of K. 
}
\item{mod}{
A list where each element refers to a K value. If you run FBEd with method = "LR" the selected variables, their test statistic and their p-value is returned. If you run it with method = "eBIC", the selected variables and the eBIC of the model with those variables are returned.
}
}

\references{
Borboudakis G. and Tsamardinos I. (2019). Forward-backward selection with early dropping. Journal of Machine Learning Research, 20(8): 1-39.

Tsagis M. (2018). Guide on performing feature selection with the R package MXM.
http://mensxmachina.org/site/wp-content/uploads/2018/04/Guide-on-performing-feature-selection-with-the-R-package-MXM.pdf
}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{ \link{fs.reg}, \link{ebic.bsreg}, \link{bic.fsreg}, \link{MMPC} }
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix( runif(100 * 50, 1, 100), ncol = 50 )
#define a simulated class variable 
target <- rt(100, 10)

a1 <- fbed.reg(target, dataset, test = "testIndReg") 
y <- rpois(100, 10)
a2 <- fbed.reg(y, dataset, test = "testIndPois") 
a3 <- MMPC(y, dataset)
}

