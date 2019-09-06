\name{Forward Backward Early Dropping selection regression with GLMM}
\alias{fbed.glmm.reg}
\title{
Forward Backward Early Dropping selection regression with GLMM
}

\description{
Forward Backward Early Dropping selection regression with GLMM.
}

\usage{
fbed.glmm.reg(target, dataset, id, reps = NULL, ini = NULL, threshold = 0.05, 
wei = NULL, K = 0, method = "LR", gam = NULL, backward = TRUE, test = "testIndGLMMReg")
}

\arguments{
\item{target}{
The class variable. This can be a numerical vector with continuous data, binary or discrete valued data. It can also be a 
factor variable with two levels only.
}
\item{dataset}{
The dataset; provide a numerical a matrix (columns = variables, rows = samples).
}
\item{id}{
This is a numerical vector of the same size as target denoting the groups or the subjects. 
}
\item{reps}{
This is a numerical vector of the same size as target denoting the groups or the subjects. If you have longitudinal data and you know the time points, then supply them here. 
}
\item{ini}{
If you already have the test statistics and the p-values of the univariate associations (the first step of FBED) supply them as a 
list with the names "stat" and "pvalue" respectively in the case of method = "LR". In the case of method = "eBIC" this is a list with an 
element called "bic".
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for asmmmbsing p-values significance. Default value is 0.05.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{K}{
How many times should the process be repated? The default value is 0. 
}
\item{method}{
Do you want the likelihood ratio test to be performed ("LR" is the default value) or perform the selection using the "eBIC" criterion (BIC is a special case). 
}
\item{gam}{
In case the method is chosen to be "eBIC" one can also specify the \eqn{gamma} parameter. The default value is "NULL", so that the value is 
automatically calculated.
}
\item{backward}{
After the Forward Early Dropping phase, the algorithm proceeds witha the usual Backward Selection phase. The default value is set to TRUE. It is advised to perform this step as maybe some variables are false positives, they were wrongly selected. 

The backward phase using likelihood ratio test and eBIc are two different functions and can be called directly by the user. SO, if you want for example to perform a backard regression with a different threshold value, just use these two functions separately. 
}
\item{test}{
This is for the type of regression to be used, "testIndGLMMReg", for Gaussian regression, "testIndGLMMLogistic for logistic regression or "testIndGLMMPois" for Poisson regression, "testIndGLMMGamma" for Gamma regression and "testIndGLMMNormLog" for Gaussian regression with log link and "testIndGLMMOrdinal" for ordinal regression.
}
%\item{ncores}{
%How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
    %}
}

\details{
The algorithm is a variation of the usual forward selection. At every step, the most significant variable enters the selected variables set. In addition, only the significant variables stay and are further examined. The non signifcant ones are dropped. This goes until no variable can enter the set. The user has the option to redo this step 1 or more times (the argument K). In the end, a backward selection is performed to remove falsely selected variables. 

Bear in mind that for the "gaussian" case, the forward phase takes place using the F test (Wald statistic calibrated against an F distribution). The backward phase though takes place using the log-likelihood ratio test.

If you specify a range of values of K it returns the results of fbed.reg for this range of values of K. For example, instead of runnning fbed.reg with K=0, K=1, K=2 and so on, you can run fbed.reg with K=0:2 and the selected variables found at K=2, K=1 and K=0 are returned. Note also, that you may specify maximum value of K to be 10, but the maximum value FBED used was 4 (for example). 
}

\value{
If K is a single number a list including:
\item{res}{
A matrix with the selected variables, their test statistic and the associated loged p-value.
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
The results of fbed.glmm.reg with the maximum value of K. 
}
\item{mod}{
A list where each element refers to a K value. If you run FBED with method = "LR" the selected variables, their test statistic and their p-value is returned. If you run it with method = "eBIC", the selected variables and the eBIC of the model with those variables are returned.
}
}

\references{
Borboudakis G. and Tsamardinos I. (2019). Forward-backward selection with early dropping. Journal of Machine Learning Research, 20(8): 1-39.
  
Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.
}


\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{ \link{fbed.gee.reg}, \link{glmm.bsreg}, \link{MMPC.glmm}, \link{fbed.reg} }
}

\examples{
\dontrun{
require(lme4)
data(sleepstudy)
reaction <- sleepstudy$Reaction
subject <- sleepstudy$Subject
x <- matrix(rnorm(180 * 20),ncol = 20) ## unrelated predictor variables
m1 <- fbed.glmm.reg(reaction, x, subject) 
m2 <- MMPC.glmm(target = reaction, group = subject, dataset = x)
}
}

