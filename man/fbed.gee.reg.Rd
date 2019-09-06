\name{Forward Backward Early Dropping selection regression with GEE}
\alias{fbed.gee.reg}
\title{
Forward Backward Early Dropping selection regression with GEE
}

\description{
Forward Backward Early Dropping selection regression with GEE.
}

\usage{
fbed.gee.reg(target, dataset, id, reps = NULL, ini = NULL, threshold = 0.05, 
wei = NULL, K = 0, test = "testIndGEEReg", correl = "exchangeable", se = "jack") 
}

\arguments{
\item{target}{
The class variable. This can be a numerical vector with continuous data, binary, discrete or an ordered factor variable. It can also be a factor variable with two levels only. Pay attention to this. If you have ordinal data, then supply an ordered factor variable.
}
\item{dataset}{
The set of candidate predictor variables provide a numerical matrix or a data.frame in case of categorical variables (columns = variables, rows = samples). In the case of ordinal regression, this can only be a numerical matrix, i.e. only continuous predictor variables.
}
\item{id}{
This is a numerical vector of the same size as target denoting the groups or the subjects. 
}
\item{reps}{
This is a numerical vector of the same size as target denoting the groups or the subjects. If you have longitudinal data and you know the time points, then supply them here. 
}
\item{ini}{
If you already have the test statistics and the logged p-values of the univariate associations (the first step of FBED) supply them as a 
list with the names "stat" and "pvalue" respectively.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for asmmmbsing p-values significance. Default value is 0.05.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
It is mentioned in the "geepack" that weights is not (yet) the weight as in sas proc genmod, and hence is not recommended to use.
}
\item{K}{
How many times should the process be repated? The default value is 0. 
}
\item{test}{
This is for the type of regression to be used, "testIndGEEReg", for Gaussian regression, 
"testIndGEEGamma" for Gamma regression, "testIndGEELogistic for logistic regression, "testIndGEENormLog" for linear regression witha log link, or "testIndGEEPois" for Poisson regression.
}
%\item{ncores}{
%How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
%}
\item{correl}{
The correlation structure. For the Gaussian, Logistic, Poisson and Gamma regression this can be either "exchangeable" (compound symmetry, suitable for clustered data) or "ar1" (AR(1) model, suitable for longitudinal data). For the ordinal logistic regression its only the "exchangeable" correlation sturcture. 
}
\item{se}{
The method for estimating standard errors. This is very important and crucial. The available options for Gaussian, Logistic, Poisson and Gamma regression are: a)  'san.se', the usual robust estimate. b) 'jack': if approximate jackknife variance estimate should be computed. 
c) 'j1s': if 1-step jackknife variance estimate should be computed and d) 'fij': logical indicating if fully iterated jackknife variance estimate should be computed. If you have many clusters (sets of repeated measurements) "san.se" is fine as it is astmpotically correct, plus jacknife estimates will take longer. If you have a few clusters, then maybe it's better to use jacknife estimates. 

The jackknife variance estimator was suggested by Paik (1988), which is quite suitable for cases when the number of subjects is small (K < 30), as in many biological studies. The simulation studies conducted by Ziegler et al. (2000) and Yan and Fine (2004) showed that the approximate jackknife estimates are in many cases in good agreement with the fully iterated ones. 
}
}

\details{
The algorithm is a variation of the usual forward selection. At every step, the most significant variable enters the selected variables set. In addition, only the significant variables stay and are further examined. The non signifcant ones are dropped. This goes until no variable can enter the set. The user has the option to redo this step 1 or more times (the argument K). In the end, a backward selection is performed to remove falsely selected variables. 
  
Since GEE are likelihood free, all significance tests take place using the Wald test, hence we decided not to have a backward phase. This 
algorithm is suitable for both clustered and longitudinal (glmm) data. 

If you specify a range of values of K it returns the results of fbed.reg for this range of values of K. For example, instead of runnning fbed.reg with K=0, K=1, K=2 and so on, you can run fbed.reg with K=0:2 and the selected variables found at K=2, K=1 and K=0 are returned. Note also, that you may specify maximum value of K to be 10, but the maximum value FBED used was 4 (for example). 

Make sure you load the library geepack first. 
}

\value{
If K is a single number a list including:
\item{res}{
A matrix with the selected variables, their test statistic and the associated p-value.
}
\item{info}{
A matrix with the number of variables and the number of tests performed (or models fitted) at each round (value of K).
}
\item{runtime}{
The runtime required.
}

If K is a sequence of numbers a list tincluding:
\item{univ}{
If you have used the log-likelihood ratio test this list contains the test statistics and the associated logged p-values of the univariate associations tests. If you have used the EBIc this list contains the eBIC of the univariate associations.
}
\item{res}{
The results of fbed.gee.reg with the maximum value of K. 
}
\item{mod}{
A list where each element refers to a K value. If you run FBEd with method = "LR" the selected variables, their test statistic and their p-value is returned. If you run it with method = "eBIC", the selected variables and the eBIC of the model with those variables are returned.
}
}

\references{
Borboudakis G. and Tsamardinos I. (2019). Forward-backward selection with early dropping. Journal of Machine Learning Research, 20(8): 1-39.
  
Liang  K.Y.  and  Zeger  S.L. (1986). Longitudinal data analysis using generalized linear models. Biometrika, 73(1): 13-22.

Prentice R.L. and Zhao L.P. (1991). Estimating equations for parameters in means and covariances
of multivariate discrete and continuous responses. Biometrics, 47(3): 825-839.

Paik M.C. (1988). Repeated measurement analysis for nonnormal data in small samples. Communications in
Statistics-Simulation and Computation, 17(4): 1155-1171.

Ziegler A., Kastner C., Brunner D. and Blettner M. (2000). Familial associations of lipid proles: A generalised
estimating equations approach. Statistics in medicine, 19(24): 3345-3357

Yan J. and Fine J. (2004). Estimating equations for association structures. Statistics in medicine, 23(6): 859-874

Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.
}


\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{ \link{fbed.glmm.reg}, \link{glmm.bsreg}, \link{MMPC.glmm}, \link{fbed.reg} }
}

\examples{
\dontrun{
require(lme4)
data(sleepstudy)
reaction <- sleepstudy$Reaction
days <- sleepstudy$Days
subject <- sleepstudy$Subject
x <- matrix(rnorm(180 * 20),ncol = 20) ## unrelated preidctor variables
m1 <- fbed.gee.reg(reaction, x, subject) 
m2 <- fbed.glmm.reg(reaction, x, subject, backward = FALSE) 
}
}

