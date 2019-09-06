\name{Fast MMPC for longitudinal and clustered data}
\alias{mmpc.glmm2}
\alias{mmpc.gee2}
\title{
mmpc.glmm2/mmpc.gee2: Fast Feature selection algorithm for identifying minimal feature subsets with correlated data
}

\description{
SES.glmm algorithm follows a forward-backward filter approach for feature selection in order to provide minimal, 
highly-predictive, statistically-equivalent, multiple feature subsets of a high dimensional dataset. 
See also Details. MMPC.glmm algorithm follows the same approach without generating multiple feature subsets. 
They are both adapted to longitudinal target variables.
}

\usage{
mmpc.glmm2(target, reps = NULL, group, dataset, max_k = 3, threshold = 0.05, 
test = NULL, ini = NULL, wei = NULL, slopes = FALSE, ncores = 1)

mmpc.gee2(target, reps = NULL, group, dataset, max_k = 3, threshold = 0.05, 
test = NULL, ini = NULL, wei = NULL, correl = "exchangeable", se = "jack", ncores = 1)
}

\arguments{
\item{target}{
The class variable. Provide a vector with continuous (normal), binary (binomial) or discrete (Poisson) data.
}
\item{reps}{
A numeric vector containing the time points of the subjects. It's length is equal to the length of the target variable. 
If you have clustered data, leave this NULL. 
}
\item{group}{
A numeric vector containing the subjects or groups. It must be of the same legnth as target. 
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables , rows = samples). Currently, only 
continuous datasets are supported. 
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). Integer, default value is 3.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
The conditional independence test to use. Default value is NULL. Currently, the only available conditional independence test 
are the \code{\link{testIndGLMMLogistic}}, \code{\link{testIndGLMMPois}}, \code{\link{testIndGLMMGamma}}, 
\code{\link{testIndGLMMNormLog}} and \code{\link{testIndGLMMReg}} which fit linear mixed models. 
}
\item{ini}{
This is a supposed to be a list. After running SES or MMPC with some hyper-parameters you might want to run SES again with 
different hyper-parameters. To avoid calculating the univariate associations (first step of SES and of MPPC) again, 
you can extract them from the first run of SES and plug them here. This can speed up the second run (and subequent runs of course) 
by 50\%. See the details and the argument "univ" in the output values. 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{slopes}{
Should random slopes for the ime effect be fitted as well? Default value is FALSE. 
}
\item{correl}{
The correlation structure. For the Gaussian, Logistic, Poisson and Gamma regression this can be either "exchangeable" 
(compound symmetry, suitable for clustered data) or "ar1" (AR(1) model, suitable for longitudinal data). 
}
\item{se}{
The method for estimating standard errors. This is very important and crucial. The available options for Gaussian, Logistic, 
Poisson and Gamma regression are: a) 'san.se': the usual robust estimate. b) 'jack': approximate jackknife variance estimate. 
c) 'j1s': if 1-step jackknife variance estimate and d) 'fij': fully iterated jackknife variance estimate. If you have many clusters 
(sets of repeated measurements) "san.se" is fine as it is asympotically correct, plus jacknife estimates will take longer. 
If you have a few clusters, then maybe it's better to use jacknife estimates. 

The jackknife variance estimator was suggested by Paik (1988), which is quite suitable for cases when the number of subjects 
is small (K < 30), as in many biological studies. The simulation studies conducted by Ziegler et al. (2000) and Yan and Fine (2004) 
showed that the approximate jackknife estimates are in many cases in good agreement with the fully iterated ones. 
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes 
and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will 
not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, 
where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in 
comparison to 1 core. Note also, that the amount of reduction is definetely not linear in the number of cores.
}
}

\details{
These are faster versions of MMPC using GLMM or GEE methodologies. See also \code{\link{mmpc2}} for more details on the algorithm.

If you want to use the GEE methodology, make sure you load the library geepack first. 
}

\value{
The output of the algorithm is an object of the class 'SES.glmm.output' for SES.glmm or 'MMPC.glmm.output' for MMPC.glmm including:
\item{selectedVars}{
The selected variables, i.e., the signature of the target variable.
}
\item{pvalues}{
For each feature included in the dataset, this vector reports the strength of its association with the target in the context of 
all other variables. Particularly, this vector reports the max p-values foudn when the association of each variable with the 
target is tested against different conditional sets. Lower values indicate higher association.
}
\item{univ}{
A vector with the logged p-values of the univariate associations. This vector is very important for subsequent runs of MMPC with 
different hyper-parameters. After running SES with some hyper-parameters you might want to run MMPCagain with different hyper-parameters. 
To avoid calculating the univariate associations (first step) again, you can take this list from the first run of SES and plug it in 
the argument "ini" in the next run(s) of MMPC. This can speed up the second run (and subequent runs of course) by 50\%. See the argument 
"univ" in the output values. 
}
\item{kapa_pval}{
A list with the same number of elements as the max$_k$. Every element in the list is a matrix. The first row is the logged p-values, 
the second row is the variable whose conditional association with the target variable was tests and the other rows are the 
conditioning variables. 
}
\item{max_k}{
The max_k option used in the current run.
}
\item{threshold}{
The threshold option used in the current run.
}
\item{n.tests}{
If you have set hash = TRUE, then the number of tests performed by SES or MMPC will be returned. If you have not set this to TRUE, 
the number of univariate associations will be returned. So be careful with this number.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the 
third element is the elapsed time.
}
\item{test}{
The character name of the statistic test used.
}
\item{slope}{
Whether random slopes for the time effects were used or not, TRUE or FALSE.
}
}

\references{
Tsagris, M., Lagani, V., & Tsamardinos, I. (2018). Feature selection for high-dimensional glmm data. BMC bioinformatics, 19(1), 17.

I. Tsamardinos, M. Tsagris and V. Lagani (2015). Feature selection for longitudinal data. Proceedings of the 10th conference 
of the Hellenic Society for Computational Biology & Bioinformatics (HSCBB15).

I. Tsamardinos, V. Lagani and D. Pappas (2012). Discovering multiple, equivalent biomarker signatures. In proceedings of the 7th 
conference of the Hellenic Society for Computational Biology & Bioinformatics - HSCBB12.

Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. 
Machine learning, 65(1), 31-78.

Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.

J. Pinheiro and D. Bates. Mixed-effects models in S and S-PLUS. Springer Science \& Business Media, 2006.

Liang K.Y. and Zeger S.L. (1986). Longitudinal data analysis using generalized linear models. Biometrika, 73(1): 13-22.

Prentice R.L. and Zhao L.P. (1991). Estimating equations for parameters in means and covariances of multivariate discrete and 
continuous responses. Biometrics, 47(3): 825-839.

Heagerty P.J. and Zeger S.L. (1996) Marginal regression models for clustered ordinal measurements. Journal of the American 
Statistical Association, 91(435): 1024-1036.

Paik M.C. (1988). Repeated measurement analysis for nonnormal data in small samples. Communications in Statistics-Simulation and 
Computation, 17(4): 1155-1171.

Ziegler A., Kastner C., Brunner D. and Blettner M. (2000). Familial associations of lipid profiles: A generalised estimating 
equations approach. Statistics in medicine, 19(24): 3345-3357

Yan J. and Fine J. (2004). Estimating equations for association structures. Statistics in medicine, 23(6): 859-874.
}

\author{
Ioannis Tsamardinos, Vincenzo Lagani

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Vincenzo Lagani <vlagani@csd.uoc.gr>  
}


\seealso{
\code{\link{CondIndTests}, \link{testIndGLMMReg} }
}

\examples{
\dontrun{
require(lme4)
data(sleepstudy)
reaction <- sleepstudy$Reaction
days <- sleepstudy$Days
subject <- sleepstudy$Subject
x <- matrix(rnorm(180 * 100), ncol = 100) ## unrelated predictor variables
m1 <- mmpc.glmm2(reaction, days, subject, x)
m2 <- MMPC.glmm(reaction, days, subject, x)
}
}

\keyword{ SES }
\keyword{ Mixed models }
\keyword{ Longitudinal data }
\keyword{ Variable Selection }
