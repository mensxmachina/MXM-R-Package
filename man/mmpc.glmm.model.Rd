\name{Generalised linear mixed models based on glmm SES and MMPC outputs}
\alias{mmpc.glmm.model}
\alias{ses.glmm.model}
\alias{mmpc.gee.model}
\alias{ses.gee.model}
\title{
Generalised linear mixed model(s) based obtained from glmm SES or MMPC
}

\description{
One or more regression models obtained from SES or MMPC, are returned. 
}

\usage{
mmpc.glmm.model(target, dataset, reps = NULL, group, slopes = FALSE, wei = NULL,
mmpcglmm.Object, test = NULL) 

ses.glmm.model(target, dataset, reps = NULL, group, slopes = FALSE, wei = NULL,
sesglmm.Object, nsignat = 1, test = NULL) 

mmpc.gee.model(target, dataset, reps = NULL, group, correl = "exchangeable", 
se = "jack", wei = NULL, mmpcgee.Object, test = NULL)

ses.gee.model(target, dataset, reps = NULL, group, correl = "exchangeable", 
se = "jack", wei = NULL, sesgee.Object, nsignat = 1, test = NULL)
}

\arguments{
\item{target}{
The class variable. Provide a vector with continuous (normal), binary (binomial) or discrete (Poisson) data.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = observations). Currently, only continuous 
datasets are supported. 
}
\item{reps}{
A numeric vector containing the time points of the subjects. Its length is equal to the length of the target variable. 
If you have clustered data, leave this NULL. 
}
\item{group}{
A numeric vector containing the subjects or groups. It must be of the same legnth as target. 
}
\item{slopes}{
Should random slopes for the ime effect be fitted as well? Default value is FALSE. 
}
\item{correl}{
The correlation structure. For the Gaussian, Logistic, Poisson and Gamma regression this can be either "exchangeable" (compound symmetry, 
suitable for clustered data) or "ar1" (AR(1) model, suitable for longitudinal data). For the ordinal logistic regression its only the 
"exchangeable" correlation sturcture.  
}
\item{se}{
The method for estimating standard errors. This is very important and crucial. The available options for Gaussian, Logistic, Poisson and 
Gamma regression are: a)  'san.se', the usual robust estimate. b) 'jack': if approximate jackknife variance estimate should be computed. 
c) 'j1s': if 1-step jackknife variance estimate should be computed and d) 'fij': logical indicating if fully iterated jackknife variance 
estimate should be computed. If you have many clusters (sets of repeated measurements) "san.se" is fine as it is astmpotically correct, 
plus jacknife estimates will take longer. If you have a few clusters, then maybe it's better to use jacknife estimates. 

The jackknife variance estimator was suggested by Paik (1988), which is quite suitable for cases when the number of subjects is small 
(K < 30), as in many biological studies. The simulation studies conducted by Ziegler et al. (2000) and Yan and Fine (2004) showed that 
the approximate jackknife estimates are in many cases in good agreement with the fully iterated ones. 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{mmpcglmm.Object}{
An object with the results of an MMPC.glmm run.
}
\item{sesglmm.Object}{
An object with the results of a SES.glmm run.
}
\item{mmpcgee.Object}{
An object with the results of an MMPC.gee run.
}
\item{sesgee.Object}{
An object with the results of an SES.gee run.
}
\item{nsignat}{
How many signatures to use. If nsignat = 1 (default value) the first set of variables will be used for the model. If you want more, 
then specify the nubmer of signatures you want. If you want the models based on all signatures, specify "all". If you put a number 
which is higher than the number of signatures, all models will be returned.
}
\item{test}{
The conditional independence test to use. Default value is NULL. Currently, the only available conditional independence test are 
the \code{\link{testIndGLMMReg}}, \code{\link{testIndGLMMLogistic}}, \code{\link{testIndGLMMPois}} and \code{\link{testIndLMM}} which fit 
linear mixed models. 
}
}

\details{
 This command is useful if you want to see all models and check for example their fitting ability, MSE in linear models for exmaple.
}

\value{
A list including:
\item{mod}{
Depending on the number of signatures requested, one or models will be returned.
}
\item{signature}{
A matrix (or just one vector if one signature only) with the variables of each signature, along with the BIC of the corresponding regression model.
}
}
\references{
Tsagris, M., Lagani, V., & Tsamardinos, I. (2018). Feature selection for high-dimensional glmm data. BMC bioinformatics, 19(1), 17.

I. Tsamardinos, M. Tsagris and V. Lagani (2015). Feature selection for longitudinal data. Proceedings of the 10th conference of 
the Hellenic Society for Computational Biology & Bioinformatics (HSCBB15)

Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.

Pinheiro J. and D. Bates. Mixed-effects models in S and S-PLUS. Springer Science \& Business Media, 2006.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{SES}, \link{MMPC}, \link{cv.ses}, \link{cv.mmpc} }
}

\examples{
\dontrun{
require(lme4)
data(sleepstudy)
reaction <- sleepstudy$Reaction
days <- sleepstudy$Days
subject <- sleepstudy$Subject
x <- matrix(rnorm(180 * 50),ncol = 50) ## unrelated predictor variables
m1 <- SES.glmm(reaction, days, subject, x)
m2 <- MMPC.glmm(reaction, days, subject, x)
mod <- mmpc.glmm.model(reaction, dataset = x, reps = days, group = subject, slopes = FALSE, 
mmpcglmm.Object = m2) 
}
}



\keyword{ Regression modelling }