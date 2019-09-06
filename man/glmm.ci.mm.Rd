\name{Symmetric conditional independence test with clustered data}
\alias{glmm.ci.mm}
\alias{gee.ci.mm}
\title{
Symmetric conditional independence test with clustered data
}

\description{
Symmetric conditional independence test with clustered data.
}

\usage{
glmm.ci.mm(ind1, ind2, cs = NULL, dat, group) 
gee.ci.mm(ind1, ind2, cs = NULL, dat, group, se = "jack") 
}

\arguments{
\item{ind1}{
The index of the one variable to be considered. 
}
\item{ind2}{
The index of the other variable to be considered. 
}
\item{cs}{
The index or indices of the conditioning set of variable(s). If you have no variables set this equal to 0.
}
\item{dat}{
A numerical matrix with data. 
}
\item{group}{
This is a numerical vector denoting the groups or the subjects.  
}
\item{se}{
The method for estimating standard errors. This is very important and crucial. The available options for Gaussian, Logistic and Poisson regression are: a) 'san.se': the usual robust estimate. b) 'jack': approximate jackknife variance estimate. c) 'j1s': if 1-step jackknife variance estimate and d) 'fij': fully iterated jackknife variance estimate. If you have many clusters (sets of repeated measurements) "san.se" is fine as it is asympotically correct, plus jacknife estimates will take longer. If you have a few clusters, then maybe it's better to use jacknife estimates. 

The jackknife variance estimator was suggested by Paik (1988), which is quite suitable for cases when the number of subjects is small (K < 30), as in many biological studies. The simulation studies conducted by Ziegler et al. (2000) and Yan and Fine (2004) showed that the approximate jackknife estimates are in many cases in good agreement with the fully iterated ones. 
}
}

\details{
Two linear random intercept models are fitted, one for each variable and the p-value of the hypothesis test that the other variable is significant 
is calculated. These two p-values are combined in a meta-analytic way. The models fitted are either linear, logistic and Poisson regression.
}

\value{
A vector including the test statistic, it's associated p-value and the relevant degrees of freedom. 
}

\references{
Tsagris M., Borboudakis G., Lagani V. and Tsamardinos I. (2018). Constraint-based Causal Discovery with Mixed Data. 
International Journal of Data Science and Analytics. 

https://link.springer.com/content/pdf/10.1007%2Fs41060-018-0097-y.pdf
 
http://nugget.unisa.edu.au/CD2017/papersonly/constraint-based-causal-r1.pdf

Paik M.C. (1988). Repeated measurement analysis for nonnormal data in small samples. Communications in Statistics-Simulation and Computation, 17(4): 1155-1171.

Ziegler A., Kastner C., Brunner D. and Blettner M. (2000). Familial associations of lipid profiles: A generalised estimating equations approach. Statistics in medicine, 19(24): 3345-3357

Yan J. and Fine J. (2004). Estimating equations for association structures. Statistics in medicine, 23(6): 859-874.

Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.
}

\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
  %%  ~~further notes~~
    %}
\seealso{
\code{ \link{pc.skel}, \link{condi}, \link{testIndGLMMReg} }
}

\examples{
## we generate two independent vectors of clustered data
s1 <- matrix(1.5, 4, 4)
diag(s1) <- 2.5
s2 <- matrix(1.5, 4, 4)
diag(s2) <- 2
x1 <- MASS::mvrnorm(10, rnorm(4), s1)  
x1 <- as.vector( t(x1) )
x2 <- MASS::mvrnorm(10, rnorm(4), s2)  
x2 <- as.vector( t(x2) )
id <- rep(1:10, each = 4)
glmm.ci.mm(1, 2, dat = cbind(x1,x2), group = id)
gee.ci.mm(1, 2, dat = cbind(x1,x2), group = id)
}
