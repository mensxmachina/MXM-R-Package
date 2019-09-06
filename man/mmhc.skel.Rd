\name{Skeleton of the max-min hill-climbing (MMHC) algorithm}
\alias{mmhc.skel}
\alias{glmm.mmhc.skel}
\alias{gee.mmhc.skel}

\title{
The skeleton of a Bayesian network as produced by MMHC
}

\description{
The skeleton of a Bayesian network produced by MMHC. No orientations are involved.
}

\usage{
mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher", type = "MMPC",
hash = FALSE, backward = TRUE, symmetry = TRUE, nc = 1, ini.pvalue = NULL) 

glmm.mmhc.skel(dataset, group, max_k = 3, threshold = 0.05, test = "testIndLMM",
type = "MMPC.glmm", hash = FALSE, symmetry = TRUE, nc = 1, ini.pvalue = NULL)

gee.mmhc.skel(dataset, group, max_k = 3, threshold = 0.05, test = "testIndGEEReg",
type = "MMPC.gee", se = "jack", hash = FALSE, symmetry = TRUE, nc = 1, 
ini.pvalue = NULL)
}

\arguments{
\item{dataset}{
A matrix with the variables. The user must know if they are continuous or if they are categorical. If you have a matrix with categorical data, i.e. 0, 1, 2, 3 where each number indicates a category, the minimum number for each variable must be 0. For the "glmm.mmhc.skel" this must be only a matrix.
}
\item{group}{
This is to be used in the "glmm.mmhc.skel" and "gee.mmhc.skel" only. It is a vector for identifying the grouped data, the correlated observations, the subjects.
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details of SES or MMPC). 
}
\item{threshold}{
Threshold ( suitable values in (0, 1) ) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
The conditional independence test to use. Default value is "testIndFisher". This procedure allows for "testIndFisher", "testIndSPearman" for continuous variables and "gSquare" for categorical variables. In case the dataset is a data.frame with mixed types of data leave this NULL and an appropriate test will be selected. See \code{\link{MMPC}} for the automatic choice of tests.

For the "glmm.mmhc.skel" the available tests in \code{\link{MMPC.glmm}}.
}
\item{type}{
The type of variable selection to take place for each variable (or node). The default (and standard) is "MMPC" and "MMPC.glmm". You can also choose to run it via "SES" and "SES.glmm" and thus allow for multiple signatures of variables to be connected to a variable.
}
\item{se}{
The method for estimating standard errors. This is very important and crucial. The available options for Gaussian, Logistic, Poisson and Gamma regression are: a) 'san.se': the usual robust estimate. b) 'jack': approximate jackknife variance estimate. c) 'j1s': if 1-step jackknife variance estimate and d) 'fij': fully iterated jackknife variance estimate. If you have many clusters (sets of repeated measurements) "san.se" is fine as it is asympotically correct, plus jacknife estimates will take longer. If you have a few clusters, then maybe it's better to use jacknife estimates. 

The jackknife variance estimator was suggested by Paik (1988), which is quite suitable for cases when the number of subjects is small (K < 30), as in many biological studies. The simulation studies conducted by Ziegler et al. (2000) and Yan and Fine (2004) showed that the approximate jackknife estimates are in many cases in good agreement with the fully iterated ones. 
}
\item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE). Default value is FALSE. If TRUE a hashObject is produced and hence more memory is required. If you want to know the number of tests executed for each variable then make it TRUE.
}
\item{backward}{
If TRUE, the backward (or symmetry correction) phase will be implemented. This removes any falsely included variables in the parents and children set of the target variable. It call the \code{\link{mmpcbackphase}} for this purpose. For perm.ses and wald.ses this is not yet applicable.
}
\item{symmetry}{
In order for an edge to be added, a statistical relationship must have been found from both directions. If you want this symmetry correction to take place, leave this boolean variable to TRUE. If you set it to FALSE, then if 
a relationship between Y and X is detected but not between X and Y, the edge is still added. 
}
\item{nc}{
How many cores to use. This plays an important role if you have many variables, say thousands or so. You can try with nc = 1 and with nc = 4 for example to see the differences. If you have a multicore machine, this is a must 
option. There was an extra argument for plotting the skeleton but it does not work with the current visualisation 
packages, hence we removed the argument. Use \code{\link{plotnetwork}} to plot the skeleton.   
}
\item{ini.pvalue}{
This is a list with the matrix of the univariate p-values. If you want to run mmhc.skel again, the univariate associations need not be calculated again.
}
}

\details{
The MMPC is run on every variable. The backward phase (see Tsamardinos et al., 2006) takes place automatically. After all variables have been used, the matrix is checked for inconsistencies and they are corrected. 

A trick mentioned in that paper to make the procedure faster is the following. In the k-th variable, the algorithm checks how many previously scanned variables have an edge with the this variable and keeps them (it discards the other variables with no edge) along with the next (unscanned) variables. 

This trick reduces time, but can lead to different results. For example, if the i-th variable is removed, the k-th node might not remove an edge between the j-th variable, simply because the i-th variable that could d-sepate them is missing. 

The user is given this option via the argument "fast", which can be either TRUE or FALSE. Parallel computation is also available. 

}

\value{
A list including:
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{density}{
The number of edges divided by the total possible number of edges, that is #edges / \eqn{n(n-1)/2}, where \eqn{n} is the number of variables.
}
\item{info}{
Some summary statistics about the edges, minimum, maximum, mean, median number of edges.
}
\item{ms}{
If you run "MMPC" for each variable this is NULL. If you run "SES" is a vector denoting which variables had more than one signature, i.e. more than one set of variables associated with them.
}
\item{ntests}{
The number of tests MMPC (or SES) performed at each variable.
}
\item{ini.pvalue}{
A matrix with the p-values of all pairwise univariate assocations.
}
\item{pvalue}{
A matrix with the final p-values. These are the maximum p-values calculated during the process. When the process finishes, the matrix is not symmetric. 
It becomes symmetric though by keeping the maximum between any two off-diagonal elements. These p-values now can be used to sort the strength of the edges. 
If you know the true adjacency matrix you can use them and create a ROC curve (see \code{\link{bn.skel.utils}} for more information).
}
\item{G}{
The adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that i and j have an edge between them.
}

Bear in mind that the values can be extracted with the $ symbol, i.e. this is an S3 class output. 
}

\references{
Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.

Brown L. E., Tsamardinos I., and Aliferis C. F. (2004). A novel algorithm for scalable and accurate Bayesian network learning. Medinfo, 711-715.

Tsamardinos, Ioannis, and Laura E. Brown. Bounding the False Discovery Rate in Local Bayesian Network Learning. AAAI, 2008.

Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.

Paik M.C. (1988). Repeated measurement analysis for nonnormal data in small samples. Communications in Statistics-Simulation and Computation, 17(4): 1155-1171.

Ziegler A., Kastner C., Brunner D. and Blettner M. (2000). Familial associations of lipid profiles: A generalised estimating equations approach. Statistics in medicine, 19(24): 3345-3357

Yan J. and Fine J. (2004). Estimating equations for association structures. Statistics in medicine, 23(6): 859-874.

Liang K.Y. and Zeger S.L. (1986). Longitudinal data analysis using generalized linear models. Biometrika, 73(1): 13-22.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{pc.skel}, \link{pc.or}, \link{corfs.network}, \link{bn.skel.utils} }
}

\examples{
# simulate a dataset with continuous data
y <- rdag2(500, p = 20, nei = 3)
x <- y$x
a <- mmhc.skel(x, max_k = 5, threshold = 0.01, test = "testIndFisher" ) 
b <- pc.skel( x, alpha = 0.01 ) 
a$runtime  
b$runtime 
}

\keyword{ MMHC algorithm }