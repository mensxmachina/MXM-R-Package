\name{Skeleton of the PC algorithm}
\alias{pc.skel}
\alias{pc.con}
\alias{pc.skel.boot}
\alias{glmm.pc.skel}
\alias{gee.pc.skel}
\title{
The skeleton of a Bayesian network produced by the PC algorithm
}

\description{
The skeleton of a Bayesian network produced by the PC algorithm. No orientations are involved. 
The pc.con is for continuous data only and it calls the same C++ as pc.skel. Hence, you are advised to use pc.skel. 
}

\usage{
pc.skel(dataset, method = "pearson", alpha = 0.01, rob = FALSE, R = 1, stat = NULL, 
ini.pvalue = NULL) 

pc.con(dataset, method = "pearson", alpha = 0.01) 

pc.skel.boot(dataset, method = "pearson", alpha = 0.01, R = 199, ncores = 1)

glmm.pc.skel(dataset, group, method = "comb.mm", alpha = 0.01, stat = NULL, 
ini.pvalue = NULL)

gee.pc.skel(dataset, group, se = "jack", method = "comb.mm", alpha = 0.01, stat = NULL, 
ini.pvalue = NULL)
}

\arguments{
\item{dataset}{
A matrix with the variables. The user must know if they are continuous or if they are categorical. 
\bold{If you have categorical data though, the user must transform the data.frame into a matrix. 
In addition, the numerical matrix must have values starting from 0. For example, 0, 1, 2, instead of "A", "B" and "C"}. 
In the case of mixed variables, continuous, binary and ordinal this must a data.frame and the non continuous variables 
must be ordered factors, even the binary variables. 

For the pc.con, pc.skel.boot and glmm.pc.skel, the dataset can only be a matrix.
}
\item{method}{
If you have continuous data, you can choose either "pearson", "spearman" or "distcor". The latter uses the distance correlation and 
should not be used with lots of observations as it is by default really slow. If you have categorical data, this must be "cat". 
If you have a mix of continuous, binary and ordinal data (we will expand the available dataset in the future) then choose 
"comb.fast" or "comb.mm". These two methods perform the symmetric test for mixed data (Tsagris et al., 2017). See details for 
more information on this. For the mixed models PC algorithm, this is the same argument, but cuurently only "comb.mm" is acecpted.
}
\item{group}{
This is to be used in the "glmm.pc.skel" and "gee.pc.skel" only. It is a vector for identifying the grouped data, the correlated 
observations, the subjects.
}
\item{se}{
The method for estimating standard errors. This is very important and crucial. The available options for Gaussian, Logistic, Poisson 
and Gamma regression are: a) 'san.se': the usual robust estimate. b) 'jack': approximate jackknife variance estimate. c) 'j1s': if 
1-step jackknife variance estimate and d) 'fij': fully iterated jackknife variance estimate. If you have many clusters (sets of 
repeated measurements) "san.se" is fine as it is asympotically correct, plus jacknife estimates will take longer. If you have a 
few clusters, then maybe it's better to use jacknife estimates. 

The jackknife variance estimator was suggested by Paik (1988), which is quite suitable for cases when the number of subjects is 
small (K < 30), as in many biological studies. The simulation studies conducted by Ziegler et al. (2000) and Yan and Fine (2004) 
showed that the approximate jackknife estimates are in many cases in good agreement with the fully iterated ones. 
}
\item{alpha}{
The significance level ( suitable values in (0, 1) ) for assessing the p-values. Default value is 0.01.
}
\item{rob}{
This is for robust estimation of the Pearson correlation coefficient. Default value is FALSE.
}

\item{R}{
The number of permutations to be conducted. This is taken into consideration for the "pc.skel" only. The Pearson correlation coefficient 
is calculated and the p-value is assessed via permutations. There was an extra argument for plotting the skeleton but it does not work 
with the current visualisation packages, hence we removed the argument. Use \code{\link{plotnetwork}} to plot the skeleton.  

In the pc.skel.boot this is the number of bootstrap resamples to draw. The PC algorithm is performed in each bootstrap sample. In the end, 
the adjacency matrix on the observed data is returned, along with another adjacency matrix produced by the bootstrap. The latter one contains 
values from 0 to 1 indicating the proportion of times an edge between two nodes was present.
}
\item{ncores}{
The number of cores to use. By default this is set to 1.
}
\item{stat}{
If you have the initial test statistics (univariate associations) values supply them here.
}
\item{ini.pvalue}{
If you have the initial p-values of the univariate associations supply them here.
}
}

\details{
The PC algorithm as proposed by Spirtes et al. (2000) is implemented. The variables must be either continuous or categorical, only. 
The skeleton of the PC algorithm is order independent, since we are using the third heuristic (Spirte et al., 2000, pg. 90). 
At every ste of the alogirithm use the pairs which are least statistically associated. The conditioning set consists of variables 
which are most statistically associated with each either of the pair of variables. 

For example, for the pair (X, Y) there can be two coniditoning sets for example (Z1, Z2) and (W1, W2). All p-values and test statistics 
and degrees of freedom have been computed at the first step of the algorithm. Take the p-values between (Z1, Z2) and (X, Y) and between 
(Z1, Z2) and (X, Y). The conditioning set with the minimum p-value is used first. If the minimum p-values are the same, use the second 
lowest p-value. In the event of 2 or more p-values being the same (with permutations for example), the test statistic divided by the 
degrees of freedom is used as a means of choosing which conditioning set is to be used first. If two or more p-values are below the 
machine epsilon (.Machine$double.eps which is equal to 2.220446e-16), all of them are set to 0. To make the comparison or the ordering 
feasible we use the logarithm of the p-value. Hence, the logarithm of the p-values is always calculated and used.

In the case of the \eqn{G^2} test of independence (for categorical data) we have incorporated a rule of thumb. I the number of samples 
is at least 5 times the number of the parameters to be estimated, the test is performed, otherwise, independence is not rejected 
(see Tsamardinos et al., 2006).

The "comb.fast" and "comb.mm" methods are used with mixed variables, continuous, binary and ordinal. The "comb.mm" performs two 
log-likelihood ratio tests. For every pair of variables each of the two variables is treated as response and the suitable regression model 
is fitted. Then, two likelihood ratio tests are performed and the 2 p-values are combined in a meta-analytic way. In the case of 
"comb.fast" one regression model is fitted, the easiest (between the two) to implement. The ordering of the "easiness" is as follows: 
linear regression > logistic regression > ordinal regression.  

The "pc.con" is a faster implementation of the PC algorithm but for continuous data only, without the robust option, unlike pc.skel 
which is more general and even for the continuous datasets slower. pc.con accepts only "pearson" and "spearman" as correlations. 

If there are missing values they are placed by their median in case of continuous data and by their mode (most frequent value) if they 
are categorical. 

The "glmm.pc.skel" and "gee.pc.skel" are designed for clustered or grouped data. It uses linear mixed models and works in the same way 
as the PC with mixed data. For each variable, a random intercepts model is fitted and the significance of the other variable is assessed. 
The two p-values are meta-analytically combined. 

For all cases, we return the maximum logged p-value of the conditional independence tests between the pairs. This can be used to order 
the strength of association between pairs of variables. In addition, one can use it to estimate the AUC. See the example in \code{\link{bn.skel.utils}}.

If you want to use the GEE methodology, make sure you load the library geepack first. 
}

\value{
A list including:
\item{stat}{
The test statistics of the univariate associations. 
}
\item{ini.pvalue}{
The initival univariate associations p-values.
}
\item{pvalue}{
The logarithm of the maximum p-values from every conditional independence test. Note also, that if there is a log p-value smaller that 
the log of the threshold value that means there is an edge. This can be used to estimate the FDR (Tsamardinos and Brown, 2008). See 
details for more information. This is also an estimate of the strength of the association between pairs of variables. 

At the moment this is not possible for method = c("pearson", "spearman", "cat") with R = 1, or R > 1 because these cases are handled in 
C++. We will add those cases in the future. 
}
\item{info}{
Some summary statistics about the edges, minimum, maximum, mean, median number of edges.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the 
third element is the elapsed time.
}
\item{kappa}{
The maximum value of k, the maximum cardinality of the conditioning set at which the algorithm stopped.
}
\item{density}{
The number of edges divided by the total possible number of edges, that is #edges / \eqn{n(n-1)/2}, where \eqn{n} is the number of variables.
}
\item{info}{
Some summary statistics about the edges, minimum, maximum, mean, median number of edges.
}
\item{G}{
The adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that i and j have an edge between them.
}
\item{sepset}{
A list with the separating sets for every value of k.
}
\item{title}{
The name of the dataset. 
}

Bear in mind that the values can be extracted with the $ symbol, i.e. this is an S3 class output. 

}
\references{
Spirtes P.,  Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, 3nd edition.

Tsamardinos, Ioannis, and Laura E. Brown. Bounding the False Discovery Rate in Local Bayesian Network Learning. AAAI, 2008.

Tsagris M., Borboudakis G., Lagani V. and Tsamardinos I. (2018). Constraint-based Causal Discovery with Mixed Data. 
International Journal of Data Science and Analytics, 6: 19-30. 

http://nugget.unisa.edu.au/CD2017/papersonly/constraint-based-causal-r1.pdf

Tsagris M. Bayesian network learning with the PC algorithm: an improved and correct variation. 
Applied Artificial Intelligence (Accepted for publication).

Sedgewick, A. J., Ramsey, J. D., Spirtes, P., Glymour, C., & Benos, P. V. (2017). Mixed Graphical Models for Causal Analysis of 
Multi-modal Variables. arXiv preprint arXiv:1704.02621.

Szekely G.J. and Rizzo, M.L. (2014). Partial distance correlation with methods for dissimilarities. The Annals of Statistics, 42(6): 
2382--2412.

Szekely G.J. and Rizzo M.L. (2013). Energy statistics: A class of statistics based on distances. Journal of Statistical Planning 
and Inference 143(8): 1249--1272.

Paik M.C. (1988). Repeated measurement analysis for nonnormal data in small samples. Communications in Statistics-Simulation and 
Computation, 17(4): 1155-1171.

Ziegler A., Kastner C., Brunner D. and Blettner M. (2000). Familial associations of lipid profiles: A generalised estimating 
equations approach. Statistics in medicine, 19(24): 3345-3357

Yan J. and Fine J. (2004). Estimating equations for association structures. Statistics in medicine, 23(6): 859-874.

Liang K.Y. and Zeger S.L. (1986). Longitudinal data analysis using generalized linear models. Biometrika, 73(1): 13-22.

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
\code{ \link{bn.skel.utils}, \link{mmhc.skel}, \link{corfs.network}, \link{local.mmhc.skel}  }
}

\examples{
# simulate a dataset with continuous data
y <- rdag2(500, p = 20, nei = 3)
ind <- sample(1:20, 20)
x <- y$x[, ind]
a <- mmhc.skel(x, max_k = 3, threshold = 0.05, test = "testIndFisher" ) 
b <- pc.skel( x, method = "pearson", alpha = 0.05 ) 
a$runtime  
b$runtime  

}

\keyword{ pc algorithm }