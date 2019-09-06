\name{Bayesian Network construction using a hybrid of MMPC and PC}
\alias{mmpc.or}
\title{
Bayesian Network construction using a hybrid of MMPC and PC
}

\description{
Bayesian Network construction using a hybrid of MMPC and PC.
}

\usage{
mmpc.or(x, max_k = 5, threshold = 0.01, test = "testIndFisher", backward = TRUE, 
symmetry = TRUE, ini.pvalue = NULL)
}

\arguments{
\item{x}{
A matrix with the variables. The user must know if they are continuous or if they are categorical. If you have a matrix with categorical data, 
i.e. 0, 1, 2, 3 where each number indicates a category, the minimum number for each variable must be 0. data.frame is also supported, as the 
dataset in this case is converted into a matrix.   
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details of SES or MMPC). 
}
\item{threshold}{
Threshold ( suitable values in (0, 1) ) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
The conditional independence test to use. Default value is "testIndFisher". This procedure allows for "testIndFisher", "testIndSPearman" for continuous variables and "gSquare" for categorical variables. 
}
\item{backward}{
If TRUE, the backward (or symmetry correction) phase will be implemented. This removes any falsely included variables in the parents and children set of the target variable. It calls the \code{link{mmpcbackphase}} 
for this purpose. For perm.ses and wald.ses this is not yet applicable.
}
\item{symmetry}{
In order for an edge to be added, a statistical relationship must have been found from both directions. If you want this symmetry correction to take place, leave this boolean variable to TRUE. If you set it to FALSE, then if 
a relationship between Y and X is detected but not between X and Y, the edge is still added. 
}
\item{ini.pvalue}{
This is a list with the matrix of the univariate p-values. If you want to run mmhc.skel again, the univariate associations need not be calculated again.
}
}

\details{
The MMPC is run on every variable. The backward phase (see Tsamardinos et al., 2006) can then take place. After all variables have been used, the matrix is checked for inconsistencies and they are corrected if you want. 
The "symmetry" argument. Do you want the edge to stay if it was discovered from both variables when they were considered as responses? 
}

\value{
A list including:
\item{ini.pvalue}{
A matrix with the p-values of all pairwise univariate assocations.
}
\item{kapa}{
The maximum number of conditioning variables ever observed. 
}
\item{ntests}{
The number of tests MMPC (or SES) performed at each variable.
}
\item{info}{
Some summary statistics about the edges, minimum, maximum, mean, median number of edges.
}
\item{density}{
The number of edges divided by the total possible number of edges, that is #edges / \eqn{n(n-1)/2}, where \eqn{n} is the number of variables.
}
\item{runtime}{
The run time of the skeleton phase of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{runtime.or}{
The run time of the PC orientation rules. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{Gini}{
The adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that i and j have an edge between them.
}
\item{G}{
The final adjaceny matrix with the orientations. If G[i, j] = 2 then G[j, i] = 3. This means that there is an arrow from node i to node j. If G[i, j] = G[j, i] = 0; there is no edge between nodes i and j. 
If G[i, j] = G[j, i] = 1; there is an (undirected) edge between nodes i and j.
}
\item{sepset}{
A list with the separating sets for every value of k.
}

Bear in mind that the values can be extracted with the $ symbol, i.e. this is an S3 class output. 
}

\references{
Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.

Spirtes P., Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, 3nd edition.
}

\author{
Michail Tsagris
  
R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
  %%  ~~further notes~~
    %}
\seealso{
\code{ \link{pc.skel}, \link{pc.or}, \link{corfs.network} }
}

\examples{
y <- rdag2(500, p = 20, nei = 3)
ind <- sample(1:20, 20)
x <- y$x[, ind]
a1 <- mmpc.or(x, max_k = 3, threshold = 0.01, test = "testIndFisher" ) 
b <- pc.skel( x, alpha = 0.01 ) 
a2 <- pc.or(b)
}

\keyword{ MMHC algorithm }