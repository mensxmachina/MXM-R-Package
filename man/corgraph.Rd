\name{Graph of unconditional associations}
\alias{corgraph}

\title{
Graph of unconditional associations
}

\description{
Calcualtes the graph of unconditional associations. If the correlation (Pearson, Spearman) or the \eqn{G^2} test of independence, between pairs of continuous or categorical variables respectively is not statistically significant, there is no edge between the two respective nodes.
}

\usage{
corgraph(dataset, test = "testIndFisher", threshold = 0.01)
}

\arguments{
\item{dataset}{
A matrix with the variables. The user must know if they are continuous or if they are categorical. If you have a matrix with categorical data, i.e. 0, 1, 2, 3 where each number indicates a category, the minimum number for each variable must be 0.
}
\item{test}{
The conditional independence test to use. Default value is "testIndFisher". This procedure allows for "testIndFisher", "testIndSPearman" for continuous variables and "gSquare" for categorical variables. 
}
\item{threshold}{
Threshold ( suitable values in (0, 1) ) for assessing p-values significance. Default value is 0.05.
}
}

%\details{
%}

\value{
A list including:
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{stat}{
A matrix with the test statistics. 
}
\item{pvalue}{
A matrix with the p-values. 
}
\item{G}{
The adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that i and j have an edge between them.
}
}

\author{
Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{pc.skel}, \link{mmhc.skel}, \link{corfs.network}, \link{bn.skel.utils} }
}

\examples{
# simulate a dataset with continuous data
y <- rdag2(500, p = 20, nei = 3)
x <- y$x
a <- mmhc.skel(x, max_k = 5, threshold = 0.01, test = "testIndFisher" ) 
b <- pc.skel( x, alpha = 0.01 ) 
d <- corgraph(x, test = "testIndFisher", threshold = 0.01) 
}

