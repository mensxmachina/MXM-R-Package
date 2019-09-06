\name{Search for triangles in an undirected graph}
\alias{triangles.search}
\title{
Search for triangles in an undirected graph
}

\description{
Search for triangles in an undirected graph.
}

\usage{
triangles.search(G) 
}

\arguments{
\item{G}{
The adjacency matrix of an undirected graph. G[i, j] = G[j, i] = 1 means there is an edge between modes i and j. Zero values indicate the absence of edges. 
}
}

\details{
The functions searches for triangles, that is for tripletes of nodes (or variables) for which X-Y, Y-Z and X-Z.
}

\value{
A matrix with thre columns. If there are no triangles, the matrix is empty. If there is at least one triangle, then each row contains three numbers, one for each node. See the examples below.
}

%\references{
%
%}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{plotnetwork}, \link{pc.skel}, \link{mmhc.skel} 
}
}

\examples{
set.seed(123)
x <- rdag2(1000, p = 20, nei = 4)$x
a <- pc.skel(x, alpha = 0.05)
plotnetwork(a$G)
triangles.search(a$G)
}

