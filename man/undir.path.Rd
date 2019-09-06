\name{Undirected path(s) between two nodes}
\alias{undir.path}
\title{
Undirected path(s) between two nodes
}

\description{
Undirected path(s) between two nodes.
}

\usage{
undir.path(G, y, x) 
}

\arguments{
\item{G}{
An adjacency matrix where G[i,j] = G[j, i] = 1 means there is an edge between nodes i and j. If G[i, j] = G[j, i] = 0 there is no edge between them.     
}
\item{y}{
A numerical value indicating the first node, it has to be a number between 1 and the maximum number of variables.
}
\item{x}{
A numerical value indicating the second node, it has to be a number between 1 and the maximum number of variables. The order of the nodes does not make a difference.
}
}

\details{
The algorithm finds all the nodes between the two nodes. It finds all paths between the two chosen nodes.
}

\value{
A vector with the two nodes and all nodes between them in the case of connecting nodes. Otherwise, a matrix with the neighbours of each node.
}

\references{
Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{SES}, \link{MMPC}, \link{pc.skel} }
}

\examples{
# simulate a dataset with continuous data
set.seed(1234)
dataset <- matrix(runif(1000 * 10, 1, 100), nrow = 1000 ) 
G <- pc.con(dataset)$G
plotnetwork(G)
undir.path(G, 3, 4)
undir.path(G, 1, 3)
}

\keyword{ Undirected path }