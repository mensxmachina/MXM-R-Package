\name{Check whether a directed graph is acyclic}
\alias{is.dag}
\title{
Check whether a directed graph is acyclic
}

\description{
Check whether a directed graph is acyclic.
}

\usage{
is.dag(dag)
}

\arguments{
\item{dag}{
A square matrix representing a directed graph which contains either 0, 1 or 0, 2, and 3. In the first canse where G[i, j] = 1, means there is an arrow from node i to node j. 
In the second case G[i, j] = 2 and G[j, i] = 3 means that there is an arrow from node i to node j, where the 2 iindcates the arrohead and the 3 inducates the arrowtail. 
}
}

\details{
The topological sort is performed. If it cannot be performed, NAs are returned. Hence, the functions checks for NAs.
}

\value{
A logical value, TRUE if the matrix represents a DAG and FALSE otherwise.
}

\references{
Chickering, D.M. (1995). A transformational characterization of equivalent Bayesian network structures. Proceedings of the 11th Conference on Uncertainty in Artificial Intelligence, Montreal, Canada, 87-98. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{topological_sort}, \link{dag2eg}, \link{pc.or} 
}
}

\examples{
# simulate a dataset with continuous data
# simulate a dataset with continuous data
G <- rdag(100, 20, 0.3)$G
is.dag(G)  ## TRUE 
}

\keyword{ Directed acyclic graph }
\keyword{ check for cycles }