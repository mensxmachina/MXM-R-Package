\name{Topological sort of a DAG}
\alias{topological_sort}
\title{
Topological sort of a DAG
}

\description{
Topological sort of a DAG.
}

\usage{
topological_sort(dag)
}

\arguments{
\item{dag}{
A square matrix representing a directed graph which contains 0s and 1s. If G[i, j] = 1 it 
means there is an arrow from node i to node j. When there is no edge between nodes i and j if G[i, j] = 0.  
}
}

\details{
The function is an R translation from an old matlab code.
}

\value{
A vector with numbers indicating the sorting. If the matrix does not correspond to a DAG, NA will be returned. 
}

\references{
Chickering, D.M. (1995). A transformational characterization of equivalent Bayesian network structures. 
Proceedings of the 11th Conference on Uncertainty in Artificial Intelligence, Montreal, Canada, 87-98. 
}

\author{
Ioannis Tsamardinos and Michail Tsagris

R implementation and documentation:  Ioannis Tsamardinos <tsamard@csd.uoc.gr> and Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{plotnetwork}, \link{nei}, \link{pc.or} 
}
}

\examples{
# simulate a dataset with continuous data
# simulate a dataset with continuous data
G = rdag(100, 10, 0.3)$G
G[G == 2] <- 1
G[G == 3] <- 0
topological_sort(G)
}

\keyword{ Topological sort }
\keyword{ DAG }