\name{Transformation of a DAG into an essential graph}
\alias{dag2eg}
\title{
Transforms a DAG into an essential graph
}

\description{
Transforms a DAG into an essential graph.
}

\usage{
dag2eg(dag, type = NULL) 
}

\arguments{
\item{dag}{
The graph matrix as produced from \code{\link{pc.or}} or any other algorithm which produces directed graphs. A DAG in general. 
}
\item{type}{
This can be either NULL or 1 or 2. type = 1 means that the matrix contains 0, 1, 2, 3 where G[i, j] = g[j, i] = 0, means there is no edge between nodes 
i and j, G[i, j] = g[j, i] = 1, there is an edge between nodes i and j and G[i, j] = 2 and G[j, i] = 3 means that there is an arrow from node i to node j. 
If type 2, the matrix contains 0 for no edge and 1 for a directed edge. In this case, G[i,j]=1 and G[j,i]=0 means that there is an arrow from node i to 
node j. If you are not sure of what you have, just leave it NULL, the function will check to which case your matrix belongs.
}
}

\details{
The function is an R translation from an old matlab code.
}

\value{
The matrix of the essential graph.
}

\references{
Chickering, D.M. (1995). A transformational characterization of equivalent Bayesian network structures. Proceedings of the 11th Conference on Uncertainty in Artificial Intelligence, Montreal, Canada, 87-98. 
}

\author{
Ioannis Tsamardinos and Michail Tsagris

R implementation and documentation: Ioannis Tsamardinos <tsamard@csd.uoc.gr> and  and Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{plotnetwork}, \link{is.dag}, \link{topological_sort} 
}
}

\examples{
# simulate a dataset with continuous data
y <- rdag(1000, 10, 0.3)
tru <- y$G 
eg <- dag2eg(tru)
par( mfrow = c(1, 2) )
plotnetwork(tru)
plotnetwork(eg)
}

\keyword{ Essential graph }
\keyword{ equivalence class }