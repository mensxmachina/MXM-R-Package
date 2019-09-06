\name{Interactive plot of an (un)directed graph}
\alias{plotnetwork}
\title{
Interactive plot of an (un)directed graph
}

\description{
Interactive plot of an (un)directed graph.
}

\usage{
plotnetwork(G, titlos)
}

\arguments{
\item{G}{
The adjacency matrix as produced from \code{\link{mmhc.skel}}, \code{\link{pc.skel}}, \code{\link{pc.con}} or any other algorithm. 
This can correspond to an undirected, partially directed or a completely directed graph. 
}
\item{titlos}{A character argument specifying the title of the graph, for example "PC network". 
}
}

\details{
This visualises the directed graph. 
}

\value{
The plot of the directed graph. This is interactive, in the sense that the user can "play" with it. Move the nodes, zoom it, strectch it etc.
}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{mmhc.skel}, \link{nei}, \link{pc.skel}, \link{mb} }
}

\examples{
# simulate a dataset with continuous data
dataset <- matrix( runif(500 * 20, 1, 100), nrow = 500 ) 
a <- mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher", 
nc = 1) 
plotnetwork(a$G)
plotnetwork(a$G, titlos = "DAG skeleton")
}

\keyword{ Network plot }
\keyword{ directed graph }
\keyword{ interactive graph }