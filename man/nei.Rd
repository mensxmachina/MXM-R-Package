\name{Neighbours of nodes in an undirected graph}
\alias{nei}
\title{
Returns the node(s) and their neighbour(s), if there are any.
}

\description{
Returns the node(s) and their neighbour(s) of one or more nodes (if there are any).
}

\usage{
nei(G, node)
}

\arguments{
\item{G}{
The adjacency matrix of an undirected graph as produced by \code{\link{mmhc.skel}}, \code{\link{pc.skel}} or any other algorithm. 
}
\item{node}{
A vector with one or more numbers indicating the seleted node(s) (or variable(s)). 
}
}

\details{
This is a way to see the network for some given nodes. It is useful if you have many nodes and the whole network is a bit difficult to see clearly. 
}

\value{
A list object called "geit" containing the neighbours of the node(s). If there are no neighbours a message appears and no plot is presented. If the "graph" argument is set to TRUE and there are neighbours, a plot will appear. 

Bear in mind that the values can be extracted with the $ symbol, i.e. this is an S3 class output. 

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
\code{ \link{mmhc.skel}, \link{SES}, \link{MMPC} }
}

\examples{
# simulate a dataset with continuous data
set.seed(1234)
dataset <- matrix(runif(500 * 20, 1, 100), nrow = 500 ) 
G <- pc.con(dataset)$G
plotnetwork(G)
dev.new()
nei( G, c(3, 4) )
nei( G, c(1, 3) )
}

\keyword{ Network plot }
\keyword{ Neighbour nodes }