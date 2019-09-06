\name{Ancestors and descendants of a node in a directed graph}
\alias{findDescendants}
\alias{findAncestors}
\title{
Returns and plots, if asked, the descendants or ancestors of one or all node(s) (or variable(s))
}

\description{
Returns and plots, if asked, the descendants or ancestors of one or all node(s) (or variable(s))
}

\usage{
findDescendants(G, node = NULL, graph = FALSE)
findAncestors(G, node = NULL, graph = FALSE)
}

\arguments{
\item{G}{
The graph matrix as produced from \code{\link{pc.or}} or any other algorithm which produces directed graphs. 
}
\item{node}{
A numerical value indicating the node (or variable) whose descendants are to be returned. 
}
\item{graph}{
A boolean variable. If TRUE the relevant graph will appear (if there are descendants). 
}
}

\details{
The functions searches for the descendants of some node. This is an S3 class output. 
}

\value{
\item{isAnc}{
A matrix of the same dimensions as the original graph matrix with 0s and 1s. isAnc[i, j] = 1 indicates that the i-th node is an ancestor of the j-th node. 
If the argument "node" is NULL, only this matrix will be returned. 
}
\item{Ganc}{
A matrix of dimensions equal to the number of descendants of the node with 0s and 1s, if the argument "node" is not NULL.  
}
\item{anc}{
The descendants of the node if the argument "node" is not NULL.
}
}

%\references{
% 
%}

\author{
 Anna Roumpelaki

R implementation and documentation: Anna Roumpelaki <anna.roumpelaki@gmail.com>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{plotnetwork}, \link{nei}, \link{mb}, \link{pc.or} 
}
}

\examples{
# simulate a dataset with continuous data
y = rdag(1000, 10, 0.3)
tru = y$G 
x = y$x
mod = pc.con(x)
G = pc.or(mod)$G
plotnetwork(G)
findDescendants(G, 4, graph = FALSE)
findAncestors(G, 4, graph = FALSE)
findAncestors(G)
}

\keyword{ Ancestors }
\keyword{ Descendants }
