\name{Transitive closure of an adjacency matrix}
\alias{transitiveClosure}
\title{
Returns the transitive closure of an adjacency matrix
}

\description{
Returns the transitive closure of an adjacency matrix.
}

\usage{
transitiveClosure(amat)
}

\arguments{
\item{amat}{
The adjacency matrix of a graph. 
}
}

\details{
A function that computes the transitive closure of a graph. The transitive closure C(G) of a graph is a graph which contains an edge between nodes u and v whenever there is a directed path from u to v (Skiena 1990, p. 203).
http://mathworld.wolfram.com/TransitiveClosure.html
}

\value{
\item{closure}{
The transititve closure of the adjacency matrix representing a graph.
}
}

\references{
Skiena S. (1990). Implementing Discrete Mathematics: Combinatorics and Graph Theory with Mathematica. Reading, MA: Addison-Wesley 
}

\author{
 Anna Roumpelaki

R implementation and documentation: Anna Roumpelaki <anna.roumpelaki@gmail.com>
}

%\note{
%%  ~~further notes~~
%}
%\seealso{
%\code{ \link{plota}, \link{nei}, \link{pc.or} 
%}
%}

\examples{
# example adjacency matrix
# simulate a dataset with continuous data
dataset <- matrix( runif(300 * 20, 1, 100), nrow = 300 ) 
test <- pc.con( dataset, method = "pearson", alpha = 0.05 )$G 
transitiveClosure(test)
}

%\keyword{ Network plot }
%\keyword{ Neighbour nodes }