\name{Check Markov equivalence of two DAGs}
\alias{equivdags}
\title{
Check Markov equivalence of two DAGs 
}

\description{
Check Markov equivalence of two DAGs.
}

\usage{
equivdags(g1, g2)
}

\arguments{
\item{g1}{
The matrix of a DAG or a partially directed graph as produced from \code{\link{pc.or}} or any other algorithm. 
}
\item{g2}{
The matrix of a DAG or a partially directed graph as produced from \code{\link{pc.or}} or any other algorithm. 
}
}

\details{
Two DAGs are Markov equivalent if a) they have the same adjancencies (regardlsee of the mark, arrowhead, tail or nothing) and  b) they have the same unshielded colliders.
}

\value{
A list including:
\item{apofasi}{
A boolean variable, TRUE of FALSE.  
}
\item{mes}{
A message specyfying the result, the dimensions of the adjacency matrices do not match for example, or the number of adjancencies 
is not the same, they do not share the same unshilded colliders, or they are Markov equivalent. 
}
}

\references{
Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.
}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{pc.or}, \link{mmhc.skel}, \link{pc.con} }
}

\examples{
y <- rdag(1000, 10, 0.3)
tru <- y$G 
x <- y$x
mod <- pc.con(x)

eg <- dag2eg(y$G) ## make it essential graph first 
est <- pc.or(mod)$G

equivdags(est, tru)
}

\keyword{ DAGs }
\keyword{ Markov equivalence }
