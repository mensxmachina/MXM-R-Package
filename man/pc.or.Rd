\name{Orientation rules for the PC algorithm}
\alias{pc.or}
\title{
The orientations part of the PC algorithm.
}

\description{
The function takes the outcome of the PC algorithm, as produced by \code{\link{pc.skel}} or \code{\link{pc.con}} and performes 
the 4 orientation rules. A graph is also possible to visualize.
}

\usage{
pc.or(mod) 
}

\arguments{
\item{mod}{
An object with the results of the PC algorithm, as produced by \code{\link{pc.skel}} or \code{\link{pc.con}}. There was an extra 
argument for plotting the skeleton but it does not work with the current visualisation packages, hence we removed the argument. 
Use \code{\link{plotnetwork}} to plot the skeleton.   
}
}

\details{
After having calculated the skeleton of the PC algorithm one may wants to perform orientations, leading to causal relationships. 
The rules as stated in Spirtes, Glymour and Scheines (2001) are
\enumerate{
\item \bold{Rule 0}. For each triple of vertices X, Y, Z such that the pair X, Y and the pair Y, Z are each adjacent in C but the pair X, Z 
are not adjacent in C, orient X - Y - Z as X -> Y <- Z if and only if Y is not in Sepset(X, Z). 
\item \bold{Rule 1}. If A -> B, B and C are adjacent, A and C are not adjacent, and there is no arrowhead at B, then orient B - C as B -> C.
\item \bold{Rule 2}. If there is a directed path from A to B, and an edge between A and B, then orient A - B as A -> B.
\item \bold{Rule 3}. If A -> B <- C, A - D - C, A and C are not adjacent, and D - B, then orient D - B as D -> B. 
}
The first rule is applied once. Rules 2-4 are applied repeatedly until no more edges can be oriented. If when a rule is applied and a cycle is 
detected, the rule is cancelled. Also, when applying Rules 1-3 we try to avoid the creation of new v-structures (X -> Y <- Z). 
}

\value{
A list including:
\item{Gini}{
The initial adjacency matrix, no orientations. This is the matrix produced by \code{\link{pc.skel}} or \code{\link{pc.con}}. 
}
\item{G}{
The final adjaceny matrix with the orientations. If G[i, j] = 2 then G[j, i] = 3. This means that there is an arrow from node i to node j. 
If G[i, j] = G[j, i] = 0; 
there is no edge between nodes i and j. If G[i, j] = G[j, i] = 1; there is an (undirected) edge between nodes i and j.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the 
third element is the elapsed time.
}
}

\references{
Spirtes P.,  Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, 3nd edition.

Zhang, Jiji. (2008). On the completeness of orientation rules for causal discovery in the presence of latent confounders and selection bias. 
Artificial Intelligence 172(16): 1873--1896.

Tsagris M. (2019). Bayesian network learning with the PC algorithm: an improved and correct variation. 
Applied Artificial Intelligence 33(2): 101-123.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{pc.con}, \link{pc.skel}, \link{mmhc.skel}, \link{is.dag}, \link{mb} 
}
}

\examples{
# simulate a dataset with continuous data
y <- rdag2(2000, p = 20, nei = 3)
ind <- sample(1:20, 20)
tru <- y$G[ind, ind] 
x <- y$x[, ind]
mod <- pc.con(x)
mod$runtime

plotnetwork(tru) 

b <- pc.or(mod)
plotnetwork(b$G)

plotnetwork( dag2eg(tru) )  ## essential graph
plotnetwork(b$G)
}

\keyword{ PC algorithm }