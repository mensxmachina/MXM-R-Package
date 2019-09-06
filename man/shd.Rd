\name{Structural Hamming distance between two partially oriented DAGs}
\alias{shd}
\title{
Structural Hamming distance between two partially oriented DAGs
}

\description{
Structural Hamming distance between two partially oriented DAGs.
}

\usage{
shd(est, true) 
}

\arguments{
\item{est}{
The first (partially oriented) DAG. This could also be the estimated DAG.
}
\item{true}{
The second (partially oriented) DAG. This could also be the equivalence class of the true DAG.  
}
}

\details{
The structural Hamming distance as proposed by Tsamardinos et al. (2006) is calculated and returned. The cases are listed below
\tabular{lll}{ 
\bold{True} \tab  \bold{Estimated}  \tab  \bold{Penalty}  \cr 
   -   \tab      \tab   1  \cr   
       \tab  -   \tab   1  \cr   
   ->  \tab      \tab   1  \cr   
       \tab  <-  \tab   1  \cr   
   ->  \tab  -   \tab   1  \cr   
   -   \tab  <-  \tab   1  \cr   
   ->  \tab  <-  \tab   1  \cr        
}
}

\value{
A list including
\item{mat}{
A table with the agreements and disagreements between the two DAGs.
}
\item{shd}{
The structural Hamming distance.
}
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
\code{ \link{pc.skel}, \link{pc.or}, \link{mmhc.skel}, \link{plotnetwork} 
}
}

\examples{
y <- rdag(1000, 20, 0.2)
tru <- y$G 
mod <- pc.skel(y$x)
a <- pc.or(mod)
shd( a$G, dag2eg(tru) )
}

\keyword{ Structurall Haming distance }
\keyword{ distance between DAGs }
