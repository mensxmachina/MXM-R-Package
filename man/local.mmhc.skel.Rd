\name{Skeleton (local) around a node of the max-min hill-climbing (MMHC) algorithm}
\alias{local.mmhc.skel}
\title{
Skeleton (local) around a node of the MMHC algorithm
}

\description{
The local skeleton of a Bayesian network around a node produced by MMHC. No orientations are involved.
}

\usage{
local.mmhc.skel(dataset, node, max_k = 3, threshold = 0.05, test = "testIndFisher") 
}

\arguments{
\item{dataset}{
A matrix with the variables. The user must know if they are continuous or if they are categorical. If you have a matrix with categorical data, i.e. 0, 1, 2, 3 where each number indicates a category, the minimum number for each variable must be 0. data.frame is also supported, as the dataset in this case is converted into a matrix.   
}
\item{node}{
A number between 1 and the number of columns of the dataset. The local network (edges only) will be built around this node. At first the parents and children of this node are identified and then their parents and children. No inconsistencies correction whastsoever is attempted. A variable detected by the node, but the node was not detected by that variable.   
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details of SES or MMPC). 
}
\item{threshold}{
Threshold ( suitable values in (0, 1) ) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
The conditional independence test to use. Default value is "testIndFisher". This procedure allows for "testIndFisher", "testIndSPearman" for continuous variables and "gSquare" for categorical variables. Or in general, if the dataset is a data.frame with different types of data, leave this NULL. See also \code{link{MMPC}} for the automatic choice of tests.
}
}

\details{
The MMPC is run on the user specific variable. The backward phase (see Tsamardinos et al., 2006) takes place automatically. Then, the MMPC is run on the parents and children of that variable. If the node variable is not detected by a variable, this variable is not removed though. 
}

\value{
A list including:
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{ntests}{
The number of tests MMPC (or SES) performed at each variable.
}
\item{res}{
A list with the parents and children of each variable. The first element is the parents and children of the node variable.
}
\item{Gloc}{
The local adjancency matrix. A value of 1 in G[i, j] may not appear appear in G[j, i] also, indicating that variable j was discovered as a possible parent or child of node i, but not the covnerse. The usual MMHC (\code{\link{mmhc.skel}}) removes the edge between them as this is an inconsistency. 
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
\code{ \link{mmhc.skel}, \link{pc.skel}, \link{pc.or}, \link{corfs.network} }
}

\examples{
# simulate a dataset with continuous data
dataset <- matrix(runif(500 * 30, 1, 100), nrow = 500 ) 
a1 <- mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher") 
a2 <- local.mmhc.skel(dataset, 4)
a1$runtime  
a2$runtime 

dataset <- rdag2(500, p = 20, nei = 3)$x
a1 <- mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher") 
a2 <- local.mmhc.skel(dataset, 5)
a1$runtime  
a2$runtime 

}

\keyword{ MMHC algorithm }