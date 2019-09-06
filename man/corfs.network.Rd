\name{Network construction using the partial correlation based forward regression or FBED}
\alias{corfs.network}
\alias{corfbed.network}
\title{
Network construction using the partial correlation based forward regression of FBED
}

\description{
Network construction using the partial correlation based forward regression or FBED. 
}

\usage{
corfs.network(x, threshold = 0.05, tolb = 2, tolr = 0.02, stopping = "BIC", 
symmetry = TRUE, nc = 1) 

corfbed.network(x, threshold = 0.05, symmetry = TRUE, nc = 1)
}

\arguments{
\item{x}{
A matrix with continuous data.
}
\item{threshold}{
Threshold ( suitable values in (0, 1) ) for assessing p-values significance. Default value is 0.05.
}
\item{tolb}{
The difference in the BIC bewtween two successive values. By default this is is set to 2. If for example, the BIC difference between two succesive models is less than 2, the process stops and the last variable, 
even though significant does not enter the model.
}
\item{tolr}{
The difference in the adjusted \eqn{R^2} bewtween two successive values. By default this is is set to 0.02. If for example, the difference between the adjusted \eqn{R^2} of two succesive models is less than 0.02, 
the process stops and the last variable, even though significant does not enter the model.
}
\item{stopping}{
The stopping rule. The "BIC" is the default value, but can change to "ar2" and in this case the adjusted \eqn{R^2} is used. If you want both of these criteria to be satisified, type "BICR2". 
}
\item{symmetry}{
In order for an edge to be added, a statistical relationship must have been found from both directions. If you want this symmetry correction to take place, leave this boolean variable to TRUE. If you set it to FALSE, then if 
a relationship between Y and X is detected but not between X and Y, the edge is still added. 
}
\item{nc}{
How many cores to use. This plays an important role if you have many variables, say thousands or so. You can try with nc = 1 and with nc = 4 for example to see the differences. If you have a multicore machine, this is a must 
option. There was an extra argument for plotting the skeleton but it does not work with the current visualisation 
packages, hence we removed the argument. Use \code{\link{plotnetwork}} to plot the skeleton.   
}
}

\details{
In the MMHC algorithm (see \code{\link{mmhc.skel}}), the MMPC or SES algorithms are run for every variable. Hence, one can use forward regression for each variable and this is what we are doing here. Partial correlation 
forward regression is very efficient, since only correlations are being calculated. 
}

\value{
A list including:
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{density}{
The number of edges divided by the total possible number of edges, that is #edges / \eqn{n(n-1)/2}, where \eqn{n} is the number of variables.
}
\item{info}{
Some summary statistics about the edges, minimum, maximum, mean, median number of edges.
}
\item{ntests}{
The number of tests MMPC (or SES) performed at each variable.
}
\item{G}{
The adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that i and j have an edge between them.
}

Bear in mind that the values can be extracted with the $ symbol, i.e. this is an S3 class output. 
}

\references{
Sanford Weisberg (2014). Applied Linear Regression. Hoboken NJ: John Wiley, 4th edition.

Draper N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition. 

Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.

Borboudakis G. and Tsamardinos I. (2017). Forward-Backward Selection with Early Dropping. 
https://arxiv.org/pdf/1705.10770.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{mmhc.skel}, \link{pc.skel} }
}

\examples{
# simulate a dataset with continuous data
dataset <- matrix(runif(600 * 30, 1, 100), nrow = 600 ) 
a1 <- mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher", 
nc = 1) 
a2 <- corfs.network(dataset, threshold = 0.05, tolb = 2, tolr = 0.02, stopping = "BIC", 
symmetry = TRUE, nc = 1) 
a1$runtime  
a2$runtime 
}

\keyword{ Network construction }
\keyword{ forward regression }