\name{Data simulation from a DAG}
\alias{rdag}
\alias{rdag2}
\alias{rmdag}
\title{
Simulation of data from DAG (directed acyclic graph)
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent 
from 'TARGET' given a conditioning set CS. This test is based on the log likelihood ratio test.
}

\usage{
rdag(n, p, s, a = 0, m, A = NULL, seed = FALSE) 
rdag2(n, A = NULL, p, nei, low = 0.1, up = 1) 
rmdag(n, A = NULL, p, nei, low = 0.1, up = 1) 
}

\arguments{
\item{n}{
A number indicating the sample size. 
}
\item{p}{
A number indicating the number of nodes (or vectices, or variables).
}
\item{nei}{
The average number of neighbours.
}
\item{s}{
A number in \eqn{(0, 1)}. This defines somehow the sparseness of the model. It is the probability that a node has an edge. 
}
\item{a}{
A number in \eqn{(0, 1)}. The defines the percentage of outliers to be included in the simulated data. If \eqn{a=0}, no outliers are generated. 
}
\item{m}{
A vector equal to the number of nodes. This is the mean vector of the normal distribution from which the data are to be generated. This is used only when \eqn{a>0} 
 so as to define the mena vector of the multivariate normal from which the outliers will be generated. 
}
\item{A}{
If you already have an an adjacency matrix in mind, plug it in here, otherwise, leave it NULL. 
}
\item{seed}{
If seed is TRUE, the simulated data will always be the same. 
}
\item{low}{
Every child will be a function of some parents. The beta coefficients of the parents will be drawn uniformly from two numbers, low and up. See details for more information on this.
}
\item{up}{
Every child will be a function of some parents. The beta coefficients of the parents will be drawn uniformly from two numbers, low and up. See details for more information on this. 
}
}

\details{
In the case where no adjacency matrix is given, an \eqn{p \times p} matrix with zeros everywhere is created. 
very element below the diagonal is is replaced by random values from a Bernoulli distribution with probability of success equal to s. 
This is the matrix B. Every value of 1 is replaced by a uniform value in \eqn{0.1, 1}. This final matrix is called A. 
The data are generated from a multivariate normal distribution with a zero mean vector and covariance matrix equal to 
\eqn{\left({\bf I}_p- A\right)^{-1}\left({\bf I}_p- A\right)}, where \eqn{{\bf I}_p} is the \eqn{p \times p} identiy matrix. 
If a is greater than zero, the outliers are generated from a multivariate normal with the same covariance matrix and mean vector the one 
specified by the user, the argument "m". The flexibility of the outliers is that you cna specifiy outliers in some variables only or in all of them. For example, m = c(0,0,5) introduces outliers in the third variable only, whereas m = c(5,5,5) introduces outliers in all variables. 
The user is free to decide on the type of outliers to include in the data.

For the "rdag2", this is a different way of simulating data from DAGs. The first variable is normally generated. Every other variable can be a function of some previous ones. Suppose now that the i-th variable is a child of 4 previous variables. We need for coefficients \eqn{b_j} to multiply the 4 variables and then generate the i-th variable from a normal with mean \eqn{\sum_{j=1}b_j X_j} and variance 1. The \eqn{b_j} will be either positive or negative values with equal probability. Their absolute values ranges between "low" and "up". The code is accessible and you can see in detail what is going on. In addition, every generated data, are standardised to avoid numerical overflow.

The "rmdag" generates data from a BN with continous, ordinal and binary data in proportions 50\%, 25\% and 25\% resepctively on average. This was used in the experiments run by Tsagris et al. (2017). If you want to generate data and then use them in the "pcalg" package with the function "ci.fast2" or "ci.mm2" you should transform the resulting data into a matrix. The factor variables must becomw numeric starting from 0. See the examples for more on this. 
}


\value{
A list including:
\item{nout}{
The number of outliers.
}
\item{G}{
The adcacency matrix used. For the "rdag" if G[i, j] = 2, then G[j, i] = 3 and this means that there is an arrow from j to i. For the "rdag2" and "rmdag" the entries are either G[i, j] = G[j, i] = 0 (no edge) or G[i, j] = 1 and G[j, i] = 0 (indicating i -> j). 
}
\item{A}{
The matrix with the with the uniform values in the interval \eqn{0.1, 1}. This is returned only by "rdag". 
}
\item{x}{
The simulated data. 
}
}

\references{
Tsagris M., Borboudakis G., Lagani V. and Tsamardinos I. (2018). Constraint-based Causal Discovery with Mixed Data. 
International Journal of Data Science and Analytics. 

https://link.springer.com/content/pdf/10.1007%2Fs41060-018-0097-y.pdf
 
http://nugget.unisa.edu.au/CD2017/papersonly/constraint-based-causal-r1.pdf

Spirtes P.,  Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, 3nd edition.

Colombo, Diego, and Marloes H. Maathuis (2014). Order-independent constraint-based causal structure learning. The Journal of Machine Learning Research 15(1): 3741--3782.
}

\author{
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}


\seealso{
\code{\link{pc.skel}, \link{pc.or}, \link{ci.mm}, \link{mmhc.skel} }
}

\examples{
y <- rdag(100, 20, 0.2)
x <- y$x
tru <- y$G 

mod <- pc.con(x)
b <- pc.or(mod)
plotnetwork(tru) 
dev.new()
plotnetwork(b$G)

}

\keyword{ Directed acylic graph }
\keyword{ simulation of random values }
