\name{Permutation based p-value for the Pearson correlation coefficient}
\alias{permcor}
\alias{permcorrels}
\title{
Permutation based p-value for the Pearson correlation coefficient
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS.
}

\usage{
permcor(x1, x2, R = 999) 
permcorrels(y, x, R = 999)
}

\arguments{
\item{x1}{
A numerical vector.
}
\item{x2}{
A numerical vector of the same size as x1.
}
\item{y}{
A vector whose length is equal to the number of rows of x.
}
\item{x}{
This is a matrix with many variables. 
}
\item{R}{
The number of permutations to be conducted; set to 999 by default. 
}
}

\details{
This is a computational non parametric (permutation based) correlation coefficient test and is advised to be used when a small sample size is available. If you want to use the Spearman correlation instead, simply provide the ranks of x or of y and x. 
}

\value{
For the case of "permcor" a vector consisting of two values, the Pearson correlation and the permutation based p-value. 
For the "permcorrels" a vector with three values, the Pearson correlation, the test statistic value and the permutation based logged p-value.
}

\references{
Legendre Pierre (2000). Comparison of permutation methods for the partial correlation and partial Mantel tests. Journal of Statistical Computation and Simulation 67(1):37-73.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pc.skel}, \link{testIndSpearman}, \link{testIndFisher}, \link{SES}, \link{CondIndTests}}
}

\examples{
MXM::permcor(iris[, 1], iris[, 2], R = 999)
x <- matrix( rnorm(50 * 100), ncol = 100)
a <- permcorrels(iris[1:50, 1], x)
}

\keyword{ Permutation test }
\keyword{ conditional independence test }
