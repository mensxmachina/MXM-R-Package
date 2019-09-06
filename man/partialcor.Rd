\name{Partial correlation between two variables}
\alias{partialcor}
\title{
Partial correlation
}

\description{
Partial correlation between two variables when a correlation matrix is given.
}

\usage{
partialcor(R, indx, indy, indz, n) 
}

\arguments{
\item{R}{
A correlation matrix.
}
\item{indx}{
The index of the first variable whose conditional correlation is to estimated.
}
\item{indy}{
The index of the second variable whose conditional correlation is to estimated.
}
\item{indz}{
The index of the conditioning variables.
}
\item{n}{
The sample size of the data from which the correlation matrix was computed.
}
}

\details{
Given a correlation matrix the function will caclulate the partial correlation between variables indx and indy conditioning on variable(s) indz. The logarithm of the p-value is also returned.  
}

\value{
The partial correlation coefficient and the logged p-value for the test of no association.
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
\code{ \link{testIndFisher}, \link{testIndSpearman}, \link{permcor}, \link{pc.con}
}
}

\examples{
r <- cor( iris[, 1:4] )
partialcor(r, 1, 2, 0, 150) 
r[1, 2]  ## the same as above

y <- as.vector( iris[, 1] )
x <- as.vector( iris[, 2] )
z <- as.vector( iris[, 3] )
e1 <- resid( lm(y ~ z) )
e2 <- resid( lm(x ~ z) )
cor(e1, e2)
partialcor(r, 1,2, 3, 150)
}

\keyword{ Partial correlation }
\keyword{ Conditioning set }
