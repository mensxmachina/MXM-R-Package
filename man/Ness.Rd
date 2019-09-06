\name{Effective sample size for G^2 test in BNs with case control data}
\alias{Ness}
\title{
Effective sample size for G^2 test in BNs with case control data
}

\description{
Effective sample size for G^2 test in BNs with case control data.
}

\usage{
Ness(propNt, N, K = 10000) 
}

\arguments{
\item{propNt}{
A numerical vector with the proportions (distribution) of the (single) selection variable.   
}
\item{N}{
The sample size of the data. 
}
\item{K}{
The number of repetitions to be used for estimating the effective sample size.
}
}

\details{
When dealing with case control data, spurious correlations or relationships arise. To deal with this one way is to adjust the sample size used in the G^2 test statistic. This function does exactly this, estimates the effective sample size as per the Borboudakis and Tsamardinos (2012) suggestion. The idea is that after learning the skeleton with the usual G^2 test, one should go to the edges and perform a conditional G^2
}

\value{
The estimated effective sample size.
}

\references{
Borboudakis G. and Tsamardinos I. (2015). Bayesian Network Learning with Discrete Case-Control Data. 31st Conference on Uncertainty in Artificial Intelligence (UAI), 151-160.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{SES}, \link{MMPC}, \link{testIndLogistic} 
}
}

\examples{
Ness(c(0.3, 0.7), N = 1000, K = 10000) 
}

\keyword{ Area under the curve }
\keyword{ receiver operating curve }
