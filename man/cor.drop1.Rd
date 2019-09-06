\name{Drop all possible single terms from a model using the partial correlation}
\alias{cor.drop1}
\title{
Drop all possible single terms from a model using the partial correlation
}

\description{
Drop all possible single terms from a model using the partial correlation.
}

\usage{
cor.drop1(y, x, logged = FALSE)
}

\arguments{
\item{y}{
A numerical vector with the response variable.
}
\item{x}{
A numerical matrix or a data.frame with the predictor variables. If is is a matrix it is internally transformed into a data.frame form, hence the user is advised to supply a data.frame in order to save some time. If the number of columns (variables) is higher than the number of rows (observations) the function will simply not work.
}
\item{logged}{
If you want the p-values be returned leave this FALSE. If it is TRUE their logarithm is returned.
}
}

\details{
This uses R's command \code{\link{drop1}} and modifies it so as to calculate the p-value using Fisher's conditional independence test.
}

\value{
A matrix with two columns, the test statistic values and its associated p-value.
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
\code{ \link{glm.bsreg}, \link{fbed.reg}, \link{mmpcbackphase} 
}
}

\examples{
y <- rnorm(500)
x <- matrix( rnorm(500 * 20), ncol = 20)
cor.drop1(y, x)
}

\keyword{ Area under the curve }
\keyword{ receiver operating curve }
