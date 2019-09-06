\name{Probability residual of ordinal logistic regreession}
\alias{ord.resid}
\title{
Probability residual of ordinal logistic regreession
}

\description{
Probability residual of ordinal logistic regreession.
}

\usage{
ord.resid(y, est)
}

\arguments{
\item{y}{
An ordered factor variable or a numerical vector.
}
\item{est}{
A matrix with the fitted values of an ordinal logistic regression model. 
}
}

\details{
The probability residual of an ordinal logistic regression model is calculated (Li and Shepherd, 2012). 
It is a vector, irrespective of how many categories there are.
}

\value{
A vector with the probability residuals.
}

\references{
Li C., & Shepherd B. E. (2012). A new residual for ordinal outcomes. Biometrika, 99(2): 473--480.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{testIndOrdinal}, \link{ordinal.reg} 
}
}

\examples{
library(MASS)
y <- factor( rbinom(400, 3, 0.6), ordered = TRUE )
x <- rnorm(400)
mod <- MASS::polr(y ~ x) 
res <- ord.resid(y, mod$fitted.values)
}

