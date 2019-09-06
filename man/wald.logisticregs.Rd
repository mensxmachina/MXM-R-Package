\name{Many Wald based tests for logistic and Poisson regressions with continuous predictors}
\alias{wald.logisticregs}
\alias{wald.poissonregs}
\title{
Many Wald based tests for logistic and Poisson regressions with continuous predictors
}

\description{
Many Wald based tests for logistic and Poisson regressions with continuous predictors.
}

\usage{
wald.logisticregs(y, x, tol = 1e-09, wei = NULL, check = FALSE, logged = FALSE, 
ncores = 1)  
wald.poissonregs(y, x, tol = 1e-09, wei = NULL, check = FALSE, logged = FALSE, 
ncores = 1)  
}

\arguments{
\item{y}{
A vector with either 0s and 1 (logistic regression) or discrete data, counts (Poisson regression).
}
\item{x}{
A data.frame, the predictor variables. If you have no categorical variables, the fucntion will still work but it's better to use the \code{\link{score.glms}} because it is faster. 
}
\item{tol}{
The tolerance value to stop the Newton-Raphson iterations. It is set to 1e-09 by default. 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{check}{
A boolean variable indicating whether to chekc for variables with identical values. The defauls is FALSE.
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
\item{ncores}{
How many to cores to useq the default value is 1. 
}
}

\details{
Instead of using R built-in function \code{\link{glm}} we implemented the newton-Raphson in order to avoid unnecessary calculations. The functions are much faster.  
}

\value{
A matrix with three columns, the test statistic, its associated (logged) p-value and the BIC of each model.   
}

\references{
Draper, N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition.

McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{univregs}, \link{perm.univregs}
}
}

\examples{
## 100 variables, hence 100 univariate regressions are to be fitted
x <- matrix( rnorm(200 * 100), ncol = 100 )
y <- rpois(200, 4)
a <- wald.poissonregs(y, x)
b <- univregs(y, x, test = testIndPois)
cor(exp(a[, 2]), exp(b$pvalue) )
}

\keyword{ Poisson regressions }
\keyword{ logistic regressions }
\keyword{ Wald test }
