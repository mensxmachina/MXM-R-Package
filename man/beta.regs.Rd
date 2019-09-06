\name{Many simple beta regressions}
\alias{beta.regs}
\alias{perm.betaregs}
\alias{wald.betaregs}
\title{
Many simple beta regressions.
}

\description{
Many simple beta regressions.
}

\usage{
beta.regs(target, dataset, wei = NULL, check = FALSE, logged = FALSE, ncores = 1) 

perm.betaregs(target, dataset, wei = NULL, check = FALSE, logged = FALSE, 
threshold = 0.05, R = 999, ncores = 1)

wald.betaregs(target, dataset, wei = NULL, check = FALSE, logged = FALSE, ncores = 1) 
}

\arguments{
\item{target}{
The target (dependent) variable. It must be a numerical vector with integers. 
}
\item{dataset}{
The indendent variable(s). It can be a vector, a matrix or a dataframe with continuous only variables, a data frame with mixed or only categorical variables.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{check}{
A boolean variable indicating whether to check for variables with identical values. The default is FALSE.
}
\item{logged}{
A boolean variable; it will return the logarithm of the p-value if set to TRUE.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{R}{
The number of permutations, set to 999 by default. There is a trick to avoind doing all permutations. As soon as the number of times the permuted test statistic is more than the observed test statistic is more than 50 (in this example case), the p-value has exceeded the signifiance level (threshold value) and hence the predictor variable is not significant. There is no need to continue do the extra permutations, as a decision has already been made. 
}
\item{ncores}{
The number of cores to use. The default value is 1.
}
}

\details{
Many simple beta regressions are fitted. 
}

\value{
A matrix with the test statistic values, their relevant (logged) p-values and the BIC values. 
}

\references{
Ferrari S.L.P. and Cribari-Neto F. (2004). Beta Regression for Modelling Rates and Proportions. Journal of Applied Statistics, 31(7): 799-815.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{beta.mod}, \link{testIndBeta}, \link{reg.fit}, \link{ridge.reg} }
}

\examples{
y <- rbeta(80, 5, 3)
x <- matrix( rnorm(80 * 7), ncol = 7)
a1 <- beta.regs(y, x)
a2 <- perm.betaregs(y, x[, 1:4], R = 499)
}

\keyword{ Regression models }
\keyword{ beta distribution }