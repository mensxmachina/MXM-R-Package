\name{Many simple zero inflated Poisson regressions}
\alias{zip.regs}
\alias{perm.zipregs}
\alias{wald.zipregs}
\title{
Many simple zero inflated Poisson regressions.
}

\description{
Many simple zero inflated Poisson regressions.
}

\usage{
zip.regs(target, dataset, wei = NULL, check = FALSE, logged = FALSE, ncores = 1) 

perm.zipregs(target, dataset, wei = NULL, check = FALSE, logged = FALSE, R = 999, 
threshold = 0.05, ncores = 1)

wald.zipregs(target, dataset, wei = NULL, check = FALSE, logged = FALSE, ncores = 1) 
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
A boolean variable indicating whether to chekc for variables with identical values. The defauls is FALSE.
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
\item{R}{
The number of permutations, set to 999 by default. There is a trick to avoind doing all permutations. As soon as the number of times the permuted test statistic is more than the observed test statistic is more than 50 (in this example case), the p-value has exceeded the signifiance level (threshold value) and hence the predictor variable is not significant. There is no need to continue do the extra permutations, as a decision has already been made. 
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{ncores}{
The number of cores to use. The default value is 1.
}
}

\details{
Many simple zero inflated Poisson regressions are fitted. The permutations based approach may cause some errors sometimes, and this is due to the nature of the distribution and its maximisation process. "nlm" and "optim" are used internally. 
}

\value{
A matrix with the test statistic values, their relevant (logged) p-values and the BIC values. 
}

\references{
Lambert D. (1992). Zero-inflated Poisson regression, with an application to defects in manufacturing. Technometrics, 34(1):1-14.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{zip.mod}, \link{testIndZIP}, \link{reg.fit}, \link{ridge.reg} }
}

\examples{
y <- rpois(50, 3)
x <- matrix( rnorm(50 * 7), ncol = 7)
y[1:10] <- 0
a1 <- zip.regs(y, x)
a2 <- perm.zipregs(y, x[, 1:3], R = 499)
}

\keyword{ Regression models }
\keyword{ zero inflated Poisson }