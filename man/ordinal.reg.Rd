\name{Generalised ordinal regression}
\alias{ordinal.reg}
\title{
Generalised ordinal regression
}

\description{
Generalised ordinal regression.
}

\usage{
ordinal.reg(formula, data)
}

\arguments{
\item{formula}{
An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. This is the usual formula use by many regression
models in R and other packages. Type "glm" or "formula" in R for more information. 
}
\item{data}{
A data.frame object carrying the relevant data.
}
}

\details{
Generalised ordinal regression is fitted. This means the instead of having the same coefficient for each predictor variable, they are allowed 
to vary. The usual, proportional odds, ordinal regression specifies that the lines do not cross. This one does not need the proportional odds assumption. The proportional odds assumption does not always hold in practice and is a rather restrictive model. Be careful though, you may 
end up qith negative probabilities. We do a tick to fix them, but in that case, you may have not found the optimal model. 
This is a problematic case unfortunately. Williams (2006) explains in a very nice way how one can fit this model by using many 
logistic regressions in an incremental way. The number of logistic regression models to be fit is the number of categories of the 
response variables - 1.  

It may be the case that the message says "problematic region". In this case the optimization was not succesful and perhaps the deviance is not 
at the global minimum. For example, with the addition of one extra variable the deviance might increase. I know, this type of ordinal regression is hard. In these difficult situations other packages return "error". Another difficult I have seen is when "NA" appear in the coefficients. In this case I do not consider these coefficients, not their corresponding variables in the calculation of the deviance. 
}

\value{
A list including:
\item{message}{
If you hit negative probabilities, the message "problematic region" will be printed. Otherwise, this is NULL.
}
\item{be}{
The regression coefficients.
}
\item{devi}{
The deviance of the model. 
}
}

\references{
Williams, R. (2006). Generalized ordered logit/partial proportional odds models for ordinal dependent variables. Stata Journal, 6(1), 58-82
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{pc.skel}, \link{ridge.reg} }
}

\examples{
y <- factor( rbinom(100, 3, 0.5) )
x <- matrix( rnorm(100 * 3), ncol = 3)
ordinal.reg(y ~ x, data = data.frame(x) )
ordinal.reg(y ~ 1, data = data.frame(x) )
}

\keyword{ Regression model }
\keyword{ ordinal regression }