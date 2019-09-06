\name{Regression models fitting}
\alias{reg.fit}
\title{
Regression modelling
}

\description{
Generic regression modelling function.
}

\usage{
reg.fit(y, dataset, event = NULL, reps = NULL, group = NULL, slopes = FALSE, 
reml = FALSE, model = NULL, wei = NULL, xnew = NULL)
}

\arguments{
\item{y}{
The target (dependent) variable. It can be a numerical variable, factor, ordinal factor, percentages, matrix, or time to event. 
If the values are proportions or percentages, i.e. strictly within 0 and 1 they are mapped into R using log( target/(1 - target) ).
If they are compositional data the additive log-ratio (multivariate logit) transformation is aplied beforehand. If the model is "clogit" 
this must be a matrix with two columns. The first column must be 0 and 1, standing for 0 = control and 1 = case. The second column is the id of the patients. A numerical variable, for example c(1,2,3,4,5,6,7,1,2,3,4,5,6,7).
}
\item{dataset}{
The indendent variable(s). It can be a vector, a matrix or a dataframe with continuous only variables, a data frame with mixed or only categorical variables.
}
\item{event}{
This is NULL unless you have time to event data (survival regression).
}
\item{reps}{
This is NULL unless you have time measurements (longitudinal data).
}
\item{group}{
This is NULL unless you have grouped (or clustered) data or longitudinal data (is the latter case the arugment reps is required). 
}
\item{slopes}{
This is for the longitudinal data only, TRUE or FALSE. Should random slopes be added or not?
}
\item{reml}{
This is for the longitudinal or grouped data only, TRUE or FALSE. If TRUE, REML will be used, otherwise ML will be used.
}
\item{model}{
The type of model you want to use. It can be specified by the user or left NULL, if other correct arguments have been passed. 
Poissible values (apart from NULL) are: "gaussian" (default), "binary", "binomial", "multinomial", "poisson", "ordinal", "gamma", "normlog", "tobit", "cox", "weibull", "exponential", "zip", "beta", "median", "negbin", "longitudinal", "grouped", "qpois" and "qbinom". The "zip" means 
that the zero part is constant, the variables are not associated with the excessive zeros. The value "grouped" refers to grouped data, but this does not have to be given if the argument "group" is given, but not the argument "reps. The "binomial" is when you have the number of successes 
and also the number of trials. "MM" stands for robust regression using MM estimation and "clogit" stands for conditional logistic regression.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL.
}
\item{xnew}{
If you have new data whose target values you want to predict put it here, otherwise leave it blank.
}
}

\details{
This is a generic regression function, which offers prediction as well. It is important that you pass the arguments with their names, for example if you have time to event data, write "event = ..." and not just put your event variable. This will avoid confusion. For the mixed models you need to specify the relevant arguments, "slopes", "reps", "reml" and "group" 
}

\value{
A list including:
\item{mod}{
The fitted model.
}
\item{pred}{
If you have new data the predicted values of the target (dependent) variable.
}
}

\references{
Almost the same as in \code{\link{CondIndTests}}.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{modeler}, \link{fbedreg.bic}, \link{mmpc.model}, \link{ridge.reg} }
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(runif(100 * 5, 1, 100), nrow = 100 ) 
#the target feature is the last column of the dataset as a vector
target <- dataset[, 1]
dataset <- dataset[, -1]
a <- reg.fit(target, dataset)
}

\keyword{ Regression model }