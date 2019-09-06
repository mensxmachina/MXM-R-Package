\name{BIC based forward selection with generalised linear models}
\alias{bic.glm.fsreg}
\alias{bic.mm.fsreg}
\alias{bic.gammafsreg}
\alias{bic.normlog.fsreg}
\title{
Variable selection in generalised linear models with forward selection based on BIC

}

\description{
Variable selection in generalised linear models with forward selection based on BIC
}

\usage{
bic.glm.fsreg( target, dataset, wei = NULL, tol = 0, ncores = 1) 

bic.mm.fsreg( target, dataset, wei = NULL, tol = 0, ncores = 1) 

bic.gammafsreg(target, dataset, wei = NULL, tol = 0, ncores = 1) 

bic.normlog.fsreg(target, dataset, wei = NULL, tol = 0, ncores = 1) 
}

\arguments{
\item{target}{
The class variable. It can be either a vector with binary data (binomial regression), counts (poisson regression). If none of these is identified, linear regression is used.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). These can be continous and or categorical.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when robust is set to TRUE.
}
\item{tol}{
The difference bewtween two successive values of BIC. By default this is is set to 2. If for example, the BIC difference between two succesive models is less than 2, the process stops and the last variable, even though significant does not enter the model. 
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cammmb it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
}

\details{
Forward selection via the BIC is implemented. A variable which results in a reduction of BIC will be included, until the reduction is below a threshold set by the user (argument "tol").
}

\value{
The output of the algorithm is S3 object including:
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{mat}{
A matrix with the variables and their latest test statistics and p-values.
}
\item{info}{
A matrix with the selected variables, and the BIC of the model with that and all the previous variables.
}
\item{ci_test}{
The conditional independence test used. 
}
\item{final}{
The final regression model.
}
}

%\references{
%Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal %relations. In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). 
%}


\author{
Michail Tsagris

R implementation and documentation: Giorgos Aathineou <athineou@csd.uoc.gr> Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{\link{fs.reg}, \link{lm.fsreg}, \link{bic.fsreg}, \link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)
dataset <- matrix( runif(500 * 50, 1, 100), ncol = 50 )
target <- 3 * dataset[, 10] + 2 * dataset[, 20] + 3 * dataset[, 30] + rnorm(500, 0, 5)
a1 <- bic.glm.fsreg(target, dataset, tol = 2, ncores = 1 ) 
a2 <- bic.glm.fsreg( round(target), dataset, tol = 2, ncores = 1 ) 
y <- target   ;   me <- median(target)  ;   y[ y < me ] <- 0   ;   y[ y >= me ] <- 1
a3 <- bic.glm.fsreg( y, dataset, tol = 2, ncores = 1 ) 
}

\keyword{ Markov Blanket }
\keyword{ Variable Selection }
