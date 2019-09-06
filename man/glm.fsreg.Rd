\name{Forward selection with generalised linear regression models}
\alias{glm.fsreg}
\alias{gammafsreg}
\alias{normlog.fsreg}
\title{
Variable selection in generalised linear regression models with forward selection
}

\description{
Variable selection in generalised linear regression models with forward selection
}

\usage{
glm.fsreg(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, tol = 2, 
ncores = 1) 

gammafsreg(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, tol = 2, 
ncores = 1) 

normlog.fsreg(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, tol = 2, 
ncores = 1) 
}

\arguments{
\item{target}{
A numerical vector or a factor variable with two levels. The Gamma regression reuqires strictly positive numbers, wehreas the normlog requires positive numbes, zero included.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). In either case, only two cases are available, either all data are continuous, or categorical. 
}
\item{ini}{
If you have a set of variables already start with this one. Otherwise leave it NULL. 
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing the p-values significance. Default value is 0.05.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when robust is set to TRUE.
}
\item{tol}{
The difference bewtween two successive values of the stopping rule. By default this is is set to 2. If for example, the BIC difference between two succesive models is less than 2, the process stops and the last variable, even though significant does not enter the model.
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cammmb it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
%\item{maxit}{
%The maximum number of iterations required by the MM estimators in the robust binary or poisson regrression. This argument is %taken into account only if the argument robust is set to TRUE.
%}
}

%\details{
%
%}

\value{
The output of the algorithm is S3 object including:
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{mat}{
A matrix with the variables and their latest test statistics and p-values. If you have logistic or Poisson regression with continuous predictor variables in a matrix form and no weights, this will not appear. In this case, a C++ code is called and the output is less.
}
\item{info}{
A matrix with the selected variables, their logged p-values and test statistics. Each row corresponds to a model which contains the variables up to that line. The BIC in the last column is the BIC of that model.
}
\item{ci_test}{
The conditional independence test used. 
}
\item{final}{
The final regression model.
}
}

%\references{
%Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal relations. In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). 
%}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{\link{fs.reg}, \link{lm.fsreg}, \link{bic.fsreg}, \link{bic.glm.fsreg}, \link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)

#simulate a dataset with continuous data
dataset <- matrix( runif(500 * 40, 1, 100), ncol = 40 )

#define a simulated class variable 
target <- rpois(500, 10)

a <- glm.fsreg(target, dataset, threshold = 0.05, tol = 2, ncores = 1 ) 
b <- MMPC(target, dataset, max_k = 3, threshold = 0.05, test = "testIndPois")
}

\keyword{ Markov Blanket }
\keyword{ Variable Selection }
