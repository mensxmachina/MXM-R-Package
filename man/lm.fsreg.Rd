\name{Forward selection with linear regression models}
\alias{lm.fsreg}
\title{
Variable selection in linear regression models with forward selection

}

\description{
Variable selection in linear regression models with forward selection
}

\usage{
lm.fsreg(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, stopping = "BIC", 
tol = 2, ncores = 1)
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). In either case, only two cases are avaialble, either all data are continuous, or categorical. 
}
\item{ini}{
If you have a set of variables already start with this one. Currently this can only be a matrix with continuous variables. In such cases, the dataset must also contain continuous variables only. Otherwise leave it NULL. 
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing the p-values significance. Default value is 0.05.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when robust is set to TRUE.
}
\item{stopping}{
The stopping rule. The BIC ("BIC") or the adjusted \eqn{R^2} ("adjrsq") can be used. 
}
\item{tol}{
The difference bewtween two successive values of the stopping rule. By default this is is set to 2. If for example, the BIC difference between two succesive models is less than 2, the process stops and the last variable, even though significant does not enter the model. If the adjusted \eqn{R^2} is used, the tol should be something like 0.01 or 0.02. 
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cammmb it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
}

\details{
Only linear regression (robust and non robust) is supported from this function.  
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
A matrix with the selected variables, their logged p-values and test statistics. Each row corresponds to a model which contains the variables up to that line. The BIC in the last column is the BIC of that model.
}
\item{models}{
The regression models, one at each step.
}
\item{ci_test}{
The conditional independence test used, "testIndReg".  
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
\code{\link{fs.reg}, \link{lm.fsreg}, \link{bic.fsreg}, \link{bic.glm.fsreg}. \link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)

#simulate a dataset with continuous data
dataset <- matrix( runif(500 * 30, 1, 100), ncol = 30 )

#define a simulated class variable 
target <- 3 * dataset[, 10] + 2 * dataset[, 20] + rnorm(500, 0, 5)
a1 <- lm.fsreg(target, dataset, threshold = 0.05, stopping = "BIC", tol = 2) 
}

\keyword{ Markov Blanket }
\keyword{ Variable Selection }
