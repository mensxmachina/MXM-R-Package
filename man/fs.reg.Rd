\name{Forward selection regression}
\alias{fs.reg}
\title{
Variable selection in regression models with forward selection
}

\description{
Variable selection in regression models with forward selection
}

\usage{
fs.reg(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, test = NULL, 
user_test = NULL, stopping = "BIC", tol = 2, ncores = 1) 
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). In either case, only two cases are avaialble, either all data are continuous, or categorical. 
}
\item{ini}{
If you have a set of variables already start with this one. Otherwise leave it NULL. 
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for asmmmbsing p-values significance. Default value is 0.05.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{test}{
The regression model to use. Available options are most of the tests for SES and MMPC. The ones NOT available are "gSquare", "censIndER", "testIndMVreg", "testIndClogit" and "testIndSpearman". 
If you have continuous predictor variables in matrix form, you can put "testIndFisher" and this is pretty fast. Instead of calcualting partial F-tests, which requires liear regression models to be fit, 
it calcualtes partial correlation coefficients and this is much more efficient. 
}
\item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, 
the "test" argument is ignored.
}
\item{stopping}{
The stopping rule. The BIC is always used for all methods. If you have linear regression though you can change this to "adjrsq" and in this case the adjusted R qaured is used.
}
\item{tol}{
The difference bewtween two successive values of the stopping rule. By default this is is set to 2. If for example, the BIC difference between two succesive models is less than 2, the process stops and the last variable, even though significant does not enter the model.
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
%\item{maxit}{
%The maximum number of iterations required by the MM estimators in the robust binary or poisson regrression. This argument is %taken into account only if the argument robust is set to TRUE.
%}
}

\details{
If the current 'test' argument is defined as NULL or "auto" and the user_test argument is NULL then the algorithm automatically selects the best test based on the type of the data. Particularly:
\itemize{
	\item if target is a factor, the multinomial or the binary logistic regression is used. If the target has two values only,   binary logistic regression will be used.
	\item if target is a ordered factor, the ordered logit regression is used. Hence, if you want to use multinomial or ordinal logistic regression, make sure your target is factor. 
	\item if target is a numerical vector and the dataset is not a matrix, but a data.frame linear regression is used. If however, the dataset is a matrix, the correlation based forward selection is used. That is, instead of partial F-tests, we do partial correlation tests.
	\item if target is discrete numerical (counts), the poisson regression conditional independence test is used. If there are only two values, the binary logistic regression is to be used.
	\item if target is a Surv object, the Survival conditional independence test is used.
}
}

\value{
In the case of test="testIndMMReg" and class(dataset) = matrix, just one matrix is returned with the index of the selected variable(s), the p-value, the test statistic and the BIC value of each model. For all other cases, the output of the algorithm is S3 object including:
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{mat}{
A matrix with the variables and their latest test statistics and p-values.
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
%Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal %relations. In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). 
%}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{\link{glm.fsreg}, \link{lm.fsreg}, \link{bic.fsreg}, \link{bic.glm.fsreg}, \link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)

#simulate a dataset with continuous data
dataset <- matrix( runif(500 * 20, 1, 100), ncol = 20 )

#define a simulated class variable 
target <- rt(500, 10)

a0 <- fs.reg(target, dataset, threshold = 0.05, stopping = "BIC", tol = 2) 

a1 <- fs.reg(target, dataset, threshold = 0.05, test = "testIndRQ", stopping = "BIC", 
tol = 2) 

require(survival, quietly = TRUE)
y <- survival::Surv(rexp(500), rep(1, 500) )
a2 <- fs.reg(y, dataset, threshold = 0.05, test = "censIndWR", stopping = "BIC", tol = 2) 
a3 <- MMPC(target, dataset)

target <- factor( rbinom(500, 1, 0.6) )
b2 <- fs.reg(target, dataset, threshold = 0.05, test = NULL, stopping = "BIC", tol = 2) 
}

\keyword{ Markov Blanket }
\keyword{ Variable Selection }
