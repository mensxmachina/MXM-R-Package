\name{The max-min Markov blanket algorithm}
\alias{mmmb}
\title{
Max-min Markov blanket algorithm
}

\description{
The MMMB algorithm follows a forward-backward filter approach for feature selection in order to provide a minimal, highly-predictive, feature subset of a high dimensional dataset. See also Details. 
}

\usage{
mmmb(target , dataset , max_k = 3 , threshold = 0.05 , test = "testIndFisher", 
user_test = NULL, ncores = 1)
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). In either case, only two cases are avaialble, either all data are continuous, or categorical. 

}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). Integer, default value is 3.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing the p-values. Default value is 0.05.
}
\item{test}{
The conditional independence test to use. Default value is "testIndFisher". See also \code{link{CondIndTests}}.
}
\item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, the "test" argument is ignored.
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cammmb it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
}

\details{
The idea is to run the MMPC algorithm at first and identify the parents and children (PC) of the target variable. As a second step, the MMPC algorithm is run on the discovered variables to return PC. The parents of the children of the target are the spouses of the target. Every variable in PCi is checked to see if it is a spouse of the target. If yes, it is included in the Markov Blanket of the target, otherwise it is thrown. If the data are continous, the Fisher correlation test is used or the Spearman correlation (more robust). If the data are categorical, the \eqn{G^2} test is used. 
}

\value{
The output of the algorithm is an S3 object including:
\item{mb}{
The Markov Blanket of the target variable. The parents and children of the target variable, along with the spouses of the target, which are the parents of the children of the target variable.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

\references{
Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal relations. In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining, pp. 673-678.

Tsagis M. (2018). Guide on performing feature selection with the R package MXM.
http://mensxmachina.org/site/wp-content/uploads/2018/04/Guide-on-performing-feature-selection-with-the-R-package-MXM.pdf
}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{\link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)
#simulate a dataset with continuous data
dataset <- matrix( runif(200 * 50, 1, 100), ncol = 50 )
#define a simulated class variable 
target <- 3 * dataset[, 10] + 2 * dataset[, 50] + 3 * dataset[, 20] + rnorm(200, 0, 5)
a1 <- mmmb(target , dataset , max_k = 3 , threshold = 0.05, test= "testIndFisher", 
ncores = 1,)
a2 <- MMPC(target, dataset, test="testIndFisher")
}

\keyword{ Markov Blanket }
\keyword{ Variable Selection }
