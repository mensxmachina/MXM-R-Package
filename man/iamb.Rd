\name{IAMB variable selection}
\alias{iamb}
\title{
IAMB variable selection
}

\description{
IAMB variable selection.
}

\usage{
iamb(target, dataset, threshold = 0.05, wei = NULL, test = NULL, user_test = NULL, 
stopping = "BIC", tol = 2, ncores = 1, back = "iambbs")
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. 
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = observations). In either case, only two cases are avaialble, either all data are continuous, or categorical. 
  }
\item{threshold}{
Threshold (suitable values in (0,1)) for assessing p-values significance. Default value is 0.05.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when robust is set to TRUE.
}
\item{test}{
The regression model to use. Available options are most of the tests for SES and MMPC. The ones NOT available are "gSquare", "censIndER", "testIndMVreg", "testIndClogit", "testIndSpearman" and "testIndFisher" and "testIndIGreg". 
}
\item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, the "test" argument is ignored.
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
\item{back}{
The backward phase. If this "iambbs" (default value) the IAMB backward phase is performed and hence the IAMB algorithm is completed. If "bs", a simple backward selection phase is performed. This way, the IAMB algorithm is slightly more general. 
}
}
\details{
IAMB stands for Incremental Association Markov Blanket. The algorithm comprises of a forward selection and a modified backward selection process. This functions does the modified 
backward selection process. In the usual backward selection, among the non singificant variabels, the one with the maximum p-value is dropped. 
So, one variable is removed at every step. In the IAMB backward phase, at aevery step, all non significant variables are removed. This makes it a lot faster. 
}

\value{
The output of the algorithm is a list of an S3 object including:
\item{vars}{
A vector with the selcted variables.
}
\item{mod}{
The output of the backward phase. In the case of no backward procedure this is the output of the forward phase. 
}
\item{mess}{
If the forward regression returned at most one variable, no backward procedure takes place and a message appears informing the user about this.
}
}

\references{
Tsamardinos, I., Aliferis, C. F., Statnikov, A. R., & Statnikov, E. (2003). Algorithms for Large Scale Markov Blanket Discovery. In FLAIRS conference, pp. 376-380.

Tsagis M. (2018). Guide on performing feature selection with the R package MXM.
http://mensxmachina.org/site/wp-content/uploads/2018/04/Guide-on-performing-feature-selection-with-the-R-package-MXM.pdf
}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
  \code{\link{glm.fsreg}, \link{lm.fsreg}, \link{bic.fsreg}, \link{bic.glm.fsreg}, \link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)
dataset <- matrix( runif(100 * 50, 1, 100), ncol = 50 )

target <- rpois(100, 10)
a1 <- iamb(target, dataset, threshold = 0.05, stopping = "BIC", tol = 0, back = "iambbs")
a2 <- iamb(target, dataset, threshold = 0.05, stopping = "BIC", tol = 0, back = "bs")
}

\keyword{ Backward regression }
\keyword{ Markov Blanket }
\keyword{ Variable Selection }
