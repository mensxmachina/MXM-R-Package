\name{Bootstrap bias correction for the performance of the cross-validation procedure}
\alias{bbc}

\title{
Bootstrap bias correction for the performance of the cross-validation procedure
}

\description{
Bootstrap bias correction for the performance of the cross-validation procedure.
}

\usage{
bbc(predictions, target, metric = "auc.mxm", conf = 0.95, B = 1000)
}

\arguments{
\item{predictions}{
A matrix with the predictived values.
}
\item{target}{
A vector with the target variable, survival object, factor (ordered or unordered) or a numerical vector.
}
\item{metric}{
The possible values are: "auc.mxm", "acc.mxm", "ord_mae.mxm", "mae.mxm", "mse.mxm", "pve.mxm", "ci.mxm", "ciwr.mxm", "fscore.mxm", 
"euclid_sens.spec.mxm", "acc_multinom.mxm", "poisdev.mxm", "binomdev.mxm". The "nbdev.mxm" (negative binomial deviance) is missing. For more information on these see \code{\link{cv.ses}}. bold{Note} that they come with "".
}
\item{conf}{
A number between 0 and 1, the confidence level.
}
\item{B}{
The number of bootstrap replicates. The default number is 1000.
}
}

\details{
Upon completion of the cross-validation, the predicted values produced by all predictive models across all folds is collected in a matrix \eqn{P} of dimensions \eqn{n \times M}, where \eqn{n} is the number of samples and \eqn{M} the number of trained models or configurations. Sampled with replacement a fraction of rows (predictions) from \eqn{P} are denoted as the in-sample values. On average, the newly created set will be comprised by 63.2\% of the original individuals (The probability of sampling, with replacement, a sample of \eqn{n} numbers from a set of \eqn{n} numbers is \eqn{1-\left(1-\frac{1}{n} \right)^n \simeq 1-\frac{1}{e}=0.632}.), whereas the rest 36.8\% will be random copies of them. The non re-sampled rows are denoted as out-of-sample values. The performance of each model in the in-sample rows is calculated and the model (or configuration) with the optimal performance is selected, followed by the calculation of performance in the out-of-sample values. This process is repeated B times and the average performance is returned. 

Note, that the only computational overhead is with the repetitive re-sampling and calculation of the predictive performance, i.e. no model is fitted nor trained. The final estimated performance usually underestimates the true performance, but this negative bias is smaller than the optimistic uncorrected performance. 

Note, that all metrics are for maximization. For this reason "mse.mxm", "mae.mxm", "ord_mae.mxm", "poisdev.mxm", "binomdev.mxm" are multiplied by -1.
}

\value{
A list including:
\item{bbc.perf}{
The bootstrap bias corrected performance of the chosen algorithm, model or configuration.
}
\item{out.perf}{
The B out sampled performances. Their mean is the "bbc.perf" given above.
}
\item{ci}{
The (1- conf)\% confidence interval of the BBC performance. It is based on the empirical or percentile method for bootstrap samples. The lower and upper 2.5\% of the "out.perf".
}
}

\references{
Ioannis Tsamardinos, Elissavet Greasidou and Giorgos Borboudakis (2018). Bootstrapping the out-of-sample predictions for efficient and accurate cross-validation. Machine Learning (To appear).  

https://link.springer.com/article/10.1007/s10994-018-5714-4
}

\author{
  R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}

\seealso{
\code{ \link{cv.ses}, \link{cv.gomp} }
}

\examples{
predictions <- matrix(rbinom(200 * 100, 1, 0.7), ncol = 100) 
target <- rbinom(200, 1, 0.5)
bbc(predictions, target, metric = "auc.mxm")
}

\keyword{ Cross validation}
\keyword{ SES }
\keyword{ MMPC }
\keyword{ parallel }

