\name{Generate random folds for cross-validation}
\alias{generatefolds}
\title{
Generate random folds for cross-validation
}

\description{
Random folds for use in a cross validation are generated. There is the option for stratified splitting as well. 
}

\usage{
generatefolds(target, nfolds = 10, stratified = TRUE, seed = FALSE) 
}

\arguments{
\item{target}{
A vector with some data, either continuous or categorical. 
}
\item{nfolds}{
The number of folds to produce.
}
\item{stratified}{
A boolean variable specifying whether stratified random (TRUE) or simple random (FALSE) sampling is to be used when producing the folds.
}
\item{seed}{
A boolean variable. If set to TRUE, the folds will always be the same. 
}
}

\details{
I was inspired by the sam command in the package \bold{TunePareto} in order to do the stratified version. 
}


\value{
A list with nfolds elements where each elements is a fold containing the indices of the data.
}

%\references{
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{\link{cv.ses} }
}

\examples{
a <- generatefolds(iris[, 5], nfolds = 5, stratified = TRUE)
table(iris[a[[1]], 5])  ## 10 values from each group
}

\keyword{ Generation of folds }

