\name{Certificate of exclusion from the selected variables set using SES or MMPC}
\alias{certificate.of.exclusion}
\alias{certificate.of.exclusion2}
\title{
Certificate of exclusion from the selected variables set using SES or MMPC
}

\description{
Information on why one ore more variables were not selected.
}

\usage{
certificate.of.exclusion(xIndex, sesObject = NULL, mmpcObject = NULL) 
certificate.of.exclusion2(xIndex, mmpc2object) 
}

\arguments{
\item{xIndex}{
A numerical vector with the  indices of the predictor variables.  
}
\item{sesObject}{
If you ran SES, wald.ses or perm.ses, give the whole SES object here, otherwise leave it NULL.
}
\item{mmpcObject}{
If you ran MMPC, wald.mmpc or prm.mmpc, give the whole MMPC object here, otherwise leave it NULL.
}
\item{mmpc2object}{
If you ran mmpc2, give the whole MMPC object here.
}
}

%\details{
%}

\value{
A list with the conditioning variables (if any), the test statistic and the logarithm of the p-value. In case a variable has been selected a message appears.
}

\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{\link{MMPC} }
}

\examples{
set.seed(123)
#simulate a dataset with continuous data
dataset <- matrix(runif(100 * 100, 1, 100), ncol = 100)
#define a simulated class variable 
target <- 3 * dataset[, 10] + 2 * dataset[, 100] + 3 * dataset[, 20] + rnorm(100, 0, 5)
# define some simulated equivalences
dataset[, 15] <- dataset[, 10] + rnorm(100, 0, 2)
dataset[, 100] <- dataset[, 100] + rnorm(100, 0, 2) 
dataset[, 20] <- dataset[, 100] + rnorm(100, 0, 2)
# run the SES algorithm
mod1 <- SES(target, dataset, max_k = 5, threshold = 0.05, test = "testIndFisher", 
hash = TRUE, hashObject = NULL);
mod2 <- MMPC(target, dataset, max_k = 5, threshold = 0.05, test = "testIndFisher", 
hash = TRUE, hashObject = NULL);
certificate.of.exclusion(c(10, 15, 30, 45, 20), mod1)
certificate.of.exclusion(c(10, 15, 30, 45, 20), NULL, mod2)
}

