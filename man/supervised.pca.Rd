\name{Supervised PCA}
\alias{supervised.pca}
\title{
Supervised PCA
}

\description{
Supervised PCA.
}

\usage{
supervised.pca(target, dataset, indices, center = TRUE, scale = TRUE, 
colours = NULL, graph = TRUE) 
}

\arguments{
\item{target}{
A numerical vector or a factor denoting the class of each sample, the response variable.
}
\item{dataset}{
A matrix with numerical data (the predictor variables).
}
\item{indices}{
A vector with indices denoting whcih variables have been selected.
}
\item{center}{
In the calculation of the PCA, should the data be centered? Default value is TRUE.
}
\item{scale}{
In the calculation of the PCA, should the data be scaled to unity variance? Default value is TRUE.
}
\item{colours}{
Should the colour of the points be defined by the target variable or do you want to pass your own colours? This must be a vector whose length is equal to the length of the target.
}
\item{graph}{
Should two graphs be returned? The scores of the frist two principal components based on all the data and based on the selected variables.
}
}

\details{
This is not exactly the standard supervised PCA as suggested by Bair et al (2006). What we do here essentially is the following: PCA on all variables and on the variables selected by a variable selection algortihm. 
}

\value{
A list including:
\item{mod.all}{
The output returned by \code{\link{prcomp}} applied to all variables.
}
\item{mode.sel}{
The output returned by \code{\link{prcomp}} applied to the selected variables.
}
\item{var.percent}{
The percentage of variance explained by the selected variables.
}
}

\references{
Bair E., Hastie T., Debashis P. and Tibshirani R. (2006). Prediction by supervised principal components. Journal of the American Statistical Association 101(473): 119--137.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{gomp}, \link{fbed.reg}, \link{MMPC} 
}
}

\examples{
x <- as.matrix(iris[, 1:4])
target <- iris[, 5]
supervised.pca(target, x, indices = 1:2)
}

\keyword{ Area under the curve }
\keyword{ receiver operating curve }
