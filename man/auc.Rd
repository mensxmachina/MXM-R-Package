\name{ROC and AUC}
\alias{auc}
\title{
ROC and AUC
}

\description{
Receiver operating curve and area under the curve.
}

\usage{
auc(group, preds, roc = FALSE, cutoffs = NULL) 
}

\arguments{
\item{group}{
A numerical vector with the predicted values of each group as 0 and 1.
}
\item{preds}{
The predicted values of each group. 
}
\item{roc}{
If you want the ROC to appear set it to TRUE.
}
\item{cutoffs}{
If you provide a vector with decreasing numbers from 1 to 0 that will be used for the ROC, otherwise, the values from 1 to 0 with a step equal to -0.01 will be used.
}
}

\details{
The ara under the curve is returned. The user has the option of getting the receiver operating curve as well.
}

\value{
A list including:
\item{cutoffs}{
The cutoff values.
}
\item{sensitivity}{
The sensitivity values for each cutoff value.
}
\item{specificity}{
The specificity value for each cutoff value.
}
\item{youden}{
The pair of of 1- specificity and sensitivity where the Youden's J appears on the graph and the Youden index which is defined as the maximum value of sensitivity - specificity + 1. 
}
\item{auc}{
The area under the curve, plus a circle with the point where Youden's J is located. If "roc" is set to FALSE, this is the only item in the list to be returned.
}
}

%\references{
%
%}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{bbc}, \link{testIndLogistic} 
}
}

\examples{
g <- rbinom(150, 1, 0.6)
f <- rnorm(150)
auc(g, f, roc = FALSE)
}

\keyword{ Area under the curve }
\keyword{ receiver operating curve }
