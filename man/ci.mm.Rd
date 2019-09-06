\name{Symmetric conditional independence test with mixed data}
\alias{ci.mm}
\alias{ci.fast}
\alias{ci.mm2}
\alias{ci.fast2}
\title{
Symmetric conditional independence test with mixed data
}

\description{
Symmetric conditional independence test with mixed data.
}

\usage{
ci.mm(ind1, ind2, cs = NULL, dat, type, rob = FALSE, R = 1) 
ci.fast(ind1, ind2, cs = NULL, dat, type, rob = FALSE, R = 1)
ci.mm2(ind1, ind2, cs = NULL, suffStat) 
ci.fast2(ind1, ind2, cs = NULL, suffStat) 
}

\arguments{
\item{ind1}{
The index of the one variable to be considered. 
}
\item{ind2}{
The index of the other variable to be considered. 
}
\item{cs}{
The index or indices of the conditioning set of variable(s). If you have no variables set this equal to 0.
}
\item{dat}{
A data.frame with numerical, binary, nominal and ordinal variables only. 
}
\item{type}{
This is obsolete basically, but we need it here, so that the functions ci.mm and ci.fast have the same signatures as in \code{\link{cat.ci}}, \code{\link{condi}} and \code{\link{dist.condi}}.
}
\item{rob}{
This is obsolete basically, but we need it here, so that the functions ci.mm and ci.fast have the same signatures as in \code{\link{cat.ci}}, \code{\link{condi}} and \code{\link{dist.condi}}.
}
\item{R}{
This is obsolete basically, but we need it here, so that the functions ci.mm and ci.fast have the same signatures as in \code{\link{cat.ci}}, \code{\link{condi}} and \code{\link{dist.condi}}.
}
\item{suffStat}{
This is a list with only the dataset and the name must be "dataset". A data.frame with numerical, binary, nominal and ordinal variables only. 
}
}

\details{
The functions "ci.mm" and "ci.fast" are general functions to be used anywhere. The functions "ci.mm2" and "ci.fast2" are designed to be accepted by the command "pc" in the package "pcalg". The functions "ci.mm2" and "ci.fast2" can be fed in the "pc" function of the "pcalg" package in order to produce a PDAG with mixed data using the PC algorithm. For more information see the relevant paper in the references.

The "ci.mm" and "ci.fast" work with linear, logistic, multinomial and ordinal regression, whereas the "ci.mm2" and "ci.fast2" work with linear, logistic and ordinal regression only.
}

\value{
A vector including the test statistic, it's associated p-value and the relevant degrees of freedom. 
}

\references{
Tsagris M., Borboudakis G., Lagani V. and Tsamardinos I. (2018). Constraint-based Causal Discovery with Mixed Data. 
International Journal of Data Science and Analytics. 

https://link.springer.com/content/pdf/10.1007%2Fs41060-018-0097-y.pdf
 
http://nugget.unisa.edu.au/CD2017/papersonly/constraint-based-causal-r1.pdf

Spirtes P.,  Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, 3nd edition.

Sedgewick, A. J., Ramsey, J. D., Spirtes, P., Glymour, C., & Benos, P. V. (2017). Mixed Graphical Models for Causal Analysis of Multi-modal Variables. arXiv preprint arXiv:1704.02621.
}

\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
  %%  ~~further notes~~
    %}
\seealso{
\code{ \link{pc.skel}, \link{mmhc.skel}, \link{cat.ci}, \link{condi} }
}

\examples{
ci.mm(1, 2, dat = iris)
ci.mm(1, 5, dat = iris)
ci.fast(1, 5, dat = iris)
x <- iris
x[, 5] <- as.numeric(x[ ,5])  ## Caution:: this will be treated as ordered variable. 
x[, 5] <- factor(x[, 5], ordered = TRUE)
## ci.mm2 and ci.fast2 do not perform multinomial regression.
ci.mm2(1, 5, suffStat = list(dataset = x) )
ci.fast2(1, 5, suffStat = list(dataset = x) )
}
