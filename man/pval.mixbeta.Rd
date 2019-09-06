\name{Fit a mixture of beta distributions in p-values}
\alias{pval.mixbeta}
\title{
Fit a mixture of beta distributions in p-values
}

\description{
Fit a mixture of beta distributions in p-values.
}

\usage{
pval.mixbeta(p)
}

\arguments{
\item{p}{
A vector of p-values.
}
}

\details{
The p-values are assumed to follow a mixture of two beta distributions. The null p-values follow Be(1, 1) and the non-null p-values follow Be(\eqn{\xi}, 1). 
In the first step, the proportion of true null values using Storey and Tibshirani (2003) is calculated and then MLE is adopted to obtain \eqn{\xi}. 
For more information on this see Triantafillou (2014).
}

\value{
A vector with the estimated \eqn{\pi_0} and \eqn{\xi} values. 
}

\references{
Triantafillou S., Tsamardinos I. and Roumpelaki A. (2014). Learning neighborhoods of high confidence in constraint-based causal discovery. 
In European Workshop on Probabilistic Graphical Models, pp. 487-502.

Storey J.D. and Tibshirani R. (2003). Statistical significance for genome-wide experiments. Proceedings of the National Academy of Sciences, 100: 9440-9445. 
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{pc.skel}, \link{mmhc.skel}, \link{corfs.network}, \link{local.mmhc.skel}, \link{conf.edge.lower}  }
}

\examples{
## simulate a dataset with continuous data
y <- rdag2(400, p = 25, nei = 3)
ind <- sample(1:25, 25)
x <- y$x[, ind]
mod <- pc.skel( x, method = "comb.fast", alpha = 0.01 ) 
pval <- as.vector( mod$pvalue[lower.tri(mod$pvalue)] )
pval.mixbeta(pval)
}
