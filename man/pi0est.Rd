\name{Estimation of the percentage of Null p-values}
\alias{pi0est}
\title{
Estimation of the percentage of Null p-values
}

\description{
Estimation of the percentage of Null p-values.
}

\usage{
pi0est(p, lambda = seq(0.05, 0.95, by = 0.01), dof = 3) 
}

\arguments{
\item{p}{
A vector of p-values.
}
\item{lambda}{
A vector of values of the tuning parameter lambda.
}
\item{dof}{
Number of degrees of freedom to use when estimating pi_0 with smoothing splines.
}
}

\details{
The estimated proporiton of null p-values is estimated the algorithm by Storey and Tibshirani (2003). 
}

\value{
The estimated proportion of non significant (null) p-values. In the paper Storey and Tibshirani mention that the estimate of pi0 is with lambda=1, but in their R code they use the highest value of lambda and thus we do the same here.  
}
\references{
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
  \code{\link{conf.edge.lower}, \link{bn.skel.utils}, \link{mmhc.skel}, \link{pc.skel}, \link{corfs.network}, \link{local.mmhc.skel}  }
}

\examples{
## simulate a dataset with continuous data
y <- rdag2(1000, p = 20, nei = 3)
ind <- sample(1:20, 20)
x <- y$x[, ind]
mod <- pc.skel( x, method = "comb.fast", alpha = 0.01 ) 
pval <- exp(mod$pvalue)
pval <- lower.tri(pval)
pi0est(pval)
}
