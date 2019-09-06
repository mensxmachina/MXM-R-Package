\name{Utilities for the skeleton of a (Bayesian) Network}
\alias{bn.skel.utils}
\alias{bn.skel.utils2}
\title{
Utilities for the skeleton of a (Bayesian) Network
}

\description{
Utilities for the skeleton of a (Bayesian) Network
}

\usage{
bn.skel.utils(mod, G = NULL, roc = TRUE, alpha = 0.01) 
bn.skel.utils2(mod, G = NULL, roc = TRUE, alpha = 0.01) 
}

\arguments{
\item{mod}{
An object as retured by pc.skel, glmm.pc.skel or mmhc.skel. 
}
\item{G}{
The true adjacency matrix with 1 indicating an edge and zero its absence. Symmetric or not is not important. If this is not available, 
leave it NULL.
}
\item{roc}{
Do you want a graph with the ROC curve be returned? Default value is TRUE.
}
\item{alpha}{
The significance level ( suitable values in (0, 1) ) for assessing the p-values. Default value is 0.01.
}
}

\details{
Given the true adjaceny matrix one can evaluate the estimated adjacency matrix, skeleton, of the Pc or the MMHC algorithm. 

The bn.skels.utils give you the area under the curve, false discovery rate and sorting of the edges based on their p-values. 

The bn.skel.utils2 estimates the confidence of each edge. The estimated proporiton of null p-values is estimated the algorithm by 
Storey and Tibshirani (2003). 
}

\value{
For the "bn.skel.utils" a list including:
\item{fdr}{
The false discovery rate as estimated using the Benjamini-Hochberg correction. 
}
\item{area}{
This is a list with the elements of the \code{\link{auc}} function. The area under the curve, the sensitivy and specificity for a range of values, the Youden index, etc.
}
\item{sig.pvalues}{
A matrix with the row and column of each significant p-value sorted in asending order. As we move down the matrix, the p-values increase and hence the strength of the associations decreases.
}
For the "bn.skel.utils2" a list including:
\item{area}{
This is a list with the elements of the \code{\link{auc}} function. The area under the curve, the sensitivy and specificity for a range of values, the Youden index, etc.
}
\item{pxy}{
A matrix with the row and column of the confidence of each p-value sorted in asending order. As we move down the matrix, the confidences decrease.
}
\item{lower}{
The lower confidcence limit of an edge as estimated by \code{\link{conf.edge.lower}}. 
}
}
\references{
Tsamardinos I. and Brown L.E. Bounding the False Discovery Rate in Local Bayesian Network Learning. AAAI, 2008.

Triantafillou S., Tsamardinos I. and Roumpelaki A. (2014). Learning neighborhoods of high confidence in constraint-based causal discovery. In European Workshop on Probabilistic Graphical Models, pp. 487-502.

Storey J.D. and Tibshirani R. (2003). Statistical significance for genome-wide experiments. Proceedings of the National Academy of Sciences, 100: 9440-9445. 

Benjamini Y. and Hochberg Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. 
Journal of the Royal Statistical Society Series B, 57(1), 289-300.

Spirtes P.,  Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, 3nd edition.
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
y <- rdag2(500, p = 25, nei = 3)
ind <- sample(1:25, 25)
x <- y$x[, ind]
mod <- pc.skel( x, method = "comb.fast", alpha = 0.01 ) 
G <- y$G[ind, ind]
G <- G + t(G)
bn.skel.utils(mod, G, roc = FALSE, alpha = 0.01) 
bn.skel.utils2(mod, G, roc = FALSE, alpha = 0.01) 
}
