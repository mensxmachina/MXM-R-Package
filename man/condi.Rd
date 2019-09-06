\name{Conditional independence tests with and without permutation p-value}
\alias{condi}
\alias{dist.condi}
\alias{cat.ci}
\title{
Conditional independence test for continuous class variables with and without permutation based p-value
}

\description{
The main task of this test is to provide a permutation based p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a 
conditioning set CS.
}

\usage{
condi(ind1, ind2, cs, dat, type = "pearson", rob = FALSE, R = 1) 
dist.condi(ind1, ind2, cs, dat, type = NULL, rob = FALSE, R = 499) 
cat.ci(ind1, ind2, cs, dat, type, rob = FALSE, R = 1)
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
A matrix with the data. In the case of "cat.ci" the minimum must be 0, i.e. the data must be like 0, 1, 2 and NOT 1, 2, 3... There is a C++ code behind and the minimum must be 0. 
}
\item{type}{
Do you want the Pearson (type = "pearson") or the Spearman (type = "spearman") correlation to be used. For "dist.condi" this is an obsolete argument but it requires to exist when it is used in 
the PC algorithm. 

For "cat.ci" this should be a vector with the levels, the number of distinct (different) values of each categorical variable. Its length is equal to the number of variables used in the test 
(2 + the number of conditioning variables). 
}
\item{rob}{
If you choose type="pearson" then you can sapecify whether you want a robust version of it. For "dist.condi" and "cat.ci" this is an obsolete argument but it requires to exist when it is used in the PC algorithm.
}
\item{R}{
If R = 1 then the asymptotic p-value is calculated. If R > 1 a permutation based p-value is returned. For the distance correlation based test, this is set to 499 by default and is used in the partial correlaiton test only. 
}
}

\details{
This test is currently designed for usage by the PC algorithm. The Fisher conditional independence test which is based on the Pearson or Spearman correlation coefficients is much faster than the distance based (partial) correlation test. 

The distance correlation can handle non linear relationships as well. The p-value for the partial distance correlation is calculated via permutations and is slow. 
}

\value{
A vector including the test statistic, it's associated p-value and the relevant degrees of freedom. In the case of a permutation based p-value, the returned test statistic is the observed test statistic divided by the relevant degrees of freedom (Pearson and Spearman correlation coefficients only). This is for the case of ties between many permutation based p-values. The PC algorithm choose a pair of variables based on the p-values. If they are equal it will use the test statistic. 
}

\references{
Hampel F. R., Ronchetti E. M., Rousseeuw P. J., and Stahel W. A. (1986). Robust statistics: the approach based on influence functions. John Wiley & Sons. 

Lee Rodgers J., and Nicewander W.A. (1988). "Thirteen ways to look at the correlation coefficient". The American Statistician 42(1): 59-66.

Shevlyakov G. and Smirnov P. (2011). Robust Estimation of the Correlation Coefficient: An Attempt of Survey. Austrian Journal of Statistics, 40(1 & 2): 147-156.

Spirtes P., Glymour C. and Scheines R. Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, second edition, January 2001.

Szekely G.J. and Rizzo, M.L. (2014). Partial distance correlation with methods for dissimilarities. The Annals of Statistics, 42(6): 2382--2412.

Szekely G.J. and Rizzo M.L. (2013).  Energy statistics: A class of statistics based on distances. Journal of Statistical Planning and Inference 143(8): 1249--1272.
}

\author{
Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{testIndFisher}, \link{testIndSpearman}, \link{pc.skel}, \link{gSquare}, \link{CondIndTests}
}
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(runif(500 * 5, 1, 100), ncol = 5 )
testIndFisher(dataset[, 1], dataset[, -1], xIndex = 1, csIndex = 2)
condi(ind1 = 1, ind2 = 2, cs = 3, dataset, R = 1)
condi(ind1 = 1, ind2 = 2, cs = 3, dataset, R = 999)
dist.condi(ind1 = 1, ind2 = 2, 0, dataset)
dist.condi(ind1 = 1, ind2 = 2, cs = 3, dataset, R = 99)
}

\keyword{ Fisher's test}
\keyword{ conditional independence test }
\keyword{ permutation based p-value }
