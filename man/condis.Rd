\name{Conditional independence tests counting the number of times a possible collider d-separates two nodes}
\alias{condis}
\title{
Many conditional independence tests counting the number of times a possible collider d-separates two nodes
}

\description{
Many conditional independence tests counting the number of times a possible collider d-separates two nodes
.
}

\usage{
condis(ind1, ind2, cs1, cs2, Var, dat, type = "pearson", rob = FALSE, max_k = 2, R = 1 )
}

\arguments{
\item{ind1}{
The index of the one variable to be considered. 
}
\item{ind2}{
The index of the other variable to be considered. 
}
\item{cs1}{
The index or indices of the conditioning set of variable(s). These are the neighbours of node ind1.
}
\item{cs2}{
The index or indices of the conditioning set of variable(s). These are the neighbours of node ind2.
}
\item{Var}{
The index of the possible collider. 
}
\item{dat}{
A numerical matrix or a data.frame with numerical, binary, nominal and ordinal variables only. 
}
\item{type}{
This is either "pearson", "spearman", "cat", "distcor", "ci.mm" or "ci.fast". 
}
\item{rob}{
In case you want robust estimation of the Pearson correlation prior to applying the conditional independence test. This is activated only when type = "pearson".
}
\item{max_k}{
The maximum number of conditioning variables to consider. It can be the case that each node ind1 and ind2 has 10 neighbours. We should try all possible combinations of the neighbours of ind1 and then of ind2. To reduce the computational cost we search the subsets with at most max_k variables. 
}
\item{R}{
This is used by most tests, except for type = "ci.mm" and type = "ci.fast". 
}
}

\details{
This is to be used in the conservative version of Rule 0 of the Pc algorithm. When one wants to know whether a variable is a possible collider, Ramsey, Spirtes and Zhang (2005) propose to perform the following action. For every unshilded triple (X, Y, Z) check all subsets of X's possible parents and of Z's poissible parents. 
a) If Y is NOT in any such set conditional on which, X and Z are independent orient X - Y - Z as X -> Y <- Z. 
b) If Y is in ALL sets conditional on which, X and Z are independent, leave X - Y - Z as it is, i.e. a non-collider.
c) Mark the triple X - Y - Z as "unfaitfull" otherwise.
This modification leads to the so called conservative PC (CPC) algorithm.

A few years later, Colombo and Maathuis (2014) suggested a modification of the previous action, called the majority rule. If Y is less than 50\% of thet sets that render X and Z independent orient X - Y - Z as X -> Y <- Z. If Y is found in exactly 50\% of sets that render X and Z independent, this triple is marked as "ambiuous". This modification leads the so called majority rule PC (MPC) algorithm.

This function we have implemented here, does exactly this. It applies tests to many subsets and returns a matrix with two columns. The first one contains 0 or 1 and the second is the p-value. A value of 0 indicates absenece of the possible collider from the set that produced that p-value, whereas a value of 1 indicates its presence in the set. 

This way, we can measure the proportion of times the possible collider Y was in a subset that rendered X and Z independent.
}

\value{
A matrix with two columns. The second one is the logarithm of the p-value. The first one contains 0s and 1s. The value of 0 means that the candidate collider was not in that set which produced the relevant p-value, whereas a value of 1 indicates that it was a member of that conditioning set. 
}

\references{
Ramsey, J., Zhang, J., Spirtes, P., 2006. Adjacency-faithfulness and conservative causal
inference. Proceedings of the 22nd Annual Conference on Uncertainty in Artificial Intelligence (UAI 2006).
https://arxiv.org/ftp/arxiv/papers/1206/1206.6843.pdf

Colombo, Diego, and Marloes H. Maathuis (2014). Order-independent constraint-based causal structure learning. The Journal of Machine Learning Research 15(1): 3741--3782.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{SES}, \link{MMPC}, \link{testIndLogistic} 
}
}

\examples{
x <- rdag2(1000, p = 10, nei = 5)
G <- x$G
dat <- x$x
cs1 <- which(G[6, ] > 0  |  G[, 6] > 0)
cs2 <- which(G[7, ] > 0  |  G[, 7] > 0)
cs1 <- setdiff( cs1, c(7, 3) )
cs2 <- setdiff( cs2, c(6, 3) ) 
condis(6, 7, cs1, cs2, 3, dat, type = "pearson", rob = FALSE, max_k = 3, R = 1 )
}

\keyword{ Area under the curve }
\keyword{ receiver operating curve }
