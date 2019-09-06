\name{Total causal effect of a node on another node}
\alias{ida}
\title{
Total causal effect of a node on another node
}

\description{
Total causal effect of a node on another node.
}

\usage{
ida(x, y, G, dataset)
}

\arguments{
\item{x}{
A number between 1 and the number of variables (nodes) indicating the node whose total causal effect we want to estimate. This is the independent variable. See details for more on this.
}
\item{y}{
A number between 1 and the number of variables (nodes) indicating the node who is the dependent variable. The goal is to estimate the total 
causal effect of x on y.
}
\item{G}{
A square matrix representing a (partially) directed graph which contains 0s and 1s. If G[i, j] = 2 it 
means there is an arrow from node i to node j. If G[i, j] = 1, there is an undirected edge between nodes i and j and there is no edge between nodes i and j if G[i, j] = 0.  
}
\item{dataset}{
The dataset. This is a numerical matrix with data.
}
}

\details{
The total causal effect defined in Pearl's do-calculus is E(Y|do(X=z+1)) - E(Y|do(X=z)). As Pearl described it, he used linear regression, hence this funciton works for continuous data which are assumed to be Gaussian.

We estimate a set of possible total causal effects using linear regression. if y is a parent, or a descendant, of x in G, the estimated causal effect of x on y is zero. If y is not a parent of x, we take the regression coefficient of x in the regression lm(y ~ x + pa(x)), where pa(x) denotes the parents of x. This is repeated using all parents of x (including the empty set) and all possible parents values of x and their combinations.

One restriction to bear in mind. If a collider is created that combination of nodes is not used in the regression. 
}

\value{
A lis tincluding:
\item{tc}{
A matrix with 4 elements. The first columns one is the estimated beta coefficient, the second is its standard error, its t-value and the p-value for testing whether this is equal to 0.
}
\item{mess}{
If the x node has no parents a message about this appears. Otherwise this is NULL.
}
}

\references{
M.H. Maathuis, M. Kalisch and P. Buehlmann (2009). Estimating high-dimensional intervention effects from observational data. Annals of Statistics 37, 3133-3164.

Pearl (2005). Causality. Models, reasoning and inference. Cambridge University Press, New York.
}

\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
  %%  ~~further notes~~
    %}
\seealso{
\code{ \link{plotnetwork}, \link{nei}, \link{pc.or} 
}
}

\examples{
dataset <- rdag2(1000, p = 20, nei = 3)$x
mod <- pc.skel(dataset ,alpha = 0.01)
G <- pc.or(mod)$G
ida(10, 15, G, dataset)
}
