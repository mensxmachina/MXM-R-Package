\name{Backward selection regression using the eBIC}
\alias{ebic.bsreg}
\title{
Backward selection regression using the eBIC
}

\description{
Backward selection regression using the eBIC
}

\usage{
ebic.bsreg(target, dataset, test = NULL, wei = NULL, gam = NULL) 
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. 
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). 
}
\item{test}{
The available tests: "testIndReg", "testIndPois", "testIndNB", "testIndLogistic", "testIndMMReg", 
"testIndBinom", "censIndCR", "censIndWR", "testIndBeta", "testIndZIP", "testIndGamma", "testIndNormLog" 
and "testIndTobit".
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when testIndMMReg is used.
}
\item{gam}{
In case the method is chosen to be "eBIC" one can also specify the \eqn{gamma} parameter. The default value is "NULL", so that the value is 
automatically calculated.
}
%\item{ncores}{
%How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
%}
}

\details{
The algorithm is a variation of the usual forward selection. At every step, the most significant variable enters the selected variables set. 
In addition, only the significant variables stay and are further examined. The non signifcant ones are dropped. This goes until no variable 
can enter the set. The user has the option to redo this step 1 or more times (the argument K). In the end, a backward selection is performed
to remove falsely selected variables.
}

\value{
A list including:
\item{runtime}{
The runtime required.
}
\item{info}{
A matrix with the number of variables and the number of tests performed (or models fitted) at each round (value of K). 
}
\item{mat}{
A matrix with the selected variables and their eBIC.
}
\item{back.rem}{
The variables removed in the backward phase.
}
\item{back.n.tests}{
The number of models fitted in the backward phase.
}
}

\references{
Borboudakis G. and Tsamardinos I. (2017). Forward-Backward Selection with Early Dropping. 
https://arxiv.org/pdf/1705.10770.pdf
}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{ \link{fs.reg}, \link{bic.fsreg}, \link{MMPC} }
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix( runif(100 * 20, 1, 100), ncol = 20 )
#define a simulated class variable 
target <- rt(100, 10)

a1 <- ebic.bsreg(target, dataset, test = "testIndReg") 
target <- rpois(100, 10)
a2 <- ebic.bsreg(target, dataset, test = "testIndPois") 
}

