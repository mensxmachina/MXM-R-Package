\name{Generalised orthogonal matching pursuit (gOMP)}
\alias{gomp}
\alias{gomp.path}
\alias{boot.gomp}
\title{
Generalised orthogonal matching pursuit (gOMP)
}

\description{
Generalised orthogonal matching pursuit.
}

\usage{
gomp(target, dataset, xstand = TRUE, tol = qchisq(0.95, 1), 
test = "testIndLogistic", method = "ar2" ) 

gomp.path(target, dataset, xstand = TRUE, tol = c(4, 5, 6), 
test = "testIndLogistic", method = "ar2" ) 

boot.gomp(target, dataset, tol = qchisq(0.95, 1), 
test = "testIndLogistic", method = "ar2", B = 500, ncores = 1) 
}

\arguments{
\item{target}{
The response variable, a numeric vector, a matrix or a Surv object.
}
\item{dataset}{
A matrix with continuous data, where the rows denote the samples and the columns are the variables. 
}
\item{xstand}{
If this is TRUE the independent variables are standardised. 
}
\item{tol}{
The tolerance value to terminate the algorithm. This is the change in the criterion value 
between two successive steps. The default value is the 95\% quantile of the \eqn{\chi^2} distribution
with 1 degree of freedom. For test = "testIndFisher" the BIC is already calculated. 

In the case of "gomp.path" this is a vector of values. For each tolerance value the result of the gOMP is returned. 
It returns the whole path of solutions. 
}
\item{test}{
This denotes the parametric model to be used each time. It depends upon the nature of the target variable. 
The possible values are "testIndFisher" (or "testIndReg" for the same purpose), "testIndLogistic", 
"testIndPois", "testIndQPois", "testIndQbinom", "testIndNormLog", "testIndMVreg", "testIndNB", "testIndBeta", 
"testIndGamma", "testIndMMReg", "testIndRQ", "testIndOrdinal", "testIndTobit", "censIndCR", "censIndWR", "censIndLLR" and 
"testIndMultinom".
}
\item{method}{
This is only for the "testIndFisher". You can either specify, "ar2" for the adjusted R-square or "sse" for the sum of squares 
of errors. The tolerance value in both cases must a number between 0 and 1. That will denote a percentage. If the percentage 
increase or decrease is less than the nubmer the algorithm stops. An alternative is "BIC" for BIC and the tolerance values 
are like in all other regression models.
}
\item{B}{
How many bootstrap samples to generate. The gOMP will be performed for each of these samples. 
}
\item{ncores}{
How many cores to use. This argument is valid only if you have a multi-threaded machine.
}
}

\value{
A list including:
\item{runtime}{
The runtime of the algorithm
}
\item{phi}{
The \eqn{phi} coefficient, returned in the quasi binomial (testIndQBinom), quasi Poisson (testIndQPois), Gamma (testIndGamma) and 
Gaussian with log link (testIndNormLog). In all other cases this is NULL. 
}
\item{res}{
For the case of "gomp" a matrix with two columns. The selected variable(s) and the criterion value at every step. 
For the case of "gomp.path" a matrix with many columns. Every column contains the selected variables for each tolerance calue, 
starting from the smallest value (which selected most variables). The final column is the deviance of the model at each step. 
For the "boot.gomp" this is a matrix with two columns. The first one is the selected variables and the second column is their
proportion of selection.
}
}

\references{
Pati Y. C., Rezaiifar R. & Krishnaprasad P. S. (1993). Orthogonal matching pursuit: Recursive function approximation with applications to wavelet decomposition. In Signals, Systems and Computers. 1993 Conference Record of The Twenty-Seventh Asilomar Conference on. IEEE.

Davis G. (1994). Adaptive Nonlinear Approximations. PhD thesis. 
http://www.geoffdavis.net/papers/dissertation.pdf

Mallat S. G. & Zhang Z. (1993). Matching pursuits with time-frequency dictionaries. IEEE Transactions on signal processing, 41(12), 3397-3415.
https://www.di.ens.fr/~mallat/papiers/MallatPursuit93.pdf

Gharavi-Alkhansari M., & Huang T. S. (1998, May). A fast orthogonal matching pursuit algorithm. In Acoustics, Speech and Signal Processing, 1998. 
Proceedings of the 1998 IEEE International Conference on (Vol. 3, pp. 1389-1392). IEEE.

Chen S., Billings S. A., & Luo W. (1989). Orthogonal least squares methods and their application to non-linear system identification. 
International Journal of control, 50(5), 1873-1896.
http://eprints.whiterose.ac.uk/78100/1/acse%20report%20343...pdf

Lozano A., Swirszcz G., & Abe N. (2011). Group orthogonal matching pursuit for logistic regression. In Proceedings of the Fourteenth International Conference on Artificial Intelligence and Statistics.

Razavi S. A. Ollila E., & Koivunen V. (2012). Robust greedy algorithms for compressed sensing. In Signal Processing 
Conference (EUSIPCO), 2012 Proceedings of the 20th European. IEEE.

Mazin Abdulrasool Hameed (2012). Comparative analysis of orthogonal matching pursuit and least angle regression. MSc thesis, Michigan State University.
https://www.google.gr/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0ahUKEwik9P3Yto7XAhUiCZoKHQ8XDr8QFgglMAA&url=https%3A%2F%2Fd.lib.msu.edu%2Fetd%2F1711%2Fdatastream%2FOBJ%2Fdownload%2FComparative_analysis_of_orthogonal_matching_pursuit_and_least_angle_regression.pdf&usg=AOvVaw2fRcSemcbteyWUiAKYi-8B

Tsagis M. (2018). Guide on performing feature selection with the R package MXM.
http://mensxmachina.org/site/wp-content/uploads/2018/04/Guide-on-performing-feature-selection-with-the-R-package-MXM.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{cor.fbed}, \link{cor.fsreg}, \link{correls}, \link{fs.reg}
}
}

\examples{
x <- matrix( rnorm(500 * 50), ncol = 50)
y <- rnorm(500)
b <- MXM::gomp(y, x, test = "testIndFisher")
}

