\name{Calculation of the constant and slope for each subject over time}
\alias{group.mvbetas}
\title{
Calculation of the constant and slope for each subject over time
}

\description{
Calculation of the constant and slope for each subject over time.
}

\usage{
group.mvbetas(x, id, reps)
}

\arguments{
\item{x}{
The dataset; provide a numerical a matrix. This should contain longitudinal data. Each variable is a column, and rows are longitudinal data.
}
\item{id}{
This is a numerical vector denoting the subjects. Its length must be equal to the number of rows of the x matrix.  
}
\item{reps}{
This is a numerical vector with the time points. Its length must be equal to the number of rows of the x matrix.  
}
}

\details{
This function is used internally in \code{\link{SES}} and \code{\link{MMPC}} and does calculations required bys the first step of the \bold{Static-Longitudinal} scenario of Tsagris, Lagani and Tsamardinos (2018). The measurements of each subject are regressed againt time. 
So, for each subject, we get the constant and interecept over time and this is repated for very feature.
}

\value{
A matrix. The first r ( = length( unique(id) ), the nubmer of subjects ) rows contain the constants and the other r rows contain the slopes.  
}

\references{
Tsagris M., Lagani V., & Tsamardinos I. (2018). Feature selection for high-dimensional glmm data. BMC bioinformatics, 19(1), 17.

McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}


\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{ \link{fbed.gee.reg}, \link{glmm.bsreg}, \link{MMPC.glmm} }
}

\examples{
## assume these are longitudinal data, each column is a variable (or feature)
x <- matrix( rnorm(100 * 100), ncol = 100 ) 
id <- rep(1:20, each = 5)  ## 20 subjects
reps <- rep( seq(4, 12, by = 2), 20)  ## 5 time points for each subject
a <- group.mvbetas(x, id, reps)
dim(a)  ## 5  200
## these are the regression coefficients of the first subject's values on the 
## reps (which is assumed to be time in this example)
a[c(1, 21), 1] 
coef( lm( x[id == 1, 1] ~ reps[1:5] ) )
}

