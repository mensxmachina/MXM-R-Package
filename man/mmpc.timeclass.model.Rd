\name{Regression models based on SES.timeclass and MMPC.timeclass outputs}
\alias{mmpc.timeclass.model}
\alias{ses.timeclass.model}
\title{
Regression model(s) obtained from SES.timeclass or MMPC.timeclass
}

\description{
One or more regression models obtained from SES.timeclass or MMPC.timeclass, are returned. 
}

\usage{
mmpc.timeclass.model(target, dataset, id, reps, wei = NULL, mmpctimeclass.Object)
ses.timeclass.model(target, dataset, id, reps, wei = NULL, sestimeclass.Object, 
nsignat = 1)
}

\arguments{
\item{target}{
The class variable. Provide a vector or a factor with discrete numbers indicating the class. Its length is equal to the number of rows of the dataset. 
}
\item{dataset}{
The dataset; provide a matrix. Currently, only continuous datasets are supported. The dataset
contains longitudinal data, where each column is a variable. The repeated measurements are the samples. 
}
\item{id}{
A numeric vector containing the subjects or groups. Its length is equal to the number of rows of the dataset.
}
\item{reps}{
A numeric vector containing the time points of the subjects. Its length is equal to the number of rows of the dataset.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{mmpctimeclass.Object}{
An object with the results of an MMPC.timeclass run.
}
\item{sestimeclass.Object}{
An object with the results of a SES.timeclass run.
}
\item{nsignat}{
How many signatures to use. If nsignat = 1 (default value) the first set of variables will be used for the model. If you want more, then specify the nubmer of signatures you want. If you want the models based on all signatures, specify "all". If you put a number which is higher than the number of signatures, all models will be returned.
}
}

\details{
This command is useful if you want to see all models and check for example their fitting ability.
}

\value{
A list including:
\item{mod}{
Depending on the number of signatures requested, one or models will be returned.
}
\item{signature}{
A matrix (or just one vector if one signature only) with the variables of each signature, along with the BIC of the corresponding regression model.
}
}

\references{
Tsagris, M., Lagani, V., & Tsamardinos, I. (2018). Feature selection for high-dimensional glmm data. BMC bioinformatics, 19(1), 17.

McCullagh P., and Nelder J.A. (1989). Generalized linear models. CRC press, USA, 2nd edition.
}

\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
  %%  ~~further notes~~
    %}
\seealso{
  \code{ \link{MMPC.timeclass} }
}

\examples{
## assume these are longitudinal data, each column is a variable (or feature)
dataset <- matrix( rnorm(400 * 50), ncol = 50 ) 
id <- rep(1:80, each = 5)  ## 80 subjects
reps <- rep( seq(4, 12, by = 2), 80)  ## 5 time points for each subject
## dataset contains are the regression coefficients of each subject's values on the 
## reps (which is assumed to be time in this example)
target <- rep(0:1, each = 200)
a <- MMPC.timeclass(target, reps, id, dataset)
mmpc.timeclass.model(target, dataset, id, reps, mmpctimeclass.Object = a)
}

\keyword{ Regression modelling }