\name{Incremental BIC values and final regression model of the FBED algorithm}
\alias{fbedreg.bic}
\title{
Incremental BIC values and final regression model of the FBED algorithm
}

\description{
Incremental BIC values and final regression model of the FBED algorithm. 
}

\usage{
fbedreg.bic(target, dataset, wei = NULL, fbedreg.object, test = NULL, graph = TRUE)
}

\arguments{
\item{target}{
A numeric vector containing the values of the target variable. It can also discrete data, binary data (as factor), nominal or ordinal data (as factor). In contrast to SES, no position of the target variable in the dataset is accepted. The target must be a numerical vector.
}
\item{dataset}{
A numeric matrix or data.frame containing the variables. Rows are samples and columns are features. If you have categorical variables, this should be a data frame.
  }
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. We suggest not to use weights if you 
choose "testIndMMReg" as weights are alreayd being used there.
}
\item{fbedreg.object}{
An object with the results of an FBED run.
}
\item{test}{
If you know the test used in SES put it here, otherwise leave it NULL. It will take this information from the SES object. 
}
\item{graph}{
If you want a graphical representation of the drop in the BIC values set this to TRUE.
}
}

\details{
This function takes the output of the FBED (\code{\link{fbed.reg}}) and fits succesive models calculating the BIC for each of them. A graph can also be returned. 
}

\value{
A list including:
\item{res}{
A matrix with the selected variables, their test statistic and p-value (taken from the fbedreg.object) along with the BIC. Each row contains a BIC, that is the BIC of the model with the variables up to that row.
}
\item{mod}{
The final model with all selected variables.
}
}

\references{
Borboudakis G. and Tsamardinos I. (2019). Forward-backward selection with early dropping. Journal of Machine Learning Research, 20(8): 1-39.
}

\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
  %%  ~~further notes~~
    %}
\seealso{
  \code{ \link{reg.fit}, \link{mmpc.model}, \link{SES}, \link{MMPC}, \link{cv.ses}, \link{cv.mmpc} }
}

\examples{
dataset <- matrix( runif(400 * 50, 1, 100), ncol = 50 )
#define a simulated class variable 
target <- rt(400, 10)
a <- fbed.reg(target, dataset, K = 10, test = "testIndFisher", method = "eBIC") 
fbedreg.bic(target, dataset, fbedreg.object = a, test = "testIndFisher")
}

