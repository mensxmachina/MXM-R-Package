\name{Regression models based on SES and MMPC outputs}
\alias{ses.model}
\alias{mmpc.model}
\alias{waldses.model}
\alias{waldmmpc.model}
\title{
Regression model(s) obtained from SES or MMPC
}

\description{
One or more regression models obtained from SES or MMPC, are returned. 
}

\usage{
ses.model(target, dataset, wei = NULL, sesObject, nsignat = 1, test = NULL)

mmpc.model(target, dataset, wei = NULL, mmpcObject, test = NULL)

waldses.model(target, dataset, wei = NULL, wald.sesObject, nsignat = 1, test = NULL)

waldmmpc.model(target, dataset, wei = NULL, wald.mmpcObject, test = NULL)
}

\arguments{
\item{target}{
A numeric vector containing the values of the target variable. If the values are proportions or percentages, 
i.e. strictly within 0 and 1 they are mapped into R using log( target/(1 - target) ). It can also discrete data, binary data (as factor), nominal or ordinal data (as factor). In contrast to SES, no position of the target variable in the dataset is accepted. The target must be a numerical vector.
}
\item{dataset}{
A numeric matrix or data.frame containing the variables. Rows are samples and columns are features. If you have categorical variables, this should be a data frame.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{sesObject}{
An object with the results of a SES run.
}
\item{mmpcObject}{
An object with the results of an MMPC run.
}
\item{wald.sesObject}{
An object with the results of a wald.ses run.
}
\item{wald.mmpcObject}{
An object with the results of an wald.mmpc run.
}
\item{nsignat}{
How many signatures to use. If nsignat = 1 (default value) the first set of variables will be used for the model. If you want more, then specify the nubmer of signatures you want. If you want the models based on all signatures, specify "all". If you put a number which is higher than the number of signatures, all models will be returned.
}
\item{test}{
If you know the test used in SES put it here, otherwise leave it NULL. It will take this information from the SEs object. If you used a robust version of a test (wherever possible), robust model(s) will be created.
}
}

\details{
This command is useful if you want to see all models and check for example their fitting ability, MSE in linear models for example.
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
Aitchison J. (1986). The Statistical Analysis of Compositional Data, Chapman & Hall; reprinted in 2003, with additional material, by The Blackburn Press.

Cox D.R. (1972). Regression models and life-tables. J. R. Stat. Soc., 34, 187-220.

Draper, N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition.

Ferrari S.L.P. and Cribari-Neto F. (2004). Beta Regression for Modelling Rates and Proportions. Journal of Applied Statistics, 31(7): 799-815.

Gutenbrunner C., Jureckova J., Koenker R. and Portnoy S. (1993). Tests of Linear Hypothesis based on Regression Rank Scores, Journal of NonParametric Statistics 2, 307-331.

Joseph M.H. (2011). Negative Binomial Regression. Cambridge University Press, 2nd edition.

Koenker R.W. (2005). Quantile Regression, Cambridge University Press.

Lagani V., Kortas G. and Tsamardinos I. (2013). Biomarker signature identification in "omics" with multiclass outcome. Computational and Structural Biotechnology Journal, 6(7): 1-7.

Lagani V. and Tsamardinos I. (2010). Structure-based variable selection for survival data. Bioinformatics Journal 16(15): 1887-1894.

Lambert, Diane (1992). Zero-inflated Poisson regression, with an application to defects in manufacturing. Technometrics 34(1)1: 1-14.

Mardia K.V., Kent J.T. and Bibby J.M. (1979). Multivariate Analysis, Academic Press, New York, USA. 

Maronna R.D. Yohai M.V. (2006). Robust Statistics, Theory and Methods. Wiley.

McCullagh P., and Nelder J.A. (1989). Generalized linear models. CRC press, USA, 2nd edition.
}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{SES}, \link{MMPC}, \link{cv.ses}, \link{cv.mmpc} }
}

\examples{
# simulate a dataset with continuous data
dataset <- matrix( runif(500 * 20, 1, 20), nrow = 500 ) 
#the target feature is the last column of the dataset as a vector
target <- dataset[, 20]
dataset <- dataset[, -20]
sesObject <- SES(target , dataset , max_k=3 , threshold = 0.05)
ses.model(target, dataset, sesObject = sesObject, nsignat = 1, test = NULL) 
ses.model(target, dataset, sesObject = sesObject, nsignat = 40, test = NULL)
mmpcObject <- MMPC(target, dataset, max_k=3, threshold = 0.05)
mmpc.model(target, dataset, mmpcObject = mmpcObject, test = NULL) 
}

\keyword{ Regression modelling }