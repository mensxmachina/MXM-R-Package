\name{Feature selection using SES and MMPC for classifiication with longitudinal data}
\alias{SES.timeclass}
\alias{MMPC.timeclass}
\title{
Feature selection using SES and MMPC for classifiication with longitudinal data
}

\description{
SES algorithm follows a forward-backward filter approach for feature selection in order to provide minimal, highly-predictive, statistically-equivalent, multiple feature subsets of a high dimensional dataset. See also Details. MMPC algorithm follows the same approach without generating multiple feature subsets.
}

\usage{
SES.timeclass(target, reps, id, dataset, max_k = 3, threshold = 0.05, 
ini = NULL, wei = NULL, hash = FALSE, hashObject = NULL, ncores = 1) 

MMPC.timeclass(target, reps, id, dataset, max_k = 3, threshold = 0.05, 
ini = NULL, wei = NULL, hash = FALSE, hashObject = NULL, ncores = 1) 
}

\arguments{
\item{target}{
The class variable. Provide a vector or a factor with discrete numbers indicating the class. Its length is equal to the number of rows of the dataset. 
}
\item{reps}{
A numeric vector containing the time points of the subjects. Its length is equal to the number of rows of the dataset.
}
\item{id}{
A numeric vector containing the subjects or groups. Its length is equal to the number of rows of the dataset.
}
\item{dataset}{
The dataset; provide a matrix. Currently, only continuous datasets are supported. The dataset
contains longitudinal data, where each column is a variable. The repeated measurements are the samples. 
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). Integer, default value is 3.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{ini}{
This is a supposed to be a list. After running SES or MMPC with some hyper-parameters you might want to run SES again with different hyper-parameters. To avoid calculating the univariate associations (first step of SES and of MPPC) again, you can extract them from the first run of SES and plug them here. This can speed up the second run (and subequent runs of course) by 50\%. See the details and the argument "univ" in the output values. 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to store the statistics calculated during SES execution in a hash-type object. Default value is FALSE. If TRUE a hashObject is produced.
}
\item{hashObject}{
A List with the hash objects generated in a previous run of SES.glmm. 
Each time SES runs with "hash=TRUE" it produces a list of hashObjects that can be re-used in order to speed up next runs of SES.

Important: the generated hashObjects should be used only when the same dataset is re-analyzed, possibly with different values of max_k and threshold.
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is definetely not linear in the number of cores.
}
}

\details{
This is SES and MMPC used in the static-longitudinal scenario of Tsagris, Lagani and Tsamardinos (2018). The idea is that you have many features of longitudinal data for many subjects. For each subject you have calculated the coefficients of a simple linear regression over time and this is repeated for each feature. In the end, assuming p features, you have p constants and p slopes for each subject, each constant and slope refers to a feature for a subject. Hence, each new feature consists of two vectors, the constants and the slopes and the feature selection takes place there.
}

\value{
The output of the algorithm is an object of the class 'SESoutput' for SES or 'MMPCoutput' for MMPC including:
\item{selectedVars}{
The selected variables, i.e., the signature of the target variable.
}
\item{selectedVarsOrder}{
The order of the selected variables according to increasing pvalues.
}
\item{queues}{
A list containing a list (queue) of equivalent features for each variable included in selectedVars. An equivalent signature can be built by selecting a single feature from each queue. 
Featured only in SES.
}
\item{signatures}{
A matrix reporting all equivalent signatures (one signature for each row). Featured only in SES.
}
\item{hashObject}{
The hashObject caching the statistic calculated in the current run.
}
\item{pvalues}{
For each feature included in the dataset, this vector reports the strength of its association with the target in the context of all other variables. 
Particularly, this vector reports the max p-values found when the association of each variable with the target is tested against different conditional sets. 
Lower values indicate higher association. Note that these are the logged p-values, natural logarithm of the pvalues, and not the p-values.
}
\item{stats}{
The statistics corresponding to "pvalues" (higher values indicates higher association).
}
\item{univ}{
This is a list with the univariate associations. The test statistics and their corresponding logged p-values, along with their flag (1 if the test was perfromed and 0 otherwise). 
This list is very important for subsequent runs of SES with different hyper-parameters. After running SES with some hyper-parameters you might want to run SES again with different 
hyper-parameters. To avoid calculating the univariate associations (first step of SES or MMPC) again, you can take this list from the first run of SES and plug it in the argument "ini" 
in the next run(s) of SES or MMPC. This can speed up the second run (and subequent runs of course) by 50\%. See the argument "univ" in the output values. 
}
\item{max_k}{
The max_k option used in the current run.
}
\item{threshold}{
The threshold option used in the current run.
}
\item{n.tests}{
If you have set hash = TRUE, then the number of tests performed by SES or MMPC will be returned. If you have not set this to TRUE, the number of univariate associations will be returned. 
So be careful with this number.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{test}{
The character name of the statistic test used.
}

Generic Functions implemented for SESoutput Object:

\item{plot(object=SESoutput, mode="all")}{
Plots the generated pvalues (using barplot) of the current SESoutput object in comparison to the threshold.

Argument mode can be either "all" or "partial" for the first 500 pvalues of the object.
}
}

\references{
Tsagris M., Lagani V., & Tsamardinos I. (2018). Feature selection for high-dimensional glmm data. BMC bioinformatics, 19(1), 17.

Vincenzo Lagani, George Kortas and Ioannis Tsamardinos (2013), Biomarker signature identification in "omics" with multiclass outcome. Computational and Structural Biotechnology Journal, 6(7):1-7.

McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Ioannis Tsamardinos, Vincenzo Lagani

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Vincenzo Lagani <vlagani@csd.uoc.gr>  
}

\seealso{
\code{\link{mmpc.timeclass.model} }
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
}

\keyword{ SES }
\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }
\keyword{ Variable Selection }
