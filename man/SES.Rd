\name{SES}
\alias{SES}
\title{
Feature selection algorithm for identifying multiple minimal, statistically-equivalent and equally-predictive feature signatures.
}
\description{
SES algorithm follows a forward-backward filter approach for feature selection in order to provide minimal, highly-predictive, statistically-equivalent, multiple feature subsets of a high dimensional dataset. See also Details.
}
\usage{
SES(target, dataset, max_k = 3, threshold = 0.05, test = NULL,
 user_test = NULL, hash = FALSE, hashObject = NULL, robust = FALSE)
}
\arguments{
  \item{target}{
The class variable. Provide either a string, an integer value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
  \item{dataset}{
The data-set; provide either a data frame or a matrix (columns = variables , rows = samples).
Alternatively, provide an ExpressionSet (in which case rows are samples and columns are features, see bioconductor for details).
}
  \item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). Integer, default value is 3.
}
  \item{threshold}{
Threshold (suitable values in [0,1]) for assessing p-values significance. Default value is 0.05.
}
  \item{test}{
The conditional independence test to use. Default value is NULL.

Available conditional independence tests:
\itemize{
  \item "testIndFisher": Fisher conditional independence test for continous targets.
  \item "testIndLogistic": Conditional Independence Test based on logistic regression for binary, categorical and ordinal targets.
  \item "gSquare": Conditional Independence test based on the G test of independence (log likelihood ratio test).
  \item "censIndLR": Conditional independence test for survival data based on the Log likelihood ratio test.
}
See also Details.
}
  \item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, the "test" argument is ignored.
}
  \item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to store the statistics calculated during SES execution in a hash-type object. Default value is FALSE. If TRUE a hashObject is produced.
}
  \item{hashObject}{
A List with the hash objects generated in a previous run of SES. 
Each time SES runs with "hash=TRUE" it produces a list of hashObjects that can be re-used in order to speed up next runs of SES.

Important: the generated hashObjects should be used only when the same dataset is re-analyzed, possibly with different values of max_k and threshold.
}

\item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. It takes more time than a non robust version but it is suggested in case of outliers. Default value is FALSE.
}
}
\details{
This function implements the Statistically Equivalent Signature (SES) algorithm as presented in "Tsamardinos, Lagani and Pappas, HSCBB 2012" 

(http://www.mensxmachina.org/publications/discovering-multiple-equivalent-biomarker-signatures/)

For faster computations in the internal SES functions, install the suggested package "\bold{gRbase}".

The max_k option: the maximum size of the conditioning set to use in the conditioning independence test. Larger values provide more accurate results, at the cost of higher computational times. When the sample size is small (e.g., <50 samples) the max_k parameter should be <=5, otherwise the conditional independence test may not be able to provide reliable results.

If the dataset contains missing (NA) values, they will automatically be replaced by the current variable (column) mean value with an appropriate warning to the user after the execution.

If the target is a single integer value or a string, it has to corresponds to the column number or to the name of the target feature in the dataset. In any other case the target is a variable that is not contained in the dataset.

If the current 'test' argument is defined as NULL or "auto" and the user_test argument is NULL then the algorithm automatically selects the best test based on the type of the data. Particularly:
\itemize{
	\item if target is a factor, the multinomial logistic test is used
	\item if target is a ordered factor, the ordered logit regression is used in the logistic test
	\item if target is a numerical vector, the fisher conditional independence test is used
	\item if target is a Surv object, the Survival conditional independence test is used
}

Conditional independence test functions to be pass through the user_test argument should have the same signature of the included test. See "?testIndFisher" for an example.

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".
}
\value{
The output of the algorithm is an object of the class 'SESoutput' including:
\item{selectedVars}{
The selected variables, i.e., the signature of the target variable.
}
\item{selectedVarsOrder}{
The order of the selected variables according to increasing pvalues.
}
\item{queues}{
A list containing a list (queue) of equivalent features for each variable included in selectedVars. An equivalent signature can be built by selecting a single feature from each queue.
}
\item{signatures}{
A matrix reporting all equivalent signatures (one signature for each row).
}
\item{hashObject}{
The hashObject caching the statistic calculted in the current run.
}
\item{pvalues}{
For each feature included in the dataset, this vector reports the strength of its association with the target in the context of all other variables. Particularly, this vector reports the max p-values foudn when the association of each variable with the target is tested against different conditional sets. Lower values indicate higher association.
}
\item{stats}{
The statistics corresponding to "pvalues" (higher values indicates higher association).
}
\item{max_k}{
The max_k option used in the current run.
}
\item{threshold}{
The threshold option used in the current run.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}

Generic Functions implemented for SESoutput Object:
\item{summary(x=SESoutput)}{
Summary view of the SESoutput object.
}
\item{plot(object=SESoutput, mode="all")}{
Plots the generated pvalues (using barplot) of the current SESoutput object in comparison to the threshold.

Argument mode can be either "all" or "partial" for the first 500 pvalues of the object.
}

}
\references{
I. Tsamardinos, V. Lagani and D. Pappas (2012) Discovering multiple, equivalent biomarker signatures. In proceedings of the 7th conference of the Hellenic Society for Computational Biology & Bioinformatics - HSCBB12.
}
\author{
Ioannis Tsamardinos, Vincenzo Lagani

R implementation and documentation: Giorgos Athineou <athineou@ics.forth.gr> Vincenzo Lagani <vlagani@ics.forth.gr>  
}
\note{
The packages required for the SES algorithm operations are: 

\bold{gRbase} : for faster computations in the internal functions

\bold{hash} : for the hash-based implementation

\bold{VGAM} : require(stats) and require(MASS) for the testIndLogistic test

\bold{survival} : for the censIndLR test

\bold{pcalg} : for the gSquare test.
}
\seealso{
\code{\link{CondIndTests}, \link{cv.ses}}
}

\examples{
set.seed(123)
#require(gRbase) #for faster computations in the internal functions
require(hash)

#simulate a dataset with continuous data
dataset <- matrix(nrow = 1000 , ncol = 300)
dataset <- apply(dataset, 1:2, function(i) runif(1, 1, 100))

#define a simulated class variable 
target = 3*dataset[,10] + 2*dataset[,200] + 3*dataset[,20] + runif(1, 0, 1);

#define some simulated equivalences
dataset[,15] = dataset[,10]
dataset[1,10] = dataset[1,10] + 0.2
dataset[,250] = dataset[,200] 
dataset[,230] = dataset[,200] 

require("hash", quietly = TRUE)
{
#run the SES algorithm
sesObject <- SES(target , dataset , max_k=5 , threshold=0.2 , test="testIndFisher", 
hash = TRUE, hashObject=NULL);
#print summary of the SES output
summary(sesObject);
#plot the SES output
plot(sesObject, mode="all");
#get the queues with the equivalences for each selected variable
sesObject@queues
#get the generated signatures
sesObject@signatures;
#get the run time
# > sesObject@runtime;
# user  system elapsed 
# 0.35    0.00    0.35


#re-run the SES algorithm with the same or different configuration 
#under the hash-based implementation of retrieving the statistics
#in the SAME dataset (!important)
hashObj <- sesObject@hashObject;
sesObject2 <- SES(target , dataset , max_k=2 , threshold=0.01 , test="testIndFisher",
hash = TRUE, hashObject=hashObj);
#retrieve the results: summary, plot, sesObject2@...)
summary(sesObject2)
#get the run time
# > sesObject2@runtime;
# user  system elapsed
# 0.01    0.00    0.01
}
}

\keyword{ SES }
\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }
\keyword{ Variable Selection }% __ONLY ONE__ keyword per line
