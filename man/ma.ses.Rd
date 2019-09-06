\name{Constraint based feature selection algorithms for multiple datasets}
\alias{ma.ses}
\alias{ma.mmpc}
\title{
ma.ses: Feature selection algorithm for identifying multiple minimal, statistically-equivalent and equally-predictive feature signatures with multiple datasets

ma.mmpc: Feature selection algorithm for identifying minimal feature subsets with multiple datasets
}

\description{
SES algorithm follows a forward-backward filter approach for feature selection in order to provide minimal, highly-predictive, statistically-equivalent, multiple feature subsets of two or more high dimensional datasets. 
See also Details. MMPC algorithm follows the same approach without generating multiple feature subsets.
}

\usage{
ma.ses(target, dataset, ina, statistic = FALSE, max_k = 3, threshold = 0.05, 
test = NULL, ini = NULL, user_test = NULL, hash = FALSE, hashObject = NULL, 
ncores = 1)
 
ma.mmpc(target, dataset, ina, statistic = FALSE, max_k = 3, threshold = 0.05, 
test = NULL, ini = NULL, user_test = NULL, hash = FALSE, hashObject = NULL, 
ncores = 1, backward = FALSE)
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
\item{dataset}{
The data-set; provide either a data frame or a matrix (columns = variables , rows = samples).
Alternatively, provide an ExpressionSet (in which case rows are samples and columns are features, see bioconductor for details).
}
\item{ina}{
A numerical vector indicating the dataset. The numbers must be 1, 2, 3,...
}
\item{statistic}{ 
A boolean variable indicating whether the test statistics (TRUE) or the p-values should be combined (FALSE). See the details about this. 
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). Integer, default value is 3.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
The conditional independence test to use. Default value is NULL. See also \code{\link{CondIndTests}}.
}
\item{ini}{
This is a supposed to be a list. After running SES or MMPC with some hyper-parameters you might want to run SES again with different hyper-parameters. To avoid calculating the univariate associations (first step of SES and of MPPC) again, you can extract them from the first run of SES and plug them here. This can speed up the second run (and subequent runs of course) by 50\%. See the details and the argument "univ" in the output values. 
}
\item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, the "test" argument is ignored.
}
  \item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to store the statistics calculated during SES execution in a hash-type object. Default value is FALSE. If TRUE a hashObject is produced.
}
  \item{hashObject}{
A List with the hash objects generated in a previous run of SES or MMPC. 
Each time SES runs with "hash=TRUE" it produces a list of hashObjects that can be re-used in order to speed up next runs of SES or MMPC.

Important: the generated hashObjects should be used only when the same dataset is re-analyzed, possibly with different values of max_k and threshold.
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
\item{backward}{
If TRUE, the backward (or symmetry correction) phase will be implemented. This removes any falsely included variables in the parents and children set of the target variable. This is an experimental stage, so do not trust too much or ven better do not use. 
}
}

\details{
This is more at an experimental stage at the present. 

The SES function implements the Statistically Equivalent Signature (SES) algorithm as presented in "Tsamardinos, Lagani and Pappas, HSCBB 2012" 
http://www.mensxmachina.org/publications/discovering-multiple-equivalent-biomarker-signatures/

The MMPC function mplements the MMPC algorithm as presented in "Tsamardinos, Brown and Aliferis. The max-min hill-climbing Bayesian network structure learning algorithm" http://www.dsl-lab.org/supplements/mmhc_paper/paper_online.pdf

The output value "univ" along with the output value "hashObject" can speed up the computations of subesequent runs of SES and MMPC. 
The first run with a specific pair of hyper-parameters (threshold and max_k) the univariate associations tests and the conditional independence tests (test statistic and logarithm of their corresponding p-values) 
are stored and returned. In the next run(s) with different pair(s) of hyper-parameters you can use this information to save time. With a few thousands of variables you will see the difference, 
which can be up to 50\%. For the non robust correlation based tests, the difference may not be significant though, because the unconditional correlation coefficients are calcualted very efficiently. 

The max_k option: the maximum size of the conditioning set to use in the conditioning independence test. Larger values provide more accurate results, at the cost of higher computational times. When the sample size is small (e.g., \eqn{<50} observations) the max_k parameter should be \eqn{ \leq 5}, otherwise the conditional independence test may not be able to provide reliable results.

If the dataset (predictor variables) contains missing (NA) values, they will automatically be replaced by the current variable (column) mean value with an appropriate warning to the user after the execution.

%If the target is a single integer value or a string, it has to corresponds to the column number or to the name of the target feature in the dataset. In any other case the target is a variable that is not contained in the dataset.

%If the current 'test' argument is defined as NULL or "auto" and the user_test argument is NULL then the algorithm automatically selects the best test based on the type of the data. Particularly:
%\itemize{
%	\item if target is a factor, the multinomial or the binary logistic regression is used. If the target has two values only, binary logistic regression will be used.
%	\item if target is a ordered factor, the ordered logit regression is used in the logistic test.
%	\item if target is a numerical vector and the dataset is a matrix or a data.frame with continuous variables, the Fisher conditional independence test is used. If the dataset is a data.frame and there are categorical variables, linear regression is used.
%	\item if target is discrete numerical (counts), the poisson regression conditional independence test is used. If there are only two values, the binary logistic regression is to be used.
%	\item if target is a Surv object, the Survival conditional independence test is used.
%	\item if target is a matrix with at least 2 columns, the multivariate linear regression is used.
%}

Conditional independence test functions to be pass through the user_test argument should have the same signature of the included test. See \code{ \link{testIndFisher} } for an example.

For all the available conditional independence tests that are currently included on the package, please see \code{\link{CondIndTests}}. 

If two or more p-values are below the machine epsilon (.Machine$double.eps which is equal to 2.220446e-16), all of them are set to 0. To make the comparison or the ordering feasible we use the logarithm of the p-value. The max-min heuristic though, requires comparison and an ordering of the p-values. Hence, all conditional independence tests calculate the logarithm of the p-value.

If there are missing values in the dataset (predictor variables) columnwise imputation takes place. The median is used for the continuous variables and the mode for categorical variables. It is a naive and not so clever method. For this reason the user is encouraged to make sure his data contain no missing values. 

If you have percentages, in the (0, 1) interval, they are automatically mapped into \eqn{R} by using the logit transformation. 
% If you set the test to \code{\link{testIndBeta}}, beta regression is used. If you have compositional data, positive multivariate data where each vector sums to 1, with NO zeros, they are also mapped into the Euclidean space using the additive log-ratio (multivariate logit) transformation (Aitchison, 1986).

If you use testIndSpearman (argument "test"), the ranks of the data calculated and those are used in the caclulations. This speeds up the whole procedure.

Currently only the testIndFisher and testIndSpearman tests are supported for use in the algorithm. 

If the argument \bold{statistic} is set to FALSE, the p-values from the hypothesis test of each dataset are combined via Fisher's meta-analytic approach, that is \eqn{ T=-2\sum_{i=1}^k\log{p_i} }
and \eqn{ T~\chi^2_{2k} }. If \bold{statistic} is TRUE, the test statistics are combined as 
\eqn{ T = \frac{{\sum_{i=1}^kt_i/se(t_i)}}{\sum_{i=1}^k 1/se(t_i)}} and \eqn{T~N(0,1)}.
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
A list containing a list (queue) of equivalent features for each variable included in selectedVars. An equivalent signature can be built by selecting a single feature from each queue. Featured only in SES.
}
\item{signatures}{
A matrix reporting all equivalent signatures (one signature for each row). Featured only in SES.
}
\item{hashObject}{
The hashObject caching the statistic calculated in the current run.
}
\item{pvalues}{
For each feature included in the dataset, this vector reports the strength of its association with the target in the context of all other variables. 
Particularly, this vector reports the max p-values (on a logarithmic scale) found when the association of each variable with the target is tested against different conditional sets. Lower values indicate higher association.
}
\item{stats}{
The statistics corresponding to "pvalues" (higher values indicates higher association).
}
\item{univ}{
This is a list with the univariate associations. The test statistics and their corresponding logged p-values, along with their flag (1 if the test was perfromed and 0 otherwise). This list is very important for subsequent runs of SES with different hyper-parameters. After running SES with some hyper-parameters you might want to run SES again with different hyper-parameters. To avoid calculating the univariate associations (first step of SES or MMPC) again, you can take this list from the first run of SES and plug it in the argument "ini" in the next run(s) of SES or MMPC. This can speed up the second run (and subequent runs of course) by 50\%. See the argument "univ" in the output values. 
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
\item{test}{
The character name of the statistic test used.
}

Generic Function implemented for SESoutput Object:
\item{plot(object=SESoutput, mode="all")}{
Plots the generated pvalues (using barplot) of the current SESoutput object in comparison to the threshold.

Argument mode can be either "all" or "partial" for the first 500 pvalues of the object.
}
}

\references{
V. Lagani, A.D. Karozou, D. Gomez-Cabrero, G. Silberberg and I. Tsamardinos (2016). A comparative evaluation of data-merging and meta-analysis methods for reconstructing gene-gene interactions, BMC Bioinformatics 17(Supplementary 5): 287-305. 

I. Tsamardinos, V. Lagani and D. Pappas (2012). Discovering multiple, equivalent biomarker signatures. In proceedings of the 7th conference of the Hellenic Society for Computational Biology & Bioinformatics - HSCBB12.

Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.
}

\author{
Ioannis Tsamardinos, Vincenzo Lagani

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr>, Vincenzo Lagani <vlagani@csd.uoc.gr>  and Michail Tsagris <mtsagris@uoc.gr>
}


\seealso{
\code{\link{SES}, \link{CondIndTests}, \link{cv.ses}}
}

\examples{
set.seed(123)

#simulate a dataset with continuous data
dataset <- matrix(runif(400 * 100, 1, 100), ncol = 100)

#define a simulated class variable 
target <- 3 * dataset[, 10] + 2 * dataset[, 20] + 3 * dataset[, 30] + rnorm(400, 0, 5)

#define some simulated equivalences
dataset[, 15] <- dataset[, 10] + rnorm(400, 0, 2)
dataset[, 10] <- dataset[ , 10] + rnorm(400, 0, 2)
dataset[, 25] <- dataset[, 20] + rnorm(400, 0, 2) 
dataset[, 23] <- dataset[, 20] + rnorm(400, 0, 2)

#run the SES algorithm
a1 <- SES(target , dataset, max_k = 5, threshold = 0.05, test = "testIndFisher", 
hash = TRUE, hashObject = NULL)

ina <- rbinom(400, 2, 0.5) + 1
a2 <- ma.ses(target , dataset, ina = ina, max_k = 5, threshold = 0.05, test = "testIndFisher", 
hash = TRUE, hashObject = NULL)
a3 <- ma.mmpc(target , dataset, ina = ina, max_k = 5, threshold = 0.05, test = "testIndFisher", 
hash = TRUE, hashObject = NULL)

#get the generated signatures
a1@signatures
a2@signatures
a3@selectedVars

}

\keyword{ SES }
\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }
\keyword{ Variable Selection }
