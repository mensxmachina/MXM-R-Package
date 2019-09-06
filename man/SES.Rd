\name{Constraint based feature selection algorithms}
\alias{SES}
\alias{MMPC}
\alias{wald.ses}
\alias{wald.mmpc}
\alias{perm.ses}
\alias{perm.mmpc}
\title{
SES: Feature selection algorithm for identifying multiple minimal, statistically-equivalent and equally-predictive feature signatures

MMPC: Feature selection algorithm for identifying minimal feature subsets
}

\description{
SES algorithm follows a forward-backward filter approach for feature selection in order to provide minimal, highly-predictive, statistically-equivalent, multiple feature subsets of a high dimensional dataset. See also Details. MMPC algorithm follows the same approach without generating multiple feature subsets.
}

\usage{
SES(target, dataset, max_k = 3, threshold = 0.05, test = NULL, ini = NULL, 
wei = NULL, user_test = NULL, hash = FALSE, hashObject = NULL,
ncores = 1, backward = FALSE)
 
MMPC(target, dataset, max_k = 3, threshold = 0.05, test = NULL, ini = NULL, 
wei = NULL, user_test = NULL, hash = FALSE, hashObject = NULL, 
ncores = 1, backward = FALSE)

wald.ses(target, dataset, max_k = 3, threshold = 0.05, test = NULL, ini = NULL, 
wei = NULL, user_test = NULL, hash = FALSE, hashObject = NULL, 
ncores = 1, backward = FALSE)

wald.mmpc(target, dataset, max_k = 3, threshold = 0.05, test = NULL, ini = NULL, 
wei = NULL, user_test = NULL, hash = FALSE, hashObject = NULL,  
ncores = 1, backward = FALSE)

perm.ses(target, dataset, max_k = 3, threshold = 0.05, test = NULL, ini = NULL, 
wei = NULL, user_test = NULL, hash=FALSE, hashObject = NULL, 
R = 999, ncores = 1, backward = FALSE)

perm.mmpc(target, dataset, max_k = 3, threshold = 0.05, test = NULL, ini = NULL, 
wei = NULL, user_test = NULL, hash=FALSE, hashObject = NULL,
R = 999, ncores = 1, backward = FALSE)
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
\item{dataset}{
The data-set; provide either a data frame or a matrix (columns = variables, rows = samples).
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
This is a supposed to be a list. After running SES or MMPC with some hyper-parameters you might want to run SES again with different hyper-parameters. To avoid calculating the univariate associations (first step of SES and of MMPC) again, you can extract them from the first run of SES and plug them here. This can speed up the second run (and subequent runs of course) by 50\%. See the details and the argument "univ" in the output values. 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
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
\item{R}{
The number of permutations, set to 999 by default. There is a trick to avoind doing all permutations. As soon as the number of times the permuted test statistic is more than the observed test statistic is more than 50 (if threshold = 0.05 and R = 999), the p-value has exceeded the signifiance level (threshold value) and hence the predictor variable is not significant. There is no need to continue do the extra permutations, as a decision has already been made. 
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
\item{backward}{
If TRUE, the backward (or symmetry correction) phase will be implemented. This removes any falsely included variables in the parents and children set of the target variable. It call the \code{link{mmpcbackphase}} for this purpose. 
}
}

\details{
The SES function implements the Statistically Equivalent Signature (SES) algorithm as presented in "Tsamardinos, Lagani and Pappas, HSCBB 2012" 
http://www.mensxmachina.org/publications/discovering-multiple-equivalent-biomarker-signatures/

The MMPC function implements the MMPC algorithm as presented in "Tsamardinos, Brown and Aliferis. The max-min hill-climbing Bayesian network structure learning algorithm" 
http://www.dsl-lab.org/supplements/mmhc_paper/paper_online.pdf

he output value "univ" along with the output value "hashObject" 
can speed up the computations of subsequent runs of SES and MMPC. The first run with a specific pair of hyper-parameters (threshold and max_k) the univariate associations tests 
and the conditional independence tests (test statistic and natural logarithm of their corresponding p-values) are stored and returned. In the next run(s) with different pair(s) of 
hyper-parameters you can use this information to save time. With a few thousands of variables you will see the difference, which can be up to 50\%. For the non robust correlation based tests, 
the difference may not be significant though, because the unconditional correlation coefficients are calculated very efficiently.s. 

The max_k option: the maximum size of the conditioning set to use in the conditioning independence test. Larger values provide more accurate results, at the cost of higher 
computational times. When the sample size is small (e.g., \eqn{<50} observations) the max_k parameter should be \eqn{ \leq 5}, otherwise the conditional independence test may 
not be able to provide reliable results.

If the dataset (predictor variables) contains missing (NA) values, they will automatically be replaced by the current variable (column) mean value with an appropriate warning 
to the user after the execution.

If the target is a single integer value or a string, it has to corresponds to the column number or to the name of the target feature in the datase. In any other case the target is a variable that is not contained in the dataset.

If the current 'test' argument is defined as NULL or "auto" and the user_test argument is NULL then the algorithm automatically selects the best test based on the type of the data. 
Particularly:
\itemize{
	\item if the target is a factor, the multinomial or the binary logistic regression is used. If the target has two values only, binary logistic regression will be used. 
	\item if target is a ordered factor, ordinal regression is used in the logistic test. Hence, if you want to use multinomial or ordinal logistic regression, make sure your 
	target is factor. 
	\item if target is a numerical vector and the dataset is a matrix or a data.frame with continuous variables, the Fisher conditional independence test is used. If the dataset is a 
	data.frame and there are categorical variables, linear regression is used.
	\item if target is discrete numerical (counts), the Poisson regression conditional independence test is used. If there are only two values, the binary logistic regression is to 
	be used.
	\item if target is a Surv object, a Survival conditional independence test is used.
	\item if target is a matrix with at least 2 columns, the multivariate linear regression is used. 
	\item if target is a 2 column matrix whose columns are the number of successes and the number of trials (first and second column respectively) the testIndBinom should be used.
}

Conditional independence test functions to be pass through the user_test argument should have the same signature of the included test. See \code{\link{testIndFisher}} for an example.

For all the available conditional independence tests that are currently included on the package, please see \code{\link{CondIndTests}}. 
If two or more p-values are below the machine epsilon (.Machine$double.eps which is equal to 2.220446e-16), all of them are set to 0. To make the comparison or the ordering feasible 
we use the logarithm of the p-value. The max-min heuristic though, requires comparison and an ordering of the p-values. Hence, all conditional independence tests calculate the logarithm 
of the p-value.

If there are missing values in the dataset (predictor variables) columnwise imputation takes place. The median is used for the continuous variables and the mode for categorical variables. It is a naive and not so clever method. For this reason the user is encouraged to make sure his data contain no missing values. 

If you have percentages, in the (0, 1) interval, they are automatically mapped into \eqn{R} by using the logit transformation. If you set the test to \code{\link{testIndBeta}}, 
beta regression is used. If you have compositional data, positive multivariate data where each vector sums to 1, with NO zeros, they are also mapped into the Euclidean space using the additive log-ratio (multivariate logit) transformation (Aitchison, 1986).

If you use testIndSpearman (argument "test"), the ranks of the data calculated and those are used in the caclulations. This speeds up the whole procedure.

As a rule of thumb you can try this. If for example you have counts and want to see which model fits best, there are two ways. Calculate the mean and the variance. If they are similar, 
use the Poisson instead of the negative binomial as it is much faster. If you are not convinced, you can either use the negative binomial or do the following simulation study.

# x <- matrix(rnorm(n * 1000), ncol = 1000)
# a <- Rfast::univglms(y, x)
# hist(a[, 2])  ## histogram of the p-values

If the histogram shows a uniform distribution, use the Poisson regression. If the histogram is not uniform, then repeat the simluation but with a negative binomial distribution. 
If the histogram is again not flat, then another model is necessary. If the data come from a Poisson or negative binomial, the histogram with a negative binomial regressino will be flat.
If the data come a from a negative binomial, the histogram with a Poisson will not be uniform. 

On the same page, if you have many zeros, try Rfast::zip.mle and see whether there are grounds to to facilitate the use of a zero inflated Poisson model. Otherwise, do a simulation 
study like before.  
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
For each feature included in the dataset, this vector reports the strength of its association with the target in the context of all other variable. 
Particularly, this vector reports the max p-values found when the association of each variable with the target is tested against different conditional sets. 
Lower values indicate higher association. Note that these are the logged p-values, natural logarithm of the pvalues, and not the p-values.
}
\item{stats}{
The statistics corresponding to "pvalues" (higher values indicates higher association).
}
\item{univ}{
This is a list with the univariate associations; the test statistics and their corresponding logged p-values. This list is very important for subsequent runs of SES with different hyper-parameters. After running SES with some hyper-parameters you might want to run SES again with different hyper-parameters. To avoid calculating the univariate associations (first step of SES or MMPC) again, you can take this list from the first run of SES and plug it in the argument "ini" in the next run(s) of SES or MMPC. This can speed up the second run (and subequent runs of course) by 50\%. See the argument "univ" in the output values. 
}
\item{max_k}{
The max_k option used in the current run.
}
\item{threshold}{
The threshold option used in the current run.
}
\item{n.tests}{
If you have set hash = TRUE, then the number of tests performed by SES or MMPC will be returned. If you have not set this to TRUE, the number of univariate associations will be returned. So be careful with this number.
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
Feature Selection with the R Package MXM: Discovering Statistically Equivalent Feature Subsets, Lagani, V. and Athineou, 
G. and Farcomeni, A. and Tsagris, M. and Tsamardinos, I. (2017). Journal of Statistical Software, 80(7).

I. Tsamardinos, V. Lagani and D. Pappas (2012). Discovering multiple, equivalent biomarker signatures. In proceedings of the 
7th conference of the Hellenic Society for Computational Biology & Bioinformatics - HSCBB12.

Tsamardinos, I., Aliferis, C. F., & Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal 
relations. In Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). ACM.

Brown, L. E., Tsamardinos, I., & Aliferis, C. F. (2004). A novel algorithm for scalable and accurate Bayesian network learning. 
Medinfo, 711-715.

Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 
65(1), 31-78.

Tsagis M. (2018). Guide on performing feature selection with the R package MXM.
http://mensxmachina.org/site/wp-content/uploads/2018/04/Guide-on-performing-feature-selection-with-the-R-package-MXM.pdf
}

\author{
Ioannis Tsamardinos, Vincenzo Lagani

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Vincenzo Lagani <vlagani@csd.uoc.gr>  
}

\note{
The packages required by the SES and MMPC algorithm operations are: 

\bold{quantreg}: for the quantile (median) regression

\bold{MASS}: for negative binomial regression and simple ordinal regression

\bold{nnet} : also require(stats) and require(MASS) for the testIndLogistic test

\bold{survival} : for the censIndCR, censIndWR and the censIndER tests

\bold{doParallel}: for parallel computations

\bold{lme4}: for (generalised) linear mixed models

\bold{Rfast}: for many fast functions.
}

\seealso{
\code{\link{CondIndTests}, \link{cv.ses}}
}

\examples{
set.seed(123)

#simulate a dataset with continuous data
dataset <- matrix(runif(1000 * 200, 1, 100), ncol = 200)

#define a simulated class variable 
target <- 3 * dataset[, 10] + 2 * dataset[, 100] + 3 * dataset[, 20] + rnorm(1000, 0, 5)

# define some simulated equivalences
dataset[, 15] <- dataset[, 10] + rnorm(1000, 0, 2)
dataset[, 150] <- dataset[, 100] + rnorm(1000, 0, 2) 
dataset[, 130] <- dataset[, 100] + rnorm(1000, 0, 2)

# run the SES algorithm
sesObject <- SES(target , dataset, max_k = 5, threshold = 0.05, test = "testIndFisher", 
hash = TRUE, hashObject = NULL);

# get the queues with the equivalences for each selected variable
sesObject@queues
#get the generated signatures
sesObject@signatures;

# re-run the SES algorithm with the same or different configuration 
# under the hash-based implementation of retrieving the statistics
# in the SAME dataset (!important)
hashObj <- sesObject@hashObject;
sesObject2 <- SES(target, dataset, max_k = 2, threshold = 0.01, test = "testIndFisher", 
hash = TRUE, hashObject = hashObj);

# get the run time
sesObject@runtime;
sesObject2@runtime;

# MMPC algorithm 
mmpcObject <- MMPC(target, dataset, max_k = 3, threshold = 0.05, test="testIndFisher");
mmpcObject@selectedVars
mmpcObject@runtime
}

\keyword{ SES }
\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }
\keyword{ Variable Selection }
