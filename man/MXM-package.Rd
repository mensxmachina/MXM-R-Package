\name{MXM-package}
\alias{MXM-package}
\docType{package}
\title{
This is an R package that currently implements feature selection methods for identifying minimal, 
statistically-equivalent and equally-predictive feature subsets. In addition, two algorithms for 
constructing the skeleton of a Bayesian network are included.
}

\description{
'MXM' stands for Mens eX Machina, meaning 'Mind from the Machine' in Latin. The package provides source code for the SES 
algorithm and for some appropriate statistical conditional independence tests. (Fisher and Spearman correlation, 
G-square test are some examples. Currently the response variable can be univariate or multivariate Euclidean, 
proportions within 0 and 1, compositional data without zeros and  ones, binary, nominal or ordinal multinomial, 
count data (handling also overdispersed and with more zeros than expected), longitudinal, clustered data, survival 
and case-control. Robust versions are also available in some cases and a K-fold cross validation is offered. 
Bayesian network related algorithms and ridge reression are also included. Read the package's help pages for more details.

MMPC and SES can handle even thousands of variables and for some tests, even many sample sizes of tens of thousands. 
The user is best advised to check his variables in the beginning. For some regressions, logistic and Poisson for example, 
we have used C++ codes for speed reasons. 

For more information the reader is addressed to Lagani V., Athineou G., Farcomeni A., Tsagris M. and Tsamardinos I. (2017). 
Feature Selection with the R Package MXM: Discovering Statistically Equivalent Feature Subsets. Journal of Statistical Software, 
80(7), doi:10.18637/jss.v080.i07.
}

\details{
\tabular{ll}{
Package: \tab MXM\cr
Type: \tab Package\cr
Version: \tab 1.4.4\cr
Date: \tab 2019-06-19\cr
License: \tab GPL-2\cr
}
}

\author{
Giorgos Athineou <athineou@csd.uoc.gr>, Vincenzo Lagani <vlagani@csd.uoc.gr>, Michail Tsagris <mtsagris@uoc.gr>, Giorgos Borboudakis <borbudak@csd.uoc.gr>, Ioannis Tsamardinos <tsamard@csd.uoc.gr>, Anna Roumpelaki <anna.roumpelaki@gmail.com>
}

\section{Maintainer}{
Michail Tsagris <mtsagris@uoc.gr>
}

\references{
Borboudakis G. and Tsamardinos I. (2019). Forward-backward selection with early dropping. Journal of Machine Learning Research, 20(8): 1-39.

Chen S., Billings S. A., and Luo W. (1989). Orthogonal least squares methods and their application to non-linear system identification. 
International Journal of control, 50(5), 1873-1896.
http://eprints.whiterose.ac.uk/78100/1/acse%20report%20343...pdf

Davis G. (1994). Adaptive Nonlinear Approximations. PhD thesis. 
http://www.geoffdavis.net/papers/dissertation.pdf

Demidenko E. (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.

Gharavi-Alkhansari M., anz Huang T. S. (1998, May). A fast orthogonal matching pursuit algorithm. In Acoustics, Speech and Signal Processing, 1998. 
Proceedings of the 1998 IEEE International Conference on (Vol. 3, pp. 1389-1392). 

Lagani V., Athineou G., Farcomeni A., Tsagris M. and Tsamardinos I. (2017). Feature Selection with the R Package MXM: Discovering Statistically 
Equivalent Feature Subsets. Journal of Statistical Software, 80(7), doi:10.18637/jss.v080.i07.

Lagani V., Kortas G. and Tsamardinos I. (2013), Biomarker signature identification in "omics" with multiclass outcome. 
Computational and Structural Biotechnology Journal, 6(7):1-7.

Liang  K.Y.  and  Zeger  S.L. (1986). Longitudinal data analysis using generalized linear models. Biometrika, 73(1): 13-22.

Mallat S. G. & Zhang Z. (1993). Matching pursuits with time-frequency dictionaries. IEEE Transactions on signal processing, 41(12), 3397-3415.
https://www.di.ens.fr/~mallat/papiers/MallatPursuit93.pdf

Paik M.C. (1988). Repeated measurement analysis for nonnormal data in small samples. Communications in
Statistics-Simulation and Computation, 17(4): 1155-1171.

Pati Y. C., Rezaiifar R. and Krishnaprasad P. S. (1993). Orthogonal matching pursuit: Recursive function approximation with applications to wavelet 
decomposition. In Signals, Systems and Computers. 1993 Conference Record of The Twenty-Seventh Asilomar Conference on. IEEE.

Prentice R.L. and Zhao L.P. (1991). Estimating equations for parameters in means and covariances
of multivariate discrete and continuous responses. Biometrics, 47(3): 825-839.

Spirtes P.,  Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, 3nd edition.

Tsagris M. (2019). Bayesian network learning with the PC algorithm: an improved and correct variation. 
Applied Artificial Intelligence 33(2): 101-123.

Tsagris, M., Papadovasilakis, Z., Lakiotaki, K. and Tsamardinos, I. (2018). Efficient feature selection on gene expression data: 
Which algorithm to use? BioRxiv preprint. 

Tsagris M. (2018). Guide on performing feature selection with the R package MXM.
http://mensxmachina.org/site/wp-content/uploads/2018/04/Guide-on-performing-feature-selection-with-the-R-package-MXM.pdf

Tsagris M., Lagani V., and Tsamardinos I. (2018). Feature selection for high-dimensional temporal data. BMC bioinformatics, 19:17.

Tsagris M., Borboudakis G., Lagani V. and Tsamardinos I. (2018). Constraint-based Causal Discovery with Mixed Data. 
International Journal of Data Science and Analytics. 

Tsamardinos I., Greasidou E. and Borboudakis G. (2018). Bootstrapping the out-of-sample predictions for efficient and accurate cross-validation. 
Machine Learning 107(12): 1895-1922.  

Tsamardinos I., Lagani V. and Pappas D. (2012) Discovering multiple, equivalent biomarker signatures. In proceedings of the 7th conference of the 
Hellenic Society for Computational Biology \& Bioinformatics - HSCBB12.

Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.

Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal relations. 
In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining p. 673-678. 

Yan J. and Fine J. (2004). Estimating equations for association structures. Statistics in medicine, 23(6): 859-874

Zhang, Jiji. (2008). On the completeness of orientation rules for causal discovery in the presence of latent confounders and selection bias. 
Artificial Intelligence 172(16): 1873--1896.

Ziegler A., Kastner C., Brunner D. and Blettner M. (2000). Familial associations of lipid profiles: A generalised
estimating equations approach. Statistics in medicine, 19(24): 3345-3357

}

\note{
Acknowledgments:
The research leading to these results has received funding from the European Research Council under the 
European Union's Seventh Framework Programme (FP/2007-2013) / ERC Grant Agreement n. 617393.

Michail Tsagris would like to express his acknowledgments to Marios Dimitriadis and Manos Papadakis, 
undergraduate students in the department of computer science, university of Crete, for their programming tips and advice. 
Their help has been very valuable. Dr David Gomez Cabrero Lopez is also acknowledged. Dr Uwe Ligges and Prof Kurt Hornik from the CRAN 
team are greatly acknowledged for their assistance. Prof Achim Zeileis is grealty acknowledged for this help with the quasi Poisson and 
quasi binomial regression models. Christina Chatzipantsiou and Kleio Maria Verrou are acknowledged for their suggestions. Margarita Rebolledo 
is acknowledged for spotting a bug. Nikolaos Pandis from the University of Bern is acknowledged for his suggestion of the AFT (regression) 
models and for his suggestions. Michail is grateful to James Dalgleish from Columbia University who suggested that we mention, in various places, 
that most algorithms return the logarithm of the p-values and not the p-values.

\bold{Disclaimer:} Professor Tsamardinos is the creator of this package. Mr Athineou built it, supervised by Dr Lagani. 
Dr Tsagris is the current maintainer.
}

\keyword{ SES }
\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }

\seealso{
\code{\link{SES}, \link{MMPC}, \link{censIndCR},\link{testIndFisher}, \link{testIndLogistic}, \link{gSquare}, \link{testIndRQ}}
}

