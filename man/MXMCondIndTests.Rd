\name{CondInditional independence tests}
\alias{CondIndTests}

\title{ MXM Conditional independence tests }
\description{

Currently the \pkg{MXM} package supports numerous tests for different types of target (dependent) and predictor (independent) variables. The target variable can be of continuous, discrete, categorical and of survival type. As for the predictor variables, they can be continuous, categorical or mixed. 

The \bold{testIndFisher} and the \bold{gSquare} tests have two things in common. They do not use a model implicitly (i.e. estimate some beta coefficients), even though there is an underlying assumed one. Secondly they are pure tests of independence (again, with assumptions required).

As for the other tests, they share one thing in common. For all of them, two parametric models must be fit. The null model containing the conditioning set of variables alone and the alternative model containing the conditioning set and the candidate variable. The significance of the new variable is assessed via a log-likelihood ratio test with the appropriate degrees of freedom. All of these tests which are available for SES and MMPC are summarized in the below table.

\tabular{llll}{ 
  \bold{Target variable} \tab	\bold{Predictor variables} \tab \bold{Available tests} \tab	\bold{Short explanation} \cr 
  
  Continuous \tab Continuous \tab testIndFisher \tab Partial correlation \cr 
  Continuous \tab Continuous \tab testIndMMFisher \tab Robust partial correlation \cr 
  Continuous \tab Continuous \tab testIndSpearman \tab Partial correlation \cr 
  Continuous \tab Mixed \tab testIndMMReg \tab MM regression \cr  
  Continuous \tab Mixed \tab testIndRQ \tab Median regression \cr  
  Proportions \tab Continuous \tab testIndFisher \tab Partial correlation \cr 
  Proportions \tab Continuous \tab testIndMMFisher \tab Robust partial correlation \cr 
  Proportions \tab Continuous \tab testIndSpearman \tab Partial correlation \cr 
  Proportions \tab Mixed \tab testIndReg \tab Linear regression \cr 
  Proportions \tab Mixed \tab testIndMMReg \tab MM regression \cr 
  Proportions \tab Mixed \tab testIndRQ \tab Median regression \cr 
  Proportions \tab Mixed \tab testIndBeta \tab Beta regression \cr 
  Proportions \tab Mixed \tab testIndQbinom \tab Quasi binomial regression \cr 
  Strictly positive \tab Mixed \tab testIndIGreg \tab Inverse Gaussian regression \cr 
  Strictly positive \tab Mixed \tab testIndGamma \tab Gamma regression \cr 
  Non negative \tab Mixed \tab testIndNormLog \tab Gaussian regression with log link \cr 
  Strictly Positive \tab Mixed \tab censIndWR \tab Weibull regression \cr
  Strictly Positive \tab Mixed \tab censIndER \tab Exponential regression \cr
  Strictly Positive \tab Mixed \tab censIndLLR \tab Log-logistic regression \cr
  Successes \& totals \tab Mixed \tab testIndBinom \tab Binomial regression \cr
  Discrete \tab	Mixed \tab testIndPois \tab Poisson regression \cr 
  Discrete \tab	Mixed \tab testIndZIP \tab Zero Inflated \cr 
  \tab	\tab  \tab Poisson regression \cr 
  Discrete \tab Mixed \tab testIndNB \tab Negative binomial regression \cr 
  Discrete \tab Mixed \tab testIndQPois \tab Quasi Poisson regression \cr 
  Factor with two \tab Mixed \tab testIndLogistic \tab Binary logistic regression \cr
  levels or binary \tab \tab \tab \cr 
  Factor with two \tab Mixed \tab testIndQBinom \tab Quasi binomial regression \cr
  levels or binary \tab \tab \tab \cr 
  Factor with more \tab Mixed \tab  testIndMultinom  \tab Multinomial logistic regression  \cr		 
  than two levels \tab  \tab  \tab  \cr		 
  (unordered) \tab \tab \tab  \cr		   
  Factor with more than \tab Mixed \tab testIndOrdinal \tab Ordinal logistic regression \cr		 
  two levels (ordered) \tab \tab \tab \cr  
  Categorical \tab Categorical \tab gSquare \tab G-squared test of independence \cr  
  Categorical \tab Categorical \tab testIndMultinom \tab Multinomial logistic regression \cr  
  Categorical \tab Categorical \tab testIndOrdinal \tab Ordinal logistic regression \cr  
  Survival	\tab Mixed \tab censIndCR \tab Cox regression \cr
  Survival	\tab Mixed \tab censIndWR \tab Weibull regression \cr
  Survival	\tab Mixed \tab censIndER \tab Exponential regression \cr
    Survival	\tab Mixed \tab censIndLLR \tab Log-logistic regression \cr
  Left censored	\tab Mixed \tab testIndTobit \tab Tobit regression \cr
  Case-control  \tab Mixed \tab testIndClogit \tab Conditional logistic regression \cr
  Multivariate continuous \tab Mixed \tab testIndMVreg  \tab Multivariate linear regression \cr 
  Compositional data \tab Mixed \tab testIndMVreg \tab Multivariate linear regression \cr 
  (no zeros) \tab  \tab after multivariate \tab  \cr 
  \tab \tab logit transformation \tab \cr
  Longitudinal/clustered	\tab Continuous \tab testIndGLMMReg \tab Linear mixed models \cr
  Clustered	\tab Continuous \tab testIndLMM   \tab Fast linear mixed models \cr
  Binary longitudinal	\tab Continuous \tab testIndGLMMLogistic \tab Logistic mixed regression \cr
  and clustered	\tab  \tab  \tab   \cr
  Count longitudinal	\tab Continuous \tab testIndGLMMPois \tab Poisson mixed regression \cr
  and clustered	\tab  \tab  \tab  \cr
  Positive longitudinal	\tab Continuous \tab testIndGLMMNormLog \tab GLMM with Gaussian regression \cr
  and clustered	\tab  \tab  \tab and log link  \cr
  Non negative longitudinal	\tab Continuous \tab testIndGLMMGamma \tab GLMM with Gamma regression \cr
  and clustered	\tab  \tab  \tab and log link  \cr
  Longitudinal/clustered	\tab Continuous \tab testIndGEEReg \tab GEE with Gaussian regression \cr
  Binary longitudinal	\tab Continuous \tab testIndGEELogistic \tab GEE with logistic regression \cr
  and clustered	\tab  \tab  \tab   \cr
  Count longitudinal	\tab Continuous \tab testIndGEEPois \tab GEE with Poisson regression \cr
  and clustered	\tab  \tab  \tab  \cr
  Positive longitudinal	\tab Continuous \tab testIndGEENormLog \tab GEE with Gaussian regression \cr
  and clustered	\tab  \tab  \tab and log link  \cr
  Non negative longitudinal	\tab Continuous \tab testIndGEEGamma \tab GEE with Gamma regression \cr
  and clustered	\tab  \tab  \tab and log link  \cr
  Clustered survival	\tab Contiunous \tab testIndGLMMCR \tab Mixed effects Cox regression \cr
  Circular \tab Continuous \tab testIndSPML \tab Circular-linear regression \cr 

}
}

\section{Log-likelihood ratio tests}{
\enumerate{
\item \bold{testIndFisher}. This is a standard test of independence when both the target and the set of predictor variables are continuous (continuous-continuous). When the joint multivariate normality of all the variables is assumed, we know that if a correlation is zero this means that the two variables are independent. Moving in this spirit, when the partial correlation between the target variable and the new predictor variable conditioning on a set of (predictor) variables is zero, then we have evidence to say they are independent as well. An easy way to calculate the partial correlation between the target and a predictor variable conditioning on some other variables is to regress the both the target and the new variable on the conditioning set. The correlation coefficient of the residuals produced by the two regressions equals the partial correlation coefficient. If the robust option is selected, the two aforementioned regression models are fitted using M estimators (Marona et al., 2006). If the target variable consists of proportions or percentages (within the (0, 1) interval), the logit transformation is applied beforehand.
    
\item \bold{testIndSpearman}. This is a non-parametric alternative to \bold{testIndFisher} test. It is a bit slower than its competitor, yet very fast and suggested when normality assumption breaks down or outliers are present. In fact, within SES, what happens is that the ranks of the target and of the dataset (predictor variables) are computed and the \bold{testIndSpearman} is aplied. This is faster than applying Fisher with M estimators as described above. If the target variable consists of proportions or percentages (within the (0, 1) interval), the logit transformation is applied beforehand.
    
\item \bold{testIndReg}. In the case of target-predictors being continuous-mixed or continuous-categorical, the suggested test is via the standard linear regression. In this case, two linear regression models are fitted. One with the conditioning set only and one with the conditioning set plus the new variable. The significance of the new variable is assessed via the F test, which calculates the residual sum of squares of the two models. The reason for the F test is because the new variable may be categorical and in this case the t test cannot be used. It makes sense to say, that this test can be used instead of the \bold{testIndFisher}, but it will be slower. If the robust option is selected, the two models are fitted using M estimators (Marona et al. 2006). If the target variable consists of proportions or percentages (within the (0, 1) interval), the logit transformation is applied beforehand.
    
\item \bold{testIndRQ}. An alternative to \bold{testIndReg} for the case of continuous-mixed (or continuous-continuous) variables is the \bold{testIndRQ}. Instead of fitting two linear regression models, which model the expected value, one can choose to model the median of the distribution (Koenker, 2005). The significance of the new variable is assessed via a rank based test calibrated with an F distribution (Gutenbrunner et al., 1993). The reason for this is that we performed simulation studies and saw that this type of test attains the type I error in contrast to the log-likelihood ratio test. The benefit of this regression is that it is robust, in contrast to the classical linear regression. If the target variable consists of proportions or percentages (within the (0, 1) interval), the logit transformation is applied beforehand.

\item \bold{testIndBeta}. When the target is proportion (or percentage, i.e., between 0 and 1, not inclusive) the user can fit a regression model assuming a beta distribution. The predictor variables can be either continuous, categorical or mixed. The procedure is the same as in the \bold{testIndReg} case. 
    
\item \bold{Alternatives to testIndBeta}. Instead of \bold{testIndBeta} the user has the option to choose all the previous to that mentioned tests by transforming the target variable with the logit transformation. In this way, the support of the target becomes the whole of R^d and then depending on the type of the predictors and whether a robust approach is required or not, there is a variety of alternative to beta regression tests. 

\item \bold{testIndIGreg}. When you have non negative data, i.e. the target variable takes positive values (including 0), a suggested regression is based on the the inverse gaussian distribution. The link function is not the inverse of the square root as expected, but the logarithm. This is to ensure that the fitted values will be always be non negative. The predictor variables can be either continuous, categorical or mixed. The significance between the two models is assessed via the log-likelihood ratio test. Alternatively, the user can use the Weibull regression (\bold{censIndWR}), gamma regression (\bold{testIndGamma}) or Gaussian regression with log link (\bold{testIndNormLog}). 

\item \bold{testIndGamma}. This is an alternative to \bold{testIndIGreg}.

\item \bold{testIndNormLog}. This is a second alternative to \bold{testIndIGreg}.

\item \bold{testIndPois}. When the target is discrete, and in specific count data, the default test is via the Poisson regression. The predictor variables can be either continuous, categorical or mixed. The procedure is the same as in all the previously regression model based tests, i.e. the log-likelihood ratio test is used to assess the conditional independence of the variable of interest. 
    
\item \bold{testIndNB}. As an alternative to the Poisson regression, we have included the Negative binomial regression to capture cases of overdispersion. The predictor variables can be either continuous, categorical or mixed. 

\item \bold{testIndQPois}. This is a better alternative for discrete target, better than the \bold{testIndPois} and than the \bold{testIndNB}, because it can capture both cases of overdispersion and undersidpesrion.

\item \bold{testIndZIP}. When the number of zeros is more than expected under a Poisson model, the zero inflated poisson regression is to be employed. The predictor variables can be either continuous, categorical or mixed. 
    
\item \bold{testIndLogistic}. When the target is categorical with only two outcomes, success or failure for example, then a binary logistic regression is to be used. Whether regression or classification is the task of interest, this method is applicable. The advantage of this over a linear or quadratic discriminant analysis is that it allows for categorical predictor variables as well and for mixed types of predictors.

\item \bold{testIndQBinom}. This is an alternative to either the \bold{testIndLogistic} or especially the \bold{testIndBeta}.
    
\item \bold{testIndMultinom}. If the target has more than two outcomes, but it is of nominal type, there is no ordering of the outcomes, multinomial logistic regression will be employed. Again, this regression is suitable for classification purposes as well and it to allows for categorical predictor variables.
    
\item \bold{testIndOrdinal}. This is a special case of multinomial regression, in which case the outcomes have an ordering, such as \bold{not satisfied}, \bold{neutral}, \bold{satisfied}. The appropriate method is ordinal logistic regression. 
    
\item \bold{testIndBinom}. When the target variable is a matrix of two columns, where the first one is the number of successes and the second one is the number of trials, binomial regression is to be used. 

\item \bold{gSquare}. If all variables, both the target and predictors are categorical the default test is the G-square test of independence. It is similar to the chi-squared test of independence, but instead of using the chi-squared metric between the observed and estimated frequencies in contingency tables, the Kullback-Leibler divergence of the observed from the estimated frequencies is used. The asymptotic distribution of the test statistic is a chi-squared distribution on some appropriate degrees of freedom. The target variable can be either ordered or unordered with two or more outcomes. 
    
\item \bold{Alternatives to gSquare}. An alternative to the \bold{gSquare} test is the \bold{testIndLogistic}. Depending on the nature of the target, binary, un-ordered multinomial or ordered multinomial the appropriate regression model is fitted. 
    
\item \bold{censIndCR}. For the case of time-to-event data, a Cox regression model is employed. The predictor variables can be either continuous, categorical or mixed. Again, the log-likelihood ratio test is used to assess the significance of the new variable. 
    
\item \bold{censIndWR}. A second model for the case of time-to-event data, a Weibull regression model is employed. The predictor variables can be either continuous, categorical or mixed. Again, the log-likelihood ratio test is used to assess the significance of the new variable. Unlike the semi-parametric Cox model, the Weibull model is fully parametric.

\item \bold{censIndER}. A third model for the case of time-to-event data, an exponential regression model is employed. The predictor variables can be either continuous, categorical or mixed. Again, the log-likelihood ratio test is used to assess the significance of the new variable. This is a special case of the Weibull model.

\item \bold{testIndClogit}. When the data come from a case-control study, the suitable test is via conditional logistic regression. 

\item \bold{testIndMVreg}. In the case of multivariate continuous targets, the suggested test is via a multivariate linear regression. 
The target variable can be compositional data as well. These are positive data, whose vectors sum to 1.  They can sum to any constant, as 
long as it the same, but for convenience reasons we assume that they are normalised to sum to 1.  In this case the additive log-ratio 
transformation (multivariate logit transformation) is applied beforehand. 

\item \bold{testIndSPML}. With a circular target, the projected bivariate normal distribution (Presnell et al., 1998) is used to perform regression. 
}
}

\section{Tests for clustered/longitudinal data}{
\enumerate{

\item \bold{testIndGLMMReg}, \bold{testIndGLMM}, \bold{testIndGLMMPois} & \bold{testIndGLMMLogistic}. In the case of a longitudinal or 
clustered targets (continuous, proportions, binary or counts), the suggested test is via a (generalised) linear mixed model. 
\bold{testIndGLMMCR} stands for mixed effects Cox regression. 

\item \bold{testIndGEEReg}, \bold{testIndGEELogistic}, \bold{testIndGEEPois}, \bold{testIndGEENormLog} and \bold{testIndGEEGamma}. 
In the case of a longitudinal or clustered targets (continuous, proportions, binary, counts, positive, strictly positive), 
the suggested test is via GEE (Generalised Estimating Equations).  
}
}

\section{Wald based tests}{

The available tests for wald.ses and wald.mmpc are listed below. Note, that only continuous predictors are allowed. 

\tabular{llll}{ 
  \bold{Target variable} \tab \bold{Available tests} \tab	\bold{Short explanation} \cr 
  
  Continuous  \tab waldMMReg \tab MM regression \cr  
  Proportions \tab waldMMReg \tab MM regression \cr 
  \tab \tab after logit transformation \cr  
  Proportions \tab waldBeta \tab Beta regression \cr 
  Non negative \tab waldIGreg \tab Inverse Gaussian regression \cr 
  Strictly positive \tab waldGamma \tab Gamma regression \cr 
  Non negative \tab waldNormLog \tab Gaussian regression with log link \cr 
  Successes \& totals \tab testIndBinom \tab Binomial regression \cr
  Discrete \tab waldPois \tab Poisson regression \cr 
  Discrete \tab waldSpeedPois \tab Poisson regression \cr 
  Discrete \tab waldZIP \tab Zero Inflated \cr 
  \tab	\tab  Poisson regression \cr 
  Discrete \tab waldNB \tab Negative binomial regression \cr 
  Factor with two \tab waldLogistic \tab Logistic regression \cr
  levels or binary \tab \tab \cr  
  Factor with more than \tab waldOrdinal \tab Ordinal logistic regression \cr		 
  two levels (ordered) \tab \tab \cr  
  Left censored \tab waldTobit \tab Tobit regression \cr 
  Case-control  \tab Mixed \tab testIndClogit \cr
  \tab  \tab  Conditional logistic regression \cr
  Survival \tab waldCR \tab Cox regression \cr
  Survival \tab waldWR \tab Weibull regression \cr
  Survival \tab waldER \tab Exponential regression \cr
  Survival \tab waldLLR \tab Log-logistic regression \cr
}
}

\section{Permutation based tests}{

The available tests for perm.ses and perm.mmpc are listed below. Note, that only continuous predictors are allowed. 

\tabular{llll}{ 
  \bold{Target variable} \tab \bold{Available tests} \tab	\bold{Short explanation} \cr 

  Continuous  \tab permFisher \tab Pearson correlation \cr  
  Continuous  \tab permMMFisher \tab Robust Pearson correlation \cr  
  Continuous  \tab permDcor \tab Distance correlation \cr  
  Continuous  \tab permReg \tab Linear regression \cr  
  Proportions \tab permReg \tab Linear regression \cr 
  \tab \tab after logit transformation \cr  
  Proportions \tab permBeta \tab Beta regression \cr 
  Non negative \tab permIGreg \tab Inverse Gaussian regression \cr 
  Strictly positive \tab permGamma \tab Gamma regression \cr 
  Non negative \tab permNormLog \tab Gaussian regression with log link \cr 
  Non negative \tab permWR \tab Weibull regression \cr
  Successes \& totals \tab permBinom \tab Binomial regression \cr
  Discrete \tab permPois \tab Poisson regression \cr 
  Discrete \tab permZIP \tab Zero Inflated \cr 
  \tab	\tab  Poisson regression \cr 
  Discrete \tab permNB \tab Negative binomial regression \cr 
  Factor with two \tab permLogistic \tab Binary logistic regression \cr
  levels or binary \tab \tab \cr  
  Factor with more than \tab permMultinom \tab Multinomial logistic regression \cr		 
  two levels (nominal) \tab \tab \cr 
  Factor with more than \tab permOrdinal \tab Ordinal logistic regression \cr		 
  two levels (ordered) \tab \tab \cr  
  Left censored \tab permTobit \tab Tobit regression \cr 
  Survival	 \tab permCR \tab Cox regression \cr
  Survival	 \tab permWR \tab Weibull regression \cr
  Survival	 \tab permER \tab Exponential regression \cr
  Survival	 \tab permLLR \tab Log-logistic regression \cr
}
}

\details{
 These tests can be called by SES, MMPC, wald.mmpc or individually by the user. In all regression cases, there is an option for weights. 
}

\references{
Aitchison J. (1986). The Statistical Analysis of Compositional Data, Chapman & Hall; reprinted in 2003, with additional material, by The Blackburn Press.

Brown P.J. (1994). Measurement, Regression and Calibration. Oxford Science Publications.

Cox D.R. (1972). Regression models and life-tables. J. R. Stat. Soc., 34, 187-220.

Demidenko E. (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.

Draper, N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition.

Fieller E.C. and Pearson E.S. (1961). Tests for rank correlation coefficients: II. Biometrika, 48(1 & 2): 29-40.

Ferrari S.L.P. and Cribari-Neto F. (2004). Beta Regression for Modelling Rates and Proportions. Journal of Applied Statistics, 31(7): 799-815.

Gail, M.H., Jay H.L., and Lawrence V.R. (1981). Likelihood calculations for matched case-control studies and survival studies with tied death times. Biometrika 68(3): 703-707.

Gutenbrunner C., Jureckova J., Koenker R. and Portnoy S. (1993). Tests of Linear Hypothesis based on Regression Rank Scores, Journal of NonParametric Statistics 2, 307-331.

Hoerl A.E. and Kennard R.W. (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1): 55-67.

Joseph M.H. (2011). Negative Binomial Regression. Cambridge University Press, 2nd edition.

Koenker R.W. (2005). Quantile Regression. Cambridge University Press.

Lagani V., Kortas G. and Tsamardinos I. (2013). Biomarker signature identification in "omics" with multiclass outcome. Computational and Structural Biotechnology Journal, 6(7): 1-7.

Lagani V. and Tsamardinos I. (2010). Structure-based variable selection for survival data. Bioinformatics Journal 16(15): 1887-1894.

Lambert D. (1992). Zero-inflated Poisson regression, with an application to defects in manufacturing. Technometrics 34(1)1: 1-14.

Liang K.Y. and Zeger S.L. (1986). Longitudinal data analysis using generalized linear models. Biometrika, 73(1): 13-22.

Mardia K.V., Kent J.T. and Bibby J.M. (1979). Multivariate Analysis. Academic Press, New York, USA. 

Maronna R.D. Yohai M.V. (2006). Robust Statistics, Theory and Methods. Wiley.

McCullagh P. and Nelder J.A. (1989). Generalized linear models. CRC press, USA, 2nd edition.

Paik M.C. (1988). Repeated measurement analysis for nonnormal data in small samples. Communications in Statistics-Simulation and Computation, 17(4): 1155-1171.

Pinheiro J., and D. Bates. Mixed-effects models in S and S-PLUS. Springer Science \& Business Media, 2006.

Prentice R.L. and Zhao L.P. (1991). Estimating equations for parameters in means and covariances of multivariate discrete and continuous responses. Biometrics, 47(3): 825-839.

Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data. Journal of the American Statistical Association, 93(443): 1068-1077.

Scholz, F. W. (2001). Maximum likelihood estimation for type I censored Weibull data including covariates. ISSTECH-96-022, Boeing Information & Support Services.

Smith, R. L. (1991). Weibull regression models for reliability data. Reliability Engineering & System Safety, 34(1), 55-76.

Spirtes P.,  Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, 3nd edition.

Szekely G.J. and Rizzo, M.L. (2014). Partial distance correlation with methods for dissimilarities. The Annals of Statistics, 42(6): 2382--2412.

Szekely G.J. and Rizzo M.L. (2013).  Energy statistics: A class of statistics based on distances. Journal of Statistical Planning and Inference 143(8): 1249--1272.

Therneau T.M., Grambsch P.M. and Pankratz V.S. (2003). Penalized Survival Models and Frailty, Journal of Computational and Graphical Statistics, 12(1):156-175.

Yan J. and Fine J. (2004). Estimating equations for association structures. Statistics in medicine, 23(6): 859-874

Ziegler A., Kastner C., Brunner D. and Blettner M. (2000). Familial associations of lipid proles: A generalised estimating equations approach. Statistics in medicine, 19(24): 3345-3357
}

\author{
Michail Tsagris <mtsagris@uoc.gr>
}

\keyword{ Conditional independence tests }
\keyword{ Regression models }