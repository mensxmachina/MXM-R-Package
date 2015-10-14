\name{CondIndTests}
\alias{CondIndTests}

\title{ MXM Conditional Independence Tests }
\description{

Currently the \pkg{MXM} package supports numerous tests for different types of target (dependent) and predictor (independent) variables. The target variable can be of continuous, discrete, categorical and of survival type. As for the predictor variables, they can be continuous, categorical or mixed. 

The \bold{testIndFisher} and the \bold{gSquare} tests have two things in common. They do not use a model implicitly (i.e. estimate some beta coefficients), even though there is an underlying assumed one. Secondly they are pure tests of independence (again, with assumptions required).

As for the other tests, they share one thing in common. For all of them, two parametric models must be fit. The null model containing the conditioning set of variables alone and the alternative model containing the conditioning set and the candidate variable. The significance of the new variable is assessed via a log-likelihood ratio test with the appropriate degrees of freedom. All of these tests are summarized in the below table.

\tabular{llll}{ 
  \bold{Target variable} \tab	\bold{Predictor variables} \tab \bold{Available tests} \tab	\bold{Short explanation} \cr 
  
  Continuous \tab Continuous \tab testIndFisher (robust) \tab Partial correlation \cr 
  Continuous \tab Mixed \tab testIndReg (robust) \tab Linear regression \cr 
  Continuous \tab Mixed \tab testIndRQ \tab Median regression \cr  
  Proportions \tab Continuous \tab testIndFisher (robust) \tab Partial correlation \cr 
  \tab \tab after logit transformation \tab \cr
  Proportions \tab Mixed \tab testIndReg(robust) \tab Linear regression \cr 
  \tab \tab after logit transformation \tab \cr  
  Proportions \tab Mixed \tab testIndRQ \tab Median regression \cr 
  \tab \tab after logit transformation \tab \cr  
  Proportions tab Mixed \tab testIndBeta \tab Beta regression \cr  
  Discrete \tab	Mixed \tab testIndPois \tab Poisson regression \cr 
  Numerical \tab \tab testIndNB \tab Negative binomial regression \cr 
  Factor with two \tab Mixed \tab testIndLogistic \tab Binary logistic regression \cr
  levels or binary \tab \tab \tab \cr  
  Factor with more than \tab Mixed \tab testIndLogistic \tab Multinomial logistic regression \cr		 
  two levels (unordered) \tab \tab \tab \cr  
  Factor with more than \tab Mixed \tab testIndLogistic \tab Ordinal logistic regression \cr		 
  two levels (ordered) \tab \tab \tab \cr  
  Categorical \tab Categorical \tab gSquare \tab G-squared test of independence \cr  
  Categorical \tab Categorical \tab testIndLogistic \tab Multinomial logistic regression \cr  
  Categorical \tab Categorical \tab testIndLogistic \tab Ordinal logistic regression \cr  
  Survival	\tab Mixed \tab censIndLR \tab Cox regression 
}

}

\section{Tests}{
  \enumerate{
    \item \bold{testIndFisher}. This is a standard test of independence when both the target and the set of predictor variables are continuous (continuous-continuous). When the joint multivariate normality of all the variables is assumed, we know that if a correlation is zero this means that the two variables are independent. Moving in this spirit, when the partial correlation between the target variable and the new predictor variable conditioning on a set of (predictor) variables is zero, then we have evidence to say they are independent as well. An easy way to calculate the partial correlation between the target and a predictor variable conditioning on some other variables is to regress the both the target and the new variable on the conditioning set. The correlation coefficient of the residuals produced by the two regressions equals the partial correlation coefficient. If the robust option is selected, the two aforementioned regression models are fitted using M estimators (Marona, et al, 2006). 
    
    \item \bold{testIndReg}. In the case of continuous-mixed or continuous-categorical, the suggested test is via a standard linear regression. In this case, two linear regression models are fitted. One with the conditioning set only and one with the conditioning set plus the new variable. The significance of the new variable is assessed via the F test, which calculates the residual sum of squares of the two models. The reason for the F test is because the new variable may be categorical and in this case the t test cannot be used. It makes sense to say, that this test can be used instead of the \bold{testIndFisher}, but it will be slower. If the robust option is selected, the two models are fitted using M estimators (Marona, et al, 2006).
    
    \item \bold{testIndRQ}. An alternative to \bold{testIndReg} for the case of continuous-mixed (or continuous-continuous) variables is the \bold{testIndRQ}. Instead of fitting two linear regression models, which model the expected value, one can choose to model the median of the distribution (Koenker, 2005). The significance of the new variable is assessed via a rank based test calibrated with an F distribution (Gutenbrunner et al., 1993). The reason for this is that we performed simulation studies and saw that this type of test attains the type I error in contrast to the log-likelihood ratio test. The benefit of this regression is that it is robust, in contrast to the classical linear regression. 
    
    \item \bold{testIndBeta}. When the target is proportion (or percentage, i.e., between 0 and 1, not inclusive) the user can fit a regression model assuming a beta distribution. The procedure is the same as in the \bold{testIndReg} case. 
    
    \item \bold{Alternatives to testIndBeta}. Instead of \bold{testIndBeta} the use has the option to choose all the previously to \bold{testIndBeta} mentioned tests by transforming the target variable with the logit transformation. In this way, the support of the target becomes the whole of R^d and then depending on the type of the predictors and whether a robust approach is required or not, there is a variety of alternative to beta regression tests. 
    
    \item \bold{testIndPois}. When the target is discrete, and in specific count data, the default test is via the Poisson regression. The procedure is the same as in all the previously regression model based tests. 
    
    \item \bold{testIndNB}. As an alternative to the Poisson regression, we have included the Negative binomial regression to capture cases of overdispersion. 
    
    \item \bold{testIndLogistic} (Binomial). When the target is categorical with only two outcomes, success or failure for example, then a binary logistic regression is to be used. Whether regression or classification is the task of interest, this method is applicable. The advantage of this over a linear or quadratic discriminant analysis is that it allows for categorical predictor variables as well.
    
    \item \bold{testIndLogistic} (Un-ordered multinomial). If the target has more than two outcomes, but it is of nominal type, there is no ordering of the outcomes, multinomial logistic regression will be employed. Again, this regression is suitable for classification purposes as well and it to allows for categorical predictor variables.
    
    \item \bold{testIndLogistic} (Ordered multinomial). This is a special case of multinomial regression, in which case the outcomes have an ordering, such as \bold{not satisfied}, \bold{neutral}, \bold{satisfied}. The appropriate method is ordinal logistic regression. 
    
    \item \bold{gSquare}. If all variables, both the target and predictors are categorical the default test is the G-square test of independence. It is similar to the chi-squared test of independence, but instead of using the chi-squared metric between the observed and estimated frequencies in contingency tables, the Kullback-Leibler divergence of the observed from the estimated frequencies is used. The asymptotic distribution of the test statistic is a chi-squared distribution on some appropriate degrees of freedom. The target variable can be either ordered or unordered with two or more outcomes. 
    
    \item \bold{Alternatives to gSquare}. An alternative to the \bold{gSquare} test is the \bold{testIndLogistic}. Depending on the nature of the target, binary, un-ordered multinomial or ordered multinomial the appropriate regression model is fitted. 
    
    \item \bold{censIndLR}. Finally, for the case of time-to-event data, a Cox regression model is employed. Again, the log-likelihood ratio test is used to assess the significance of the new variable. 
    
  }
}

\details{
 These tests can be called by SES or individually by the user.
}

\references{
Gutenbrunner C, Jureckova J, Koenker R and S Portnoy (1993). Tests of Linear Hypothesis based on Regression Rank Scores, Journal of NonParametric Statistics 2, 307-331.

Koenker R. W. (2005). Quantile Regression, Cambridge U. Press.

Maronna R. D. Martin V. Yohai (2006). Robust Statistics, Theory and Methods. Wiley.
}

\author{
Michail Tsagris <mtsagris@yahoo.gr>
}

\keyword{MXMCondIndTests}
\keyword{Conditional Independence Tests, Statistical Tests, Tests}
\keyword{Statistical Tests}
\keyword{Tests}