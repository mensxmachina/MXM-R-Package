### A generic regression model accepting many regression models
reg.fit <- function(y, dataset, event = NULL, reps = NULL, group = NULL, slopes = FALSE, 
                    reml = FALSE, model = NULL, wei = NULL, xnew = NULL) {
  ## possible models are "gaussian" (default), "binary", "binomial", "multinomial", "poisson",
  ## "ordinal", "Cox", "Weibull", "exponential", "zip", "beta", "median", "negbin", "gamma", "normlog",
  ##  "longitudinal", "grouped", "tobit", "qpois", "qbinom", "MM" or  "clogit".
  ## robust is either TRUE or FALSE
  ## y is the y variable, can be a numerical variable, a matrix, a factor, ordinal factor, percentages, or time to event
  ## dataset is the indendent variable(s). It can be a vector, a matrix or a dataframe with continuous only variables, 
  ## a data frame with mixed or only categorical variables
  ## event is NULL unless you have time to event data (survival regression)
  ## reps is NULL unless you have time measurements (longitudinal data)
  ## group is NULL unless you have grouped (or clustered) data or longitudinal data (needs reps)
  ## slopes is for the longitudinal data only, TRUE or FALSE 
  ## reml is for mixed models. If TRUE, REML will be used, otherwise ML will be used
  x <- as.data.frame(dataset)  ## just in case
  if ( is.null( colnames(x) ) )  colnames(x) <- paste("X", 1:ncol(x), sep = "")
  la <- length( unique(y) )

  if ( is.null(model) ) {
    ## linear regression 
    if ( sum( class(y) == "numeric" ) == 1 & is.null(event) & is.null(reps)  &  is.null(group) )  model <- "gaussian"  
    ## multivariate data
    if ( sum( class(y) == "matrix" ) == 1 ) { 
      if ( min(y) > 0 &  Rfast::Var(Rfast::rowsums(y) == 0 ) )  y <- log(y[, -1] / y[, 1])  ## compositional data
      model <- "gaussian"
    }
    ## surival data
    if ( !is.null(event) ) {
      y <- survival::Surv(time = y, event = event)
      model <- "Cox"
    }
    ## longitudinal data
    if ( !is.null(reps) & !is.null(group) )  model <- "longitudinal"
    ## grouped data
    if ( is.null(reps) & !is.null(group) )  model <- "grouped"
    ## binary data
    if ( la == 2 )  model <- "binary"   
    ## ordinal, multinomial or perhaps binary data
    if ( is.factor(y) ) {
      if ( !is.ordered(y) ) {
        if ( la == 2 ) {
          y <- as.vector(y)
          model <- "binary"
        } else  model <- "multinomial"
      } else {
        if ( la == 2 ) {
          y <- as.vector(y)
          model <- "binary"
        } else  model <- "ordinal"    
      }
    }
    ## count data
    if ( sum( is.vector(y) ) == 1 ) {
      if ( ( sum( floor(y) - y ) == 0  &  la > 2 ) )  model <- "poisson"
    }
  
  }  
  ##### model checking
 
     ## univariate gaussian model
  if ( model == "gaussian"  &  is.vector(y) ) {
     mod <- lm(y ~ ., data = x, weights = wei )

  } else if ( model == "MM" ) {
    mod <- MASS::rlm(y ~., data = x, maxit = 2000, method = "MM")
    
  } else if ( model == "gaussian"  &  is.matrix(y) ) {
    mod <- lm(y ~ ., data = x )
 
  } else if ( model == "binomial" &  is.matrix(y) ) {
    mod <- glm(y[, 1] / y[, 2] ~ ., data = x, weights = y[, 2], family = binomial )
    
    ## median (quantile) regression 
  } else if ( model == "median" ) {
    mod <- quantreg::rq(y ~ ., data = x, weights = wei ) 

  } else if ( model == "gamma" ) {
    mod <- glm(y ~ ., data = x, weights = wei, family = Gamma(log) ) 
    
  } else if ( model == "normlog" ) {
    mod <- glm(y ~ ., data = x, weights = wei, family = gaussian(log) ) 
    
  } else if ( model == "tobit" ) {
    mod <- survival::survreg(y ~ ., data = x, weights = wei, dist = "gaussian" )
    
  } else if ( model == "binary" ) {
    mod <- glm(y ~ ., data = x, binomial, weights = wei)

  } else if ( model == "multinomial" ) {
    mod <- nnet::multinom(y ~ ., data = x, trace = FALSE, weights = wei)

  } else if ( model == "ordinal" ) {
    mod <- ordinal::clm(y ~ ., data = x, weights = wei )

  } else if ( model == "poisson" ) {
    mod <- glm(y ~ ., data = x, poisson, weights = wei)

  } else if ( model == "qpois" ) {
    mod <- glm(y ~ ., data = x, quasipoisson, weights = wei)
    
  } else if ( model == "qbinom" ) {
    mod <- glm(y ~ ., data = x, quasibinomial, weights = wei)
    
  } else if ( model == "negbin" ) {
    mod <- MASS::glm.nb(y ~ ., data = x, weights = wei )

  } else if ( model == "zip" ) {
    mod <- zip.mod(y, x, wei = wei)

  } else if ( model == "beta" ) {
    mod <- beta.mod(y, x, wei = wei )
  
  } else if ( model == "cox" ) {
    mod <- survival::coxph(y ~ ., data = x, weights = wei )

  } else if ( model == "weibull" ) {
    mod <- survival::survreg(y ~ ., data = x, weights = wei )

  } else if ( model == "exponential" ) {
    mod <- survival::survreg(y ~ ., data = x, weights = wei, dist = "exponential" )

  } else if ( model == "clogit" ) {
    case <- y[, 1]  ## case control, 0 is the control 
    id <- y[, 2] #the patient id
    mod <- survival::clogit( case ~ . + strata(id), data = x)
                
  } else if ( model == "longitudinal" ) {
    if ( is.ordered(y) ) {
      if (slopes) {
        mod <- ordinal::clmm( y ~ . -group + (reps|group), data =  cbind(reps, x), weights = wei ) 
      } else  mod <- ordinal::clmm( y ~ . -group  + (1|group), data = x, weights = wei )

    } else if ( la > 2  &  sum( round(y) - y ) != 0  ) {
      if ( slopes ) {
        mod <- lme4::lmer( y ~ . -group + (reps|group), REML = reml, data = cbind(reps, x), weights = wei ) 
      } else  mod <- lme4::lmer( y ~ . -group + (1|group), REML = reml, data = x, weights = wei )
	  
    } else if ( la > 2  &  sum( round(y) - y ) == 0 ) {
      if ( slopes ) {
        mod <- lme4::glmer( y ~ . - group + (reps|group), REML = reml, family = poisson,  data = cbind(reps, x), weights = wei ) 
      } else  mod <- lme4::glmer( y ~ . -group  + (1|group), REML = reml, family = poisson, data = x, weights = wei )
	  
    } else  if ( la == 2 ) {
      y <- as.vector(y)   
      if ( slopes ) {
        mod <- lme4::glmer( y ~ . - group + (reps|group), REML = reml, family = binomial, data = cbind(reps, x), weights = wei ) 
      } else  mod <- lme4::glmer( y ~ . -group + (1|group), REML = reml, family = binomial, data = cbind(reps, x), weights = wei )
    }

  } else if ( model == "grouped" ) {
    if ( is.ordered(y) ) {
      mod <- ordinal::clmm( y ~ . -group  + (1|group), data = x, weights = wei )
    } else if ( sum( round(y) - y ) != 0  &  la > 2 ) {
      mod <- lme4::lmer( y ~ . -group + (1|group), REML = reml, data = x, weights = wei )
    } else if ( sum( round(y) - y ) == 0  &  la > 2) {
      mod <- lme4::glmer( y ~ . -group + (1|group), REML = reml, family = poisson, data = x, weights = wei )
    } else if ( la == 2 ) {
      mod <- lme4::glmer( y ~ . -group + (1|group), REML = reml, family = binomial, data = x, weights = wei )
    }
  }
  
  pred <- NULL
  if ( !is.null(xnew) )  {
    xnew <- as.data.frame(xnew)  
    colnames(xnew) <- colnames(x)
    pred <- predict(mod, xnew) 	
    result <- list(mod = mod, pred = pred)
  } else result <- list(mod = mod, pred = pred) 

  result 

}













    
      

