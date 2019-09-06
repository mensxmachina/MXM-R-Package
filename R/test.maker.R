test.maker <- function(test) {
  
  if (test == "testIndReg") {   ## It uMMPC the F test
    test <- testIndReg;

  } else if (test == "testIndFisher") {   ## It uMMPC the F test
    test <- testIndFisher;
    
  } else if (test == "testIndSpearman") {   ## It uMMPC the F test
    test <- testIndSpearman;
    
  } else if (test == "testIndMMFisher") {   ## It uMMPC the F test
    test <- testIndMMFisher;
  
  } else if(test == "testIndBeta") {
    test <- testIndBeta;
  
  } else if(test == "testIndRQ")  {   ## quantile regression
    test <- testIndRQ;
  
  } else if (test == "testIndIGreg") { ## Inverse Gaussian regression
    test <- testIndIGreg;
  
  } else if (test == "testIndMMReg") { ## Inverse Gaussian regression
    test <- testIndMMReg;
   
  } else if (test == "testIndPois") {  ## Poisson regression
    test <- testIndPois;
  
  } else if (test == "testIndNB") {  ## Negative binomial regression
    test <- testIndNB;
  
  } else if (test == "testIndGamma") {  ## Gamma regression
    test <- testIndGamma;
  
  } else if (test == "testIndNormLog") { ## Normal regression with a log link
    test <- testIndNormLog;
  
  } else if (test == "testIndZIP") {  ## Zero inflated Poisson regression
    test <- testIndZIP;
  
  } else if (test == "testIndTobit") { ## Tobit regression
    test <- testIndTobit;
  
  } else if (test == "censIndCR") {
    test <- censIndCR;
  
  } else if (test == "censIndWR") {
    test <- censIndWR;
  
  } else if (test == "censIndER") {
    test <- censIndER;
    
  } else if (test == "censIndLLR") {
    test <- censIndLLR;
  
  } else if (test == "testIndClogit") {
    test <- testIndClogit;
  
  } else if (test == "testIndBinom") {
    test <- testIndBinom;
  
  } else if (test == "testIndLogistic") {
    test <- testIndLogistic;
  
  } else if (test == "testIndMultinom") {
    test <- testIndMultinom;	
  
  } else if (test == "testIndOrdinal") {
    test <- testIndOrdinal;	
  
  } else if (test == "testIndQBinom") {
    test <- testIndQBinom;
  
  } else if (test == "testIndQPois") {
    test <- testIndQPois;
  
  } else if (test == "gSquare") {
    test <- gSquare;
    
  } else if (test == "testIndSPML") {
    test <- testIndSPML;
    
    ######################## 
  
  } else if (test == "testIndGLMMReg") {
    test <- testIndGLMMReg
    
  } else if (test == "testIndLMM") {
    test <- testIndLMM
    
  } else if (test == "testIndGLMMPois") {
    test <- testIndGLMMPois
    
  } else if (test == "testIndGLMMLogistic") {
    test <- testIndGLMMLogistic
    
  } else if (test == "testIndGLMMGamma") {
    test <- testIndGLMMGamma
    
  } else if (test == "testIndGLMMNormlog") {
    test <- testIndGLMMNormLog
    
  } else if (test == "testIndGLMMOrdinal") {
    test <- testIndGLMMOrdinal
    
  } else if (test == "testIndGLMMCR") {
    test <- testIndGLMMCR
    
    ##########################
    
  } else if (test == "testIndGEEReg") {
    test <- testIndGEEReg
    
  } else if (test == "testIndGEELogistic") {
    test <- testIndGEELogistic
    
  } else if (test == "testIndGEEPois") {
    test <- testIndGEEPois
    
  } else if (test == "testIndGEEGamma") {
    test <- testIndGEEGamma
    
  } else if (test == "testIndGEENormLog") {
    test <- testIndGEENormLog
  }    
  
  test
}