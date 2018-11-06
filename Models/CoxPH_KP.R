#### File Information #####################################################################################################################
#File Name: coxPH_KP.R
#Date Created: May 26, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#This file is used to run the Cox proportional hazards model and then use the  Kalbfleisch-Prentice method to attain survival probabilities.
#Additionally, the elastic net cox is also implemented here.

### Functions #############################################################################################################################

## Function 1: CoxPH_KP(training, testing,ElasticNet =F, numFolds =5)

#Inputs:
#   training:   The training dataset (after normalization and imputation).
#   testing:    The testing dataset (after normalization and imputation).
#   ElasticNet: Boolean indicating to use or not use elastic net cox.
#   numFolds:   Number of folds for internal cross validation (elastic-net cox only).

# Output: A list of 4 items:(1) TestCurves - The survival curves for the testing set.
#                           (2) TestData - The censor/death indicator and event time for the testing set. 
#                           (3) TrainData - The censor/death indicator and event time for the training set. 
#                           (4) TrainCurves - The survival curves for the training set.

# Usage: Train and evaluate the Cox-KP/CoxEN-KP model.

## Function 2: KPEstimator(lp,lpTime,censorStatus){

# Inputs:
#   lp:           The linear predictors from the coxEN model (for training data).
#   lpTime:       The training times.
#   censorStatus: A vector indicating the censor status of the patients.

# Output: The baseline survival function via the kalbfleisch-prentice estimator.

# Usage: This will build the baseline survival function for patients (because cox only considers hazards we need to use a method
#         of estimating the baseline survival function). (Helpher function for CoxPH_KP).

### Code #############################################################################################################################
#Library Dependencies:
#survival is needed to get survfit and the implimentation of the KP estimator.
library(survival)
#For EN-Cox
library(fastcox)
#For sindex
library(prodlim)

CoxPH_KP = function(training, testing,ElasticNet=F, numFolds = 5){
  if(ElasticNet){
    timeInd = which(names(training) == "time")
    deltaInd = which(names(training) == "delta")
    alpha = NULL
    lambda = NULL
    bestError = Inf
    #Try 6 values for alpha and then do repeated (10 times) k-fold cross validation for selecting best hyper parameters.
    for(a in c(0.01,.2,.4,.6,.8,1)){
      errors = NULL
      #Taken from https://stats.stackexchange.com/questions/97777/variablity-in-cv-glmnet-results
      for(i in 1:10){
        model =  cv.cocktail(as.matrix(training[,-c(timeInd, deltaInd)]),training[,timeInd], training[,deltaInd],
                             alpha = a,nfolds = numFolds)
        errors = cbind(errors, model$cvm)
      }
      rownames(errors) <- model$lambda
      avgErrors = rowMeans(errors)
      lambda.min <- as.numeric(names(which.min(avgErrors)))
      modelError = min(avgErrors)
      if(modelError < bestError){
        alpha = a
        lambda = lambda.min
        bestError = modelError
      }
    }
    coxModel = cocktail(as.matrix(training[,-c(timeInd, deltaInd)]),training[,timeInd], training[,deltaInd],alpha = alpha,lambda = lambda)
    linearPredictionsTraining = predict(coxModel,as.matrix(training[,-c(timeInd, deltaInd)]),type = "link")
    linearPredictionsTesting = predict(coxModel,as.matrix(testing[,-c(timeInd, deltaInd)]),type = "link")
    survivalEstimate = KPEstimator(linearPredictionsTraining, training$time,training$delta)
    survCurvs = t(sapply(survivalEstimate[[2]], function(x) x^exp(linearPredictionsTesting)))
    survCurvsTraining = t(sapply(survivalEstimate[[2]], function(x) x^exp(linearPredictionsTraining)))
    
    survivalCurves = list(time = survivalEstimate[[1]], surv = survCurvs)
    survivalCurvesTrain = list(time = survivalEstimate[[1]], surv = survCurvsTraining)
    
  }
  else{
    #Sometimes the coxPH-KP failed to converge so we catch that here.
    tryCatch({
      coxModel = coxph(Surv(time,delta)~., data = training,singular.ok = T)
      survivalCurves = survfit(coxModel, testing, type = "kalbfleisch-prentice")
      survivalCurvesTrain = survfit(coxModel, training, type = "kalbfleisch-prentice")
      
    },
    error = function(e) {
      message(e)
      warning("Cox-PH failed to converge.")
    })
    if(!exists("coxModel") | !exists("survivalCurves")){
      return(NA)
    }
  }
  #If 0 wasnt included in the timepoints we would like to manually add it with a survival probability of 1.
  if(0 %in% survivalCurves$time){
    timePoints = survivalCurves$time
    probabilities = survivalCurves$surv
    
    probabilitiesTrain = survivalCurvesTrain$surv
  } else{
    timePoints = c(0,survivalCurves$time)
    probabilities = rbind(1,survivalCurves$surv)
    
    probabilitiesTrain = rbind(1,survivalCurvesTrain$surv)
  }
  curvesToReturn = cbind.data.frame(time = timePoints, probabilities)
  trainingCurvesToReturn = cbind.data.frame(time = timePoints, probabilitiesTrain)
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  return(list(TestCurves = curvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
  
}

#Currently this estimator uses parameters which assume new ties. Future implementations should adjust the code for this (though estimates)
#are nearly identical when using the ties-adjusted model.
KPEstimator = function(lp,lpTime,censorStatus){
  indexToKeep = sindex(sort(lpTime), unique(sort(lpTime)))
  orderLPTime = order(lpTime)
  cumHaz = rev(cumsum(rev(exp(lp[orderLPTime]))))
  alpha = ((1-(exp(lp[orderLPTime])/cumHaz)))^exp(-lp[orderLPTime])
  survivalFunc = cumprod(alpha^censorStatus[orderLPTime])
  return(list(time = lpTime[orderLPTime][indexToKeep],surv = survivalFunc[indexToKeep]))
}

