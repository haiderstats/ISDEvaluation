#File Name: coxPH_KP.R
#Date Created: May 26, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file is used to run the Cox proportional hazards model and then use the  Kalbfleisch-Prentice method to attain survival probabilities.
#Here we will take in a training and a testing set. The function will then train on the training set and return a list containing (1) a 
#matrix of survival curves, where the first column is time values and all following columns are survival probabilities of test subjects and 
#(2) the true death times and event indicator (i.e. time and delta) of the test subjects.

#Input: Survival Dataset post normalization and imputation.
#Output: A list of (1) matrix of survival curves, and (2) the true death times.
############################################################################################################################################
#Library Dependencies
#survival is needed to get survfit and the implimentation of the KP estimator.
library(survival)
#For EN-Cox
library(fastcox)
#For sindex
library(prodlim)

CoxPH_KP = function(training, testing,ElasticNet=F){
  if(ElasticNet){
    timeInd = which(names(training) == "time")
    deltaInd = which(names(training) == "delta")
    alpha = NULL
    lambda = NULL
    bestError = Inf
    for(i in c(0.01,.2,.4,.6,.8,1)){
      model =  cv.cocktail(as.matrix(training[,-c(timeInd, deltaInd)]),training[,timeInd], training[,deltaInd],alpha = i)
      modelBestLambdaIndex = which(model$lambda == model$lambda.min)
      modelError = model$cvm[modelBestLambdaIndex]
      if(modelError < bestError){
        alpha = i
        lambda = model$lambda.min
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

#not considering ties.
KPEstimator = function(lp,lpTime,censorStatus){
  indexToKeep = sindex(sort(lpTime), unique(sort(lpTime)))
  orderLPTime = order(lpTime)
  cumHaz = rev(cumsum(rev(exp(lp[orderLPTime]))))
  alpha = ((1-(exp(lp[orderLPTime])/cumHaz)))^exp(-lp[orderLPTime])
  survivalFunc = cumprod(alpha^censorStatus[orderLPTime])
  return(list(time = lpTime[orderLPTime][indexToKeep],surv = survivalFunc[indexToKeep]))
}

