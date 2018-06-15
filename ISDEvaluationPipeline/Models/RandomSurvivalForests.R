#File Name: RandomSurvivalForests.R
#Date Created: May 28, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file is used to run the random survival forests (RSFs) to generate individual survival curves. Note that RSFs already have an 
#implimentaiton of individual survival curves so we should use that as opposed to making up our own scheme. The implimentation 
#given in the package trains the forest and then drops a test subject down every tree. Each terminal node of every tree is associated with
#a Kaplan Meier (KM) curve. Then these KM curves are averaged pointwise (with no weighting) to produce a final, individual survival curve.
#Here we will take in a training and a testing set. The function will then train on the training set and return a list containing (1) a 
#matrix of survival curves, where the first column is time values and all following columns are survival probabilities of test subjects and 
#(2) the true death times and event indicator (i.e. time and delta) of the test subjects.

#Input 1: Survival Dataset post normalization and imputation.
#Input 2: Number of Trees
#Output: A list of (1) matrix of survival curves, and (2) the true death times.
############################################################################################################################################
#Library Dependencies
#The random forest implimentation is given in randomForestSRC.
library(randomForestSRC)

RSF = function(training, testing, ntree){
  rsfMod = rfsrc(Surv(time,delta)~., data = training, ntree = ntree)
  survivalCurves = predict(rsfMod, testing)
  trainingTimes = survivalCurves$time.interest
  if(0 %in% trainingTimes){
    times = trainingTimes
    probabilities = t(survivalCurves$survival)
  } else{
    times = c(0,trainingTimes)
    probabilities = rbind.data.frame(1,t(survivalCurves$survival))
  }
  curvesToReturn = cbind.data.frame(time = times, probabilities)
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  return(list(curvesToReturn, timesAndCensTest,timesAndCensTrain))  
}











