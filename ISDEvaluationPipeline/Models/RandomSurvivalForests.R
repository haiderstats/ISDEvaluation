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

RSF = function(training, testing,ntree = NULL){
  params = internalCV_RSF(training)
  ntree = params$ntree
  nodesize = params$nodesize
  rsfMod = rfsrc(Surv(time,delta)~., data = training, ntree = ntree, nodesize = nodesize)
  survivalCurves = predict(rsfMod, testing)
  survivalCurvesTrain = predict(rsfMod, training)
  trainingTimes = survivalCurves$time.interest
  if(0 %in% trainingTimes){
    times = trainingTimes
    testProbabilities = t(survivalCurves$survival)
    trainProbabilities = t(survivalCurvesTrain$survival)
  } else{
    times = c(0,trainingTimes)
    testProbabilities = rbind.data.frame(1,t(survivalCurves$survival))
    trainProbabilities = rbind.data.frame(1,t(survivalCurvesTrain$survival))
  }
  curvesToReturn = cbind.data.frame(time = times, testProbabilities)
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  trainingCurvesToReturn = cbind.data.frame(time = times, trainProbabilities)
  return(list(TestCurves = curvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}



internalCV_RSF = function(training, numFolds =5){
  foldedData = createFoldsOfData(training, numFolds)[[2]] 
  
  resultsMatrix = matrix(rep(0,numFolds*9), ncol = 9,nrow =numFolds)
  for(i in 1:numFolds){
    trainingFold = foldedData[[1]][[i]]
    testingFold = foldedData[[2]][[i]]
    
    resultVec = c()
    numParam = ncol(trainingFold) -2
    for(ntree in c(500,1000,2000)){
        for(nodesize in c(1,2,3)){
          rsfMod = rfsrc(Surv(time,delta)~., data = training, ntree = ntree,nodesize = nodesize)
          gc()
          error = predict(rsfMod, testingFold)$err.rate[ntree]
          resultVec = c(resultVec,error)
          }
        }
    resultsMatrix[i,] = resultVec
  }
  meanResults = apply(resultsMatrix, 2, mean)
  print(meanResults)
  bestRes = which.min(meanResults)
  bestNtree = c(500,1000,2000)[ceiling(bestRes/3)]
  bestSize = c(1,2,3)[ifelse(bestRes%%3 ==0, 3,bestRes%%3)]
  return(list(ntree = bestNtree, nodesize = bestSize))
}
