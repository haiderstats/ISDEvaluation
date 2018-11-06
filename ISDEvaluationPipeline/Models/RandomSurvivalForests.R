#### File Information #####################################################################################################################
#File Name: RandomSurvivalForests.R
#Date Created: May 28, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################

#This file is used to run the random survival forests (RSFs) to generate individual survival curves. The implimentation 
#given in the package trains the forest and then drops a test subject down every tree. Each terminal node of every tree is associated with
#a Kaplan Meier (KM) curve. Then these KM curves are averaged pointwise (with no weighting) to produce a final, individual survival curve.

### Functions #############################################################################################################################

## Function 1: RSF(training, testing, numFolds = 5)

#Inputs:
#   training: The training dataset (after normalization and imputation).
#   testing: The testing dataset (after normalization and imputation).
#   numFolds: Number of folds for internal cross validation.

# Output: A list of 4 items:(1) TestCurves - The survival curves for the testing set.
#                           (2) TestData - The censor/death indicator and event time for the testing set. 
#                           (3) TrainData - The censor/death indicator and event time for the training set. 
#                           (4) TrainCurves - The survival curves for the training set.

# Usage: Train and evaluate the RSF model. Select the number of trees and mtry via internal cross validation. Other parameters use their
#        default value given by rsfrc().


## Function 2: internalCV_RSF(training, numFolds)

# Inputs: See RSF().

# Output: The optimal hyper-parameters (ntree, mtry).

# Usage: Perform internal cross validation for selecting ntree and mtry. The selection of these parameters is based on the minimized
#        err.rate given by rsfrc(). (Helper function for RSF).

### Code #############################################################################################################################
#Library Dependencies:
#The random forest implimentation is given in randomForestSRC.
library(randomForestSRC)

RSF = function(training, testing, numFolds = 5){
  params = internalCV_RSF(training, numFolds)
  ntree = params$ntree
  nodesize = params$nodesize
  rsfMod = rfsrc(Surv(time,delta)~., data = training, ntree = ntree, nodesize = nodesize)
  
  survivalCurves = predict(rsfMod, testing)
  survivalCurvesTrain = predict(rsfMod, training)
  trainingTimes = survivalCurves$time.interest
  
  #If 0 wasn't included in the timepoints we would like to manually add it with a survival probability of 1.
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

internalCV_RSF = function(training, numFolds){
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
