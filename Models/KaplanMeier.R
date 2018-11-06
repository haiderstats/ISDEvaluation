#### File Information #####################################################################################################################
#File Name: KaplanMeier.R
#Date Created: May 28, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#This file is used to run the Kaplan Meier to attain individual survival probabilities. While typically people use the implimentation given
#in the survival package, we instead use prodlim as prodlim has a predict function built in whereas survival does not.

### Functions #############################################################################################################################
## Function 1: KM(training, testing)

#Inputs:
#   training: The training dataset (after normalization and imputation).
#   testing: The testing dataset (after normalization and imputation).

# Output: A list of 4 items:(1) TestCurves - The survival curves for the testing set.
#                           (2) TestData - The censor/death indicator and event time for the testing set. 
#                           (3) TrainData - The censor/death indicator and event time for the training set. 
#                           (4) TrainCurves - The survival curves for the training set.

# Usage: Train and evaluate the KM model.

### Code #############################################################################################################################
#Library Dependencies:
#prodlim has a Kaplan Meier implementation that includes a predict function.
library(prodlim)

KM = function(training, testing){
  kmMod = prodlim(Surv(time,delta)~1, data = training)
  trainingTimes = c(sort(unique(training$time)))
  #If 0 wasnt included in the timepoints we would like to manually add it.
  if(0 %in% trainingTimes){
    times = trainingTimes
  } else{
    times = c(0,trainingTimes)
  }
  probabilities = predict(kmMod,times)
  curvesToReturn = cbind.data.frame(time = times,matrix(rep(probabilities,nrow(testing)),ncol = nrow(testing)))
  trainingCurvesToReturn = cbind.data.frame(time = times,matrix(rep(probabilities,nrow(training)),ncol = nrow(training)))
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  return(list(TestCurves = curvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn)) 
}
