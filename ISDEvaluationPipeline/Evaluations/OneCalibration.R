#File Name: OneCalibration.R
#Date Created: May 29th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file was created to implement one calibration as an evaluation metric for individual survival curves.
#Input 1: A list of (1) a matrix of survival curves, and (2) the true death times and event/censoring indicator (delta =1 implies death/event).
#Output: The desired L-measure value.
##############################################################################################################################################
#Library Dependencies
#We use this for the sindex function.
library(prodlim)
#We use this for ldply, a combiner of lists.
library(dplyr)
#We use this for hoslem.test
library(ResourceSelection)
#Helper Functions: predictProbabilityFromCurve(survivalCurve,predictedTimes, timeToPredict)
source("Evaluations/EvaluationHelperFunctions.R")

#type %in% c("Uncensored","Fractional","BucketKM")
OneCalibration = function(survMod, timeOfInterest = NULL, type, numBuckets){
  predictedTimes = survMod[[1]][,1]
  survivalCurves = survMod[[1]][-1]
  trueDeathTimes = survMod[[2]][,1]
  censorStatus = survMod[[2]][,2]
  #TODO: I shouldn't use the test data...
  if(is.null(timeOfInterest))
    timeOfInterest = median(trueDeathTimes)
  predictions = unlist(lapply(seq_along(trueDeathTimes),
                              function(index) predictProbabilityFromCurve(survivalCurves[,index],
                                                                          predictedTimes,
                                                                          timeOfInterest)))
  orderPredictions = order(predictions)
  #https://math.stackexchange.com/questions/199690/divide-n-items-into-m-groups-with-as-near-equal-size-as-possible
  numberOfBucketsWithExtra = length(predictions) %% numBuckets
  numberOfBucketsWithoutExtra = numBuckets - numberOfBucketsWithExtra
  bucketSizes = sample(c(rep(floor(length(predictions)/numBuckets)+1, numberOfBucketsWithExtra),
                         rep(floor(length(predictions)/numBuckets), numberOfBucketsWithoutExtra)))
  bucketLabel = rep(1:numBuckets, times = bucketSizes)
  orderedPredictionFrame = data.frame(prob = predictions[orderPredictions], label = bucketLabel, delta = censorStatus[orderPredictions])
  timeFrame = data.frame(time = trueDeathTimes[orderPredictions], label = bucketLabel, delta = censorStatus[orderPredictions])
  if(type == "Uncensored"){
    bucketSurvival = orderedPredictionFrame  %>% group_by(label) %>% filter(delta ==1) %>% summarise(survivalPrediction = mean(prob))
    numDied = timeFrame  %>% group_by(label) %>% filter(delta ==1) %>% summarise(deadCount = length(which(time < timeOfInterest)))
    observed = numDied$deadCount
    expected = bucketSizes*(1-bucketSurvival$survivalPrediction)
    HLStat = sum((observed-expected)^2/(bucketSizes*(1-bucketSurvival$survivalPrediction)*bucketSurvival$survivalPrediction))
    pval = 1-pchisq(HLStat, numBuckets-2)
    return(pval)
  }
}











