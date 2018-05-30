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
  predictedTimes = survMod[[1]]$time
  survivalCurves = survMod[[1]][-1] #removes predicted times.
  trueDeathTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  #If time is null, we will select the median from the training instances (of those who were uncensored.)
  if(is.null(timeOfInterest))
    timeOfInterest = median(survMod[[3]][survMod[[3]]$delta]$time)
  predictions = unlist(lapply(seq_along(trueDeathTimes),
                              function(index) predictProbabilityFromCurve(survivalCurves[,index],
                                                                          predictedTimes,
                                                                          timeOfInterest)))
  orderPredictions = order(predictions)
  #https://math.stackexchange.com/questions/199690/divide-n-items-into-m-groups-with-as-near-equal-size-as-possible
  numberOfBucketsWithExtra = length(predictions) %% numBuckets
  numberOfBucketsWithoutExtra = numBuckets - numberOfBucketsWithExtra
  set.seed(42)
  bucketSizes = sample(c(rep(floor(length(predictions)/numBuckets)+1, numberOfBucketsWithExtra),
                         rep(floor(length(predictions)/numBuckets), numberOfBucketsWithoutExtra)))
  bucketLabel = rep(1:numBuckets, times = bucketSizes)
  orderedPredictionFrame = data.frame(prob = predictions[orderPredictions],
                                      label = bucketLabel,
                                      delta = censorStatus[orderPredictions],
                                      time = trueDeathTimes[orderPredictions])
  timeFrame = data.frame(time = trueDeathTimes[orderPredictions],
                         label = bucketLabel,
                         delta = censorStatus[orderPredictions])
  pval = switch(type,
                Uncensored = {
                  bucketSurvival = orderedPredictionFrame  %>%
                    group_by(label) %>% 
                    filter(delta ==1 | (delta ==0 & time >= timeOfInterest)) %>%
                    summarise(survivalPrediction = mean(prob))
                  numDied = timeFrame  %>%
                    group_by(label) %>%
                    filter(delta ==1| (delta ==0 & time >= timeOfInterest)) %>%
                    summarise(deadCount = length(which(time < timeOfInterest)))
                  observed = numDied$deadCount
                  expected = bucketSizes*(1-bucketSurvival$survivalPrediction)
                  HLStat = sum((observed-expected)^2/(bucketSizes*(1-bucketSurvival$survivalPrediction)*bucketSurvival$survivalPrediction))
                  pval = 1-pchisq(HLStat, numBuckets-2)
                },
                Fractional = {
                  KM = prodlim(Surv(time,delta)~1, data = survMod[[3]])
                  bucketSurvival = orderedPredictionFrame  %>% 
                    group_by(label) %>%
                    summarise(survivalPrediction = mean(prob))
                  numDied = timeFrame  %>%
                    #Case 1: Uncensored and died before time of interest, cost = 1.
                    #Case 2: Censored before time of interest, use fractional person method for cost. 
                    #Case 3: Either uncensored but died after time of interest or censored at or after time of interest, cost = 0.
                    mutate(contribution = ifelse(delta == 1 & time <= timeOfInterest, 1,
                                                 ifelse(delta == 0 & time < timeOfInterest,
                                                        predict(KM,timeOfInterest)/predict(KM,time,level.chaos=2), 0))) %>%
                    group_by(label) %>%
                    summarise(deadCount = sum(contribution))
                  observed = numDied$deadCount
                  expected = bucketSizes*(1-bucketSurvival$survivalPrediction)
                  HLStat = sum((observed-expected)^2/(bucketSizes*(1-bucketSurvival$survivalPrediction)*bucketSurvival$survivalPrediction))
                  pval = 1-pchisq(HLStat, numBuckets-2)
                },
                BucketKM = {
                  bucketSurvival = orderedPredictionFrame  %>%
                    group_by(label) %>%
                    summarise(survivalPrediction = mean(prob))
                  numDied = timeFrame  %>%
                    group_by(label) %>%
                    summarise(deadCount = 1 - predict(prodlim(Surv(time, delta)~1),timeOfInterest))
                  observed = numDied$deadCount
                  expected = 1-bucketSurvival$survivalPrediction
                  HLStat = sum((bucketSizes*(observed-expected)^2)/((1-bucketSurvival$survivalPrediction)*bucketSurvival$survivalPrediction))
                  pval = 1-pchisq(HLStat, numBuckets-1)
                })
  return(pval)
}











