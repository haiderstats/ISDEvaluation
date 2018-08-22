#File Name: OneCalibration.R
#Date Created: May 29th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file was created to implement one calibration as an evaluation metric for individual survival curves. There are three different
#implementations of one calibration: Uncensored, Fractional, and BucketKM.
#Uncensored: This method places all indiviuals into decile buckets and then removes all censored individuals.
#DN: This method is known as the D'Agostino-Nam method and builds a KM curve within each bucket and uses the KM estimate of t* as the
#value for the observed number of events. For numBuckets <= 15 this follows a chi.square statistic with Numbucket-1 degrees of freedom, 
#For numBuckets > 15 it follows chi-squared with numBuckets -2 degrees of freedom.
#Input 1: A list of (1) a matrix of survival curves, and (2) the true death times and event/censoring indicator (delta =1 implies death/event).
#Output: The desired L-measure value.
##############################################################################################################################################
#Library Dependencies
#We use this for the Surv function.
library(survival)

#We use this for group_by, summarise, and mutate.
library(dplyr)

#Helper Functions: predictProbabilityFromCurve(survivalCurve,predictedTimes, timeToPredict)
source("Evaluations/EvaluationHelperFunctions.R")

OneCalibration = function(survMod, timeOfInterest = c(), type = "DN", numBuckets = 10){
  #Being passed an empty model.
  if(is.null(survMod)) return(NULL)
  #Being passed a model that failed.
  suppressWarnings(if(is.na(survMod[[1]])) return(NULL))
  if(!type %in% c("Uncensored","DN"))
    stop("Please enter one of 'Uncensored','DN' for type.")
  predictedTimes = survMod[[1]]$time
  survivalCurves = survMod[[1]][-1] #removes predicted times.
  trueDeathTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  trainingDeathTimes = survMod[[3]]$time
  trainingCensorStatus = survMod[[3]]$delta
  #If time is null, we will select the median from the KM curve generated from training instances.
  if(is.null(timeOfInterest)){
    timeOfInterest = unname(quantile(trueDeathTimes,c(.1,.25,.5,.75,.9)))
  }
  predictions = unlist(lapply(seq_along(trueDeathTimes),
                              function(index) predictProbabilityFromCurve(survivalCurves[,index],
                                                                          predictedTimes,
                                                                          timeOfInterest)))
  pval = c()
  for(times in 1:length(timeOfInterest)){
    #If we have multiple times of interest we need to pull out the predictions for a specific time.
    predictionsSpecified = predictions[seq(times,length(timeOfInterest)*length(trueDeathTimes),length(timeOfInterest))]
    pval = c(pval,binItUp(trueDeathTimes,censorStatus, predictionsSpecified, type, numBuckets,timeOfInterest[times]))
  }
  return(pval)
}

OneCalibrationCumulative = function(listOfSurvivalModels, timeOfInterest = c(), type = "DN", numBuckets = 10){
  if(length(listOfSurvivalModels) ==0) return(NULL)
  suppressWarnings(if(any(unlist(lapply(listOfSurvivalModels, is.na)))) return(NA))
  if(!type %in% c("Uncensored","DN"))
    stop("Please enter one of 'Uncensored','DN' for type.")
  
  predictedTimes = lapply(seq_along(listOfSurvivalModels), function(x) listOfSurvivalModels[[x]][[1]]$time)
  survivalCurves = lapply(seq_along(listOfSurvivalModels), function(x) listOfSurvivalModels[[x]][[1]][-1])
  trueDeathTimes = lapply(seq_along(listOfSurvivalModels), function(x) listOfSurvivalModels[[x]][[2]]$time)
  censorStatus   = lapply(seq_along(listOfSurvivalModels), function(x) listOfSurvivalModels[[x]][[2]]$delta)
  #If time is null, we will select the median from the training instances (of those who were uncensored.)
  if(is.null(timeOfInterest)){
    allTimes = unlist(trueDeathTimes)
    timeOfInterest = unname(quantile(allTimes,c(.1,.25,.5,.75,.9)))
  }
  pval = c()
  for(times in 1:length(timeOfInterest)){
    predictions = unlist(lapply(seq_along(listOfSurvivalModels), function(model) lapply(seq_along(trueDeathTimes[[model]]),
                                                                                        function(index) predictProbabilityFromCurve(survivalCurves[[model]][,index],
                                                                                                                                    predictedTimes[[model]],
                                                                                                                                    timeOfInterest[times]))))
    pval = c(pval,binItUp(unlist(trueDeathTimes),unlist(censorStatus), predictions, type, numBuckets,timeOfInterest[times]))
  }
  return(pval)
}


binItUp = function(trueDeathTimes,censorStatus, predictions, type, numBuckets,timeOfInterest){
  #We need to divide the number of predictions into as equal size buckets as possible. The formula for this can be found here:
  #https://math.stackexchange.com/questions/199690/divide-n-items-into-m-groups-with-as-near-equal-size-as-possible
  numberOfBucketsWithExtra = length(predictions) %% numBuckets
  numberOfBucketsWithoutExtra = numBuckets - numberOfBucketsWithExtra
  #We want buckets to randomly be assigned their respective number of predictions, but be consistent with the same input.
  #E.g. for 3 buckets with size 5 and 2 buckets of size 4 we don't necasarrily want 5,5,5,4,4. We want them to be random. We use sample().
  set.seed(42)
  bucketSizes = sample(c(rep(floor(length(predictions)/numBuckets)+1, numberOfBucketsWithExtra),
                         rep(floor(length(predictions)/numBuckets), numberOfBucketsWithoutExtra)))
  bucketLabel = rep(1:numBuckets, times = bucketSizes)
  
  #We are going to put the survival probabilities into deciles so we get the order of indexs here.
  orderPredictions = order(predictions)
  orderedPredictionFrame = data.frame(prob = predictions[orderPredictions],
                                      label = bucketLabel,
                                      delta = censorStatus[orderPredictions],
                                      time = trueDeathTimes[orderPredictions])
  timeFrame = data.frame(time = trueDeathTimes[orderPredictions],
                         label = bucketLabel,
                         delta = censorStatus[orderPredictions])
  pval = switch(type,
                Uncensored = {
                  #Read the following code as taking our predicted probabilities, placing them into buckets, 
                  #taking only uncensored individuals at the time of interest and then taking their mean. All other following piping operations should
                  #be read in a similar fashion.
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
                #For a good discusion of this method, also known as the D'Agostino-Nam tranlation/method see
                #Section 2.2 of this thesis paper: Hosmer-Lemeshow goodness-of-fit test: Translations to the Cox Proportional Hazards Model
                DN = {
                  bucketSurvival = orderedPredictionFrame  %>%
                    group_by(label) %>%
                    summarise(survivalPrediction = mean(prob))
                  numDied = timeFrame  %>%
                    group_by(label) %>%
                    mutate(slope = (1-min(prodlim(Surv(time, delta)~1)$surv))/(0-max(time))) %>%
                    summarise(deadCount = 1-ifelse(is.na(predict(prodlim(Surv(time, delta)~1),timeOfInterest)),max(0,1+slope*timeOfInterest),
                                                     predict(prodlim(Surv(time, delta)~1),timeOfInterest)))
                  observed = numDied$deadCount
                  expected = 1-bucketSurvival$survivalPrediction
                  HLStat = (bucketSizes[bucketSizes!=0]*(observed-expected)^2)/((1-bucketSurvival$survivalPrediction)*bucketSurvival$survivalPrediction)
                  HLStat = sum(ifelse(is.nan(HLStat),0,HLStat))
                  #See comments in file header for reasoning of degree of freedom choice.
                  DoF = ifelse(numBuckets > 15, numBuckets -2, numBuckets-1)
                  pval = 1-pchisq(HLStat, DoF)
                })
  return(pval)
}








