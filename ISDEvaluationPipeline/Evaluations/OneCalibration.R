#File Name: OneCalibration.R
#Date Created: May 29th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file was created to implement one calibration as an evaluation metric for individual survival curves. There are three different
#implementations of one calibration: Uncensored, Fractional, and BucketKM.
#Uncensored: This method places all indiviuals into decile buckets and then removes all censored individuals.
#Fractional: This method builds a KM curve from the training data and then uses the KM predictions of the time of interest for observed number
#of events for censored patients. Specifically we have that P(T>t* | T>c) = P(T>t* ^ T >c)/P(T>c) = P(T>t*)/P(T>c). P is generated using KM.
#BucketKM: This method is known as the D'Agostino-Nam method and builds a KM curve within each bucket and uses the KM estimate of t* as the
#value for the observed number of events. For numBuckets <= 15 this follows a chi.square statistic with Numbucket-1 degrees of freedom, 
#For numBuckets > 15 it follows chi-squared with numBuckets -2 degrees of freedom.
#Input 1: A list of (1) a matrix of survival curves, and (2) the true death times and event/censoring indicator (delta =1 implies death/event).
#Output: The desired L-measure value.
##############################################################################################################################################
#Library Dependencies
#We use this for the Surv function.
library(survival)
#We use this for the sindex function.
library(prodlim)
#We use this for group_by, summarise, and mutate.
library(dplyr)

#Helper Functions: predictProbabilityFromCurve(survivalCurve,predictedTimes, timeToPredict)
source("Evaluations/EvaluationHelperFunctions.R")

OneCalibration = function(survMod, timeOfInterest = NULL, type = "BucketKM", numBuckets = 10){
  #Being passed an empty model.
  if(is.null(survMod)) return(NULL)
  #Being passed a model that failed.
  suppressWarnings(if(is.na(survMod[[1]])) return(NA))
  if(!type %in% c("Uncensored","Fractional","BucketKM"))
    stop("Please enter one of 'Uncensored','Fractional','BucketKM' for type.")
  predictedTimes = survMod[[1]]$time
  survivalCurves = survMod[[1]][-1] #removes predicted times.
  trueDeathTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  trainingDeathTimes = survMod[[3]]$time
  trainingCensorStatus = survMod[[3]]$delta
  #If time is null, we will select the median from the KM curve generated from training instances.
  if(is.null(timeOfInterest)){
    KMCurve = prodlim(Surv(trainingDeathTimes, trainingCensorStatus)~1)
    timeOfInterest = tryCatch({
      quantile(KMCurve,.5)$quantiles.survival$quantile
    },
    #Catch the case where there is too many censored patients to generate a KM Curve to the median line.
    #In this case we use our linear extension to find the median.
    error = function(e){
      slope = (1-min(KMCurve$surv))/(0 - max(KMCurve$time))
      timeOfInterest = (-0.5)/slope
    })
  }
  #Additionally sometimes the quantile function will return NA instead of an error if it still able to produce a lower bound. Here we 
  #do the same thing as above when we catch the error.
  if(is.na(timeOfInterest)){
    slope = (1-min(KMCurve$surv))/(0 - max(KMCurve$time))
    timeOfInterest = (-0.5)/slope
  }
  predictions = unlist(lapply(seq_along(trueDeathTimes),
                              function(index) predictProbabilityFromCurve(survivalCurves[,index],
                                                                          predictedTimes,
                                                                          timeOfInterest)))
  pval = binItUp(trueDeathTimes, censorStatus, predictions, type, numBuckets,timeOfInterest)
  return(pval)
}

OneCalibrationCumulative = function(listOfSurvivalModels, timeOfInterest = NULL, type = "BucketKM", numBuckets = 10){
  if(length(listOfSurvivalModels) ==0) return(NULL)
  suppressWarnings(if(any(unlist(lapply(listOfSurvivalModels, is.na)))) return(NA))
  if(!type %in% c("Uncensored","Fractional","BucketKM"))
    stop("Please enter one of 'Uncensored','Fractional','BucketKM' for type.")
  
  predictedTimes = lapply(seq_along(listOfSurvivalModels), function(x) listOfSurvivalModels[[x]][[1]]$time)
  survivalCurves = lapply(seq_along(listOfSurvivalModels), function(x) listOfSurvivalModels[[x]][[1]][-1])
  trueDeathTimes = lapply(seq_along(listOfSurvivalModels), function(x) listOfSurvivalModels[[x]][[2]]$time)
  censorStatus   = lapply(seq_along(listOfSurvivalModels), function(x) listOfSurvivalModels[[x]][[2]]$delta)
  #If time is null, we will select the median from the training instances (of those who were uncensored.)
  if(is.null(timeOfInterest)){
    allTimes = unlist(trueDeathTimes)
    allCensorStatus = unlist(censorStatus)
    KMCurve = prodlim(Surv(allTimes, allCensorStatus)~1)
    timeOfInterest = tryCatch({
      quantile(KMCurve,.5)$quantiles.survival$quantile
    },
    #Catch the case where there is too many censored patients to generate a KM Curve. In this case we use our linear extension to find the 
    #Median.
    error = function(e){
      slope = (1-min(KMCurve$surv))/(0 - max(KMCurve$time))
      timeOfInterest = (-0.5)/slope
    })
  }
  #Additionally sometimes the quantile function will return NA instead of an error if it still able to produce a lower bound. Here we 
  #do the same thing as above when we catch the error.
  if(is.na(timeOfInterest)){
    slope = (1-min(KMCurve$surv))/(0 - max(KMCurve$time))
    timeOfInterest = (-0.5)/slope
  }
  predictions = unlist(lapply(seq_along(listOfSurvivalModels), function(model) lapply(seq_along(trueDeathTimes[[model]]),
                              function(index) predictProbabilityFromCurve(survivalCurves[[model]][,index],
                                                                          predictedTimes[[model]],
                                                                          timeOfInterest))))
  pval = binItUp(unlist(trueDeathTimes),unlist(censorStatus), predictions, type, numBuckets,timeOfInterest)
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
                  #taking only uncensored individuals at the time of interest and then taking there mean. All other piping operations should
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
                Fractional = {
                  KM = prodlim(Surv(time,delta)~1, data = survMod[[3]])
                  bucketSurvival = orderedPredictionFrame  %>% 
                    group_by(label) %>%
                    summarise(survivalPrediction = mean(prob))
                  #The Fractional method has three different cases to consider for the observed quantity of an individual.
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
                #For a good discusion of this method, also known as the D'Agostino-Nam tranlation/method see
                #Section 2.2 of this thesis paper: Hosmer-Lemeshow goodness-of-fit test: Translations to the Cox Proportional Hazards Model
                BucketKM = {
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








