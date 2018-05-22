#Date: May 22nd, 2018
#Author: Humza Haider

#Purpose and General Comments:
#This file was originally created to impliment an integrated Brier score for the use in the PSSP paper.
#The reason we need to impliment our own Brier score is that other implimentations have assumed inputs of a single 
#survival curve whereas with PSSP and other ISD models we have specific survivor curves for each patient. 
#Given N patients to evaluate, the input to this evaluation metric is N survival curves, the event time for each patient
#and the censoring status of the patient (1 = event occured, 0 = censored).

#####################################################################################################################################
#Library dependencies:
library(survival)
#prodlim gives a faster KM implimentation and also gives a predict function for KM
library(prodlim)

#Code:
#First we do a single time point brier score and follow with the integrated Brier Score.

singleBrier = function(listOfSurvivalCurves, eventTimes, censorStatus, bTime){
  inverseCensor = 1-censorStatus
  invProbCensor = prodlim(Surv(eventTimes,inverseCensor)~1)
  orderOfTimes = order(eventTimes)
  weightCat1 = (eventTimes[orderOfTimes] <= bTime & censorStatus[orderOfTimes])*predict(invProbCensor, eventTimes, level.chaos = 1)
  weightCat2 = (eventTimes[orderOfTimes] > bTime)*predict(invProbCensor, bTime, level.chaos = 1)
  weightCat2[is.na(weightCat2)] = 0
  
  #TODO: Figure out structure of a "survival curve". Is it a DF with two columns, times and probabilities or...?
  #For now we assume it is a DF with columns time and prob. Further, we will assume a nice curve as opposed to the KM curve.
  predictions = unlist(lapply(listOfSurvivalCurves, function(z) predictFromCurve(z,bTime)))
  bScore = mean(predictions^2*weightCat1 + (1-predictions)^2*weightCat2)
  return(bScore)
}


#integrated Brier Score based on the single time point brier score, found above.
integratedBrier = function(listOfSurvivalCurves, eventTimes, censorStatus, bTime, numPoints){
  inverseCensor = 1-censorStatus
  invProbCensor = prodlim(Surv(eventTimes,inverseCensor)~1)
  orderOfTimes = order(eventTimes)
  points = seq(bTime[1], bTime[2], length.out = numPoints)
  bsPointsMat = matrix(rep(points, length(eventTimes)), nrow = length(eventTimes),byrow = T)
  
  #Each column represents the indicator for a single brier score 
  weightCat1Mat = (eventTimes[orderOfTimes] <= bsPointsMat & censorStatus[orderOfTimes])*predict(invProbCensor, eventTimes, level.chaos = 1)
  weightCat2Mat = (eventTimes[orderOfTimes] > bsPointsMat)*predict(invProbCensor, points,level.chaos = 1)
  #Catch if points goes over max event time, i.e. predict gives NA
  weightCat2Mat[is.na(weightCat2Mat)] = 0
  
  #TODO: Figure out structure of a "survival curve". Is it a DF with two columns, times and probabilities or...?
  #For now we assume it is a DF with columns time and prob. Further, we will assume a nice curve as opposed to the KM curve.
  predictions = matrix(unlist(lapply(listOfSurvivalCurves, function(z) lapply(points, function(point) predictFromCurve(z,point)))),
                       ncol = numPoints, byrow = T)
  bscores = apply(predictions^2*weightCat1Mat + (1-predictions)^2*weightCat2Mat, 2, mean)
  indexs = 2:length(points)
  trapezoidIntegralVal = diff(points) %*% ((bscores[indexs - 1] + bscores[indexs])/2)
  intBScore = trapezoidIntegralVal/(bTime[2] - bTime[1])
  return(intBScore)
}

integratedBrierV2 = function(listOfSurvivalCurves, eventTimes, censorStatus, bTime){
  inverseCensor = 1-censorStatus
  invProbCensor = prodlim(Surv(eventTimes,inverseCensor)~1)
  orderOfTimes = order(eventTimes)
  sortedEvents = sort(eventTimes)
  points = sortedEvents[sortedEvents >= bTime[1] & sortedEvents <=bTime[2]]
  bsPointsMat = matrix(rep(points, length(eventTimes)), nrow = length(eventTimes), byrow = T)
  
  #Each column represents the indicator for a single brier score 
  weightCat1Mat = (eventTimes[orderOfTimes] <= bsPointsMat & censorStatus[orderOfTimes])*predict(invProbCensor, eventTimes, level.chaos = 1)
  weightCat2Mat = (eventTimes[orderOfTimes] > bsPointsMat)*predict(invProbCensor, points,level.chaos = 1)
  #Catch if points goes over max event time, i.e. predict gives NA
  weightCat2Mat[is.na(weightCat2Mat)] = 0
  
  #TODO: Figure out structure of a "survival curve". Is it a DF with two columns, times and probabilities or...?
  #For now we assume it is a DF with columns time and prob. Further, we will assume a nice curve as opposed to the KM curve.
  predictions = matrix(unlist(lapply(listOfSurvivalCurves, function(z) lapply(points, function(point) predictFromCurve(z,point)))),
                       ncol = length(points), byrow = T)
  bscores = apply(predictions^2*weightCat1Mat + (1-predictions)^2*weightCat2Mat, 2, mean)
  indexs = 2:length(points)
  trapezoidIntegralVal = diff(points) %*% ((bscores[indexs - 1] + bscores[indexs])/2)
  intBScore = trapezoidIntegralVal/diff(range(points))
  return(intBScore)
}

#We need some type of predict function for survival curves - here is one assuming the relationship between two time
#points is linear in survival probability.
predictFromCurve = function(survivalCurve, timeToPredict){
  times = sort(survivalCurve$time)
  probs = sort(survivalCurve$probs, decreasing = T)
  timesLess = which(times <= timeToPredict)
  timesGreater = which(times >= timeToPredict)
  
  #These ifelse statements catch if someone inputs a time greater than or less than the beginnning and ending times of the curve.
  indexA = ifelse(length(timesLess) > 0, max(timesLess), 1)
  indexB = ifelse(length(timesGreater) > 0, min(timesGreater),length(times))

  timeA = times[indexA]
  timeB = times[indexB]
  probA = probs[indexA]
  probB = probs[indexB]
  run = timeB - timeA
  rise = probB - probA
  
  if(indexA == indexB){
    #Catches if the timeToPredict is exactly one of the predicted points
    return(probA)
  }

  else{
    toReturn = probA + ((timeToPredict - timeA)*(probB - probA))/(timeB - timeA)
    return(toReturn)
  }
}

















