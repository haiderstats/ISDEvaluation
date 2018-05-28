#File Name: EvaluationHelperFunctions.R
#Date Created: May 28th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file should be sourced into the different evaluation methods, it contains functions which may be useful for multiple evaluation methods.
##############################################################################################################################################
#Library Dependencies
#We use prodlim for the sindex function.
library(prodlim)


#We need some type of predict function for survival curves - here is one assuming the relationship between two time
#points is linear in survival probability.
predictProbabilityFromCurve = function(survivalCurve,predictedTimes, timeToPredict){
  lowerTimeIndex = sindex(predictedTimes, timeToPredict)
  upperTimeIndex = ifelse(lowerTimeIndex < length(predictedTimes),lowerTimeIndex+1,lowerTimeIndex)
  
  timeA = predictedTimes[lowerTimeIndex]
  timeB = predictedTimes[upperTimeIndex]
  
  probA = survivalCurve[lowerTimeIndex]
  probB = survivalCurve[upperTimeIndex]
  #Point on a line formula since we assume survival probability is linear between time points.
  #Also, we have an ifelse to catch the event where timeA == timeB, i.e. we are on the last time point, in this case we should add 0.
  toReturn = probA + ifelse(timeB != timeA, ((timeToPredict - timeA)*(probB - probA))/(timeB - timeA),0)
  return(toReturn)
}

#We calculate the mean and median survival times here one for a KM way (stepwise functions) and the other as a more linear approach between
#predicted time points.
predictMeanSurvivalTimeKM = function(survivalCurve, predictedTimes){
  differences = diff(predictedTimes)
  area = sum(differences*survivalCurve)
  return(area)
}

predictMedianSurvivalTimeKM = function(survivalCurve, predictedTimes){
  medianIndex = sindex(as.vector(t(survivalCurve)), 0.5, comp = "greater")+1
  medianTime = predictedTimes[medianIndex]
  return(ifelse(is.na(medianTime), max(predictedTimes), medianTime))
}

predictMeanSurvivalTimeLinear = function(survivalCurve, predictedTimes){
  differences = diff(predictedTimes)
  idx = 1:nrow(survivalCurve)
  #Here we take the area of a right trapezoid. See http://mathworld.wolfram.com/RightTrapezoid.html. Also note we remove the last value (it's
  #NA because we go past the last row of the survivalCurve.)
  area = sum(0.5*differences*(survivalCurve[idx,] + survivalCurve[idx+1,])[-nrow(survivalCurve)])
  return(area)
}

predictMedianSurvivalTimeLinear = function(survivalCurve, predictedTimes){
  medianIndexLower = sindex(as.vector(t(survivalCurve)), 0.5, comp = "greater")
  medianIndexHigher = medianIndexLower +1
  if(is.na(predictedTimes[medianIndexHigher]))
    return(max(predictedTimes))
  else{
    timeA = predictedTimes[medianIndexLower]
    timeB = predictedTimes[medianIndexHigher]
    
    probA = survivalCurve[medianIndexLower,]
    probB = survivalCurve[medianIndexHigher,]
    #Point on a line formula since we assume survival probability is linear between time points.
    #Also, we have an ifelse to catch the event where timeA == timeB, i.e. we are on the last time point, in this case we should add 0.
    toReturn = timeA + (0.5 - probA)*(timeB - timeA)/(probB - probA)
    return(toReturn)
  }
}





