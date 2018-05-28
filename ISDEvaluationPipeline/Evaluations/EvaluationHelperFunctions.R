#File Name: EvaluationHelperFunctions.R
#Date Created: May 28th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file should be sourced into the different evaluation methods, it contains functions which may be useful for multiple evaluation methods.
##############################################################################################################################################

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

