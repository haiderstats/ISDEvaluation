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

#We need some type of predict function for survival curves - here we build a spline to fit the survival model curve. This spline is 
#the montotone spline using the Fritsch-Carlson method, see https://en.wikipedia.org/wiki/Monotone_cubic_interpolation. Also see 
#help(splinefun), specifically for method = "monoH.FC". Note that we make an alteration to the method because if the last two time points
#have the same probability (y value) then the spline is constant outside of the training data. We need this to be a decreasing function
#outside the training data so instead we take the linear fit of (0,1) and the last time point we have (p,t*) and then apply this linear
#function to all points outside of our fit.
predictProbabilityFromCurve = function(survivalCurve,predictedTimes, timeToPredict){
  spline = splinefun(predictedTimes, survivalCurve, method = "hyman")
  maxTime = max(predictedTimes)
  slope = (1-spline(maxTime))/(0 - max(predictedTimes))
  predictedProbabilities = rep(0, length(timeToPredict))
  linearChange = which(timeToPredict > maxTime)
  if(length(linearChange) > 0){
    predictedProbabilities[linearChange] = pmax(1 + timeToPredict[linearChange]*slope,0)
    predictedProbabilities[-linearChange] = spline(timeToPredict[-linearChange])
  }
  else{
    predictedProbabilities = spline(timeToPredict)
  }
  return(predictedProbabilities)
}

#We calculate the mean and median survival times assuming a monotone spline fit of the survival curve points.
predictMeanSurvivalTimeSpline = function(survivalCurve, predictedTimes){
  #If all the predicted probabilities are 1 the integral will be infinite. For this reason we slightly decrease the 
  #last value.
  if(all(survivalCurve==1)){
    survivalCurve[length(survivalCurve)] = 1 - 1e-5
  }
  spline = splinefun(predictedTimes, survivalCurve, method = "hyman")
  maxTime = max(predictedTimes)
  slope = (1-spline(maxTime))/(0 - max(predictedTimes))
  zeroProbabilitiyTime = maxTime + (0-spline(maxTime))/slope
  splineWithLinear = function(time) ifelse(time < maxTime, spline(time),1 + time*slope)
  area = integrate(splineWithLinear,0, zeroProbabilitiyTime,subdivisions = 1000,rel.tol = .001)[[1]]
  return(area)
}

predictMedianSurvivalTimeSpline = function(survivalCurve, predictedTimes){
  #If all the predicted probabilities are 1 the integral will be infinite. For this reason we slightly decrease the 
  #last value.
  if(all(survivalCurve==1)){
    survivalCurve[length(survivalCurve)] = 1 - 1e-5
  }
  spline = splinefun(predictedTimes, survivalCurve, method = "hyman")
  minProb = min(spline(predictedTimes))
  if(minProb < .5){
    maximumSmallerThanMedian = predictedTimes[min(which(survivalCurve <.5))]
    minimumGreaterThanMedian = predictedTimes[max(which(survivalCurve >.5))]
    splineInv = splinefun(spline(seq(minimumGreaterThanMedian, maximumSmallerThanMedian, length.out = 1000)),
                          seq(minimumGreaterThanMedian, maximumSmallerThanMedian, length.out = 1000))
    medianProbabilityTime = splineInv(0.5)
  }
  else{
    slope = (1-spline(maxTime))/(0 - max(predictedTimes))
    medianProbabilityTime = maxTime + (0.5-spline(maxTime))/slope
  }
  return(medianProbabilityTime)
}



#Currently these functions are not being used but I have kept them here for the time being incase we ever want to use them as opposed
#to the spline model.

#We calculate the mean and median survival times assuming a stepwise function (the KM type curve)
predictMeanSurvivalTimeKM = function(survivalCurve, predictedTimes){
  differences = diff(predictedTimes)
  area = sum(differences*survivalCurve[-length(survivalCurve)])
  return(area)
}

predictMedianSurvivalTimeKM = function(survivalCurve, predictedTimes){
  medianIndex = sindex(survivalCurve, 0.5, comp = "greater")+1
  medianTime = predictedTimes[medianIndex]
  return(ifelse(is.na(medianTime), max(predictedTimes), medianTime))
}


#We calculate the mean and median survival times assuming a linear function between time points.
predictMeanSurvivalTimeLinear = function(survivalCurve, predictedTimes){
  differences = diff(predictedTimes)
  idx = 1:length(survivalCurve)
  #Here we take the area of a right trapezoid. See http://mathworld.wolfram.com/RightTrapezoid.html. Also note we remove the last value (it's
  #NA because we go past the last row of the survivalCurve.)
  area = sum(0.5*differences*(survivalCurve[idx] + survivalCurve[idx+1])[-length(survivalCurve)])
  return(area)
}

predictMedianSurvivalTimeLinear = function(survivalCurve, predictedTimes){
  medianIndexLower = sindex(survivalCurve, 0.5, comp = "greater")
  medianIndexHigher = medianIndexLower +1
  if(is.na(predictedTimes[medianIndexHigher]))
    return(max(predictedTimes))
  else{
    timeA = predictedTimes[medianIndexLower]
    timeB = predictedTimes[medianIndexHigher]
    
    probA = survivalCurve[medianIndexLower]
    probB = survivalCurve[medianIndexHigher]
    #Point on a line formula since we assume survival probability is linear between time points.
    #Also, we have an ifelse to catch the event where timeA == timeB, i.e. we are on the last time point, in this case we should add 0.
    toReturn = timeA + (0.5 - probA)*(timeB - timeA)/(probB - probA)
    return(toReturn)
  }
}





