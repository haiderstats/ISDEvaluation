#File Name: BrierScore.R
#Date Created: May 28th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file was created to implement a single time and integrated time Brier score. The integrated Brier Score is evaluted numerically 
#by taking brier scores at a number of time points and then integrating across different time points by taking the trapezoidal 
#approximation (see: https://en.wikipedia.org/wiki/Trapezoidal_rule).

#Input 1: A list of (1) a matrix of survival curves, and (2) the true death times and event/censoring indicator (delta =1 implies death/event).
#Input 2: The time for the Brier score to be evaluated. Either a vector of length 2 e.g. c(0,100) or a single numeric value, e.g. 42.
#Input 3: The number of points to be used for the integrated Brier Score. In the event that nothing is passed in for the number of points
#AND the time is a vector of length two, the time points will be only those event times which occured between the first and second 
#time point. 
#For example, if the death times of the test cases were 5, 10, 42, and 100 and the time passed in was c(6, 75) then only 2 time points
#would be used, namely, 10 and 42. 
#Output: The desired type of Brier Score.
############################################################################################################################################
#Library dependencies:
#prodlim gives a faster KM implementation and also gives a predict function for KM
library(prodlim)
#Helper functions:
source("Evaluations/EvaluationHelperFunctions.R")

BrierScore = function(survMod, BrierTime, numberBrierPoints = NULL){
  score = ifelse(length(BrierTime) ==1,
                 singleBrier(survMod, BrierTime),
                 integratedBrier(survMod, BrierTime, numberBrierPoints))
  return(score)
}

singleBrier = function(survMod, BrierTime){
  eventTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  inverseCensor = 1-censorStatus
  invProbCensor = prodlim(Surv(eventTimes,inverseCensor)~1)
  #Here we are ordering event times and then using predict with level.chaos = 1 which returns predictions ordered by time.
  orderOfTimes = order(eventTimes)
  #Category one is individuals with event time lower than the time of interest and were NOT censored.
  weightCat1 = (eventTimes[orderOfTimes] <= BrierTime & censorStatus[orderOfTimes])*predict(invProbCensor,
                                                                                            eventTimes,
                                                                                            level.chaos = 1)
  #Category 2 is individuals whose time was greater than the time of interest (BrierTime) - both censored and uncensored individuals.
  weightCat2 = (eventTimes[orderOfTimes] > BrierTime)*predict(invProbCensor,
                                                              BrierTime,
                                                              level.chaos = 1)
  #predict returns NA if the passed in time is greater than any of the times used to build the inverse probability of censoring model.
  weightCat2[is.na(weightCat2)] = 0
  
  predictedTimes = survMod[[1]][,1]
  #Take the survival curves, remove the times column, and then order the curves by the order in which they died. We have to order them to line up
  #with the weight vectors.
  survivalCurvesOrdered = survMod[[1]][,-1][orderOfTimes]
  predictions = apply(survivalCurvesOrdered,2, function(z) predictProbabilityFromCurve(z,predictedTimes,BrierTime))
  bScore = mean(predictions^2*weightCat1 + (1-predictions)^2*weightCat2)
  return(bScore)
}


#integrated Brier Score based on the single time point brier score, found above.
integratedBrier = function(survMod, BrierTime, numPoints = NULL){
  eventTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  inverseCensor = 1-censorStatus
  invProbCensor = prodlim(Surv(eventTimes,inverseCensor)~1)
  orderOfTimes = order(eventTimes)
  if(is.null(numPoints)){
    sortedEvents = sort(eventTimes)
    points = sortedEvents[sortedEvents >= BrierTime[1] & sortedEvents <=BrierTime[2]]
  }
  else{
    points = seq(BrierTime[1], BrierTime[2], length.out = numPoints)
  } 
  bsPointsMat = matrix(rep(points, length(eventTimes)), nrow = length(eventTimes),byrow = T)
  
  #Each column represents the indicator for a single brier score 
  weightCat1Mat = (eventTimes[orderOfTimes] <= bsPointsMat & censorStatus[orderOfTimes])*predict(invProbCensor, eventTimes, level.chaos = 1)
  weightCat2Mat = (eventTimes[orderOfTimes] > bsPointsMat)*predict(invProbCensor, points,level.chaos = 1)
  #Catch if points goes over max event time, i.e. predict gives NA
  weightCat2Mat[is.na(weightCat2Mat)] = 0
  
  predictedTimes = survMod[[1]][,1]
  #Take the survival curves, remove the times column, and then order the curves by the order in which they died. We have to order them to line up
  #with the weight matricies.
  survivalCurvesOrdered = survMod[[1]][,-1][orderOfTimes]
  predictions = t(apply(survivalCurvesOrdered,2, function(curve) predictProbabilityFromCurve(as.vector(t(curve)),predictedTimes,points)))
  bscores = apply(predictions^2*weightCat1Mat + (1-predictions)^2*weightCat2Mat, 2, mean)
  indexs = 2:length(points)
  trapezoidIntegralVal = diff(points) %*% ((bscores[indexs - 1] + bscores[indexs])/2)
  if(is.null(numPoints)){
    intBScore = trapezoidIntegralVal/diff(range(points))
  }
  else{
    intBScore = trapezoidIntegralVal/(BrierTime[2] - BrierTime[1])
  }
  return(intBScore)
}

















