#File Name: BrierScore.R
#Date Created: May 28th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file was created to implement a single time and integrated time Brier score. The integrated Brier Score is evaluted numerically 
#two different ways. One option is to use the time points at the event times in the testing data. This option would correspond to 
#basedOnEvents =T. In this case, the function integratedBrier() is used. In this event, numerical integration is doen using the trapezoid rule(
#see: https://en.wikipedia.org/wiki/Trapezoidal_rule). This is what is does when calling the sbrier function. 
#The other option is to simply use the built in integrate() function which is an implementation of a numerical integration technique. 
#This evaluates enough points in the interval to produce a numerical integral with low error. 
#This function uses the singleBrierMultiplePoints() function. If only one point is evaluated the singleBrier() function is used.
#Ideally, singleBrier and singleBrierMultiplePoints would be built into one function but this has not been done yet. 

#survMod: A list of (1) a matrix of survival curves, and (2) the true death times and event/censoring indicator (delta =1 implies death/event)
#of the TESTING data, and (3) the true death times and event/censoring indicator (delta =1 implies death/event) of the TRAINING data.

#type: A string indicating whether the Brier score will be "Single" or "Integrated"

#singleBrierTime: The time for the Single Brier score to be evaluated, a single numeric value, e.g. 42.

#integratedBrierTimes: A vector of length two giving the limits of integration for the integrated Brier Score.

#numPoints: The number of points to numerically evaluate the Brier score -- the larger the number of points, the more accurate the integral.

#Output: The desired type of Brier Score (single or integrated).
###############################################################################################################################################
#Library dependencies:
#prodlim gives a faster KM implementation and also gives a predict function for KM
library(prodlim)
#Helper functions:
source("Evaluations/EvaluationHelperFunctions.R")

BrierScore = function(survMod, type = "Integrated", singleBrierTime = NULL, integratedBrierTimes = NULL, numPoints=NULL){
  #Being passed an empty model.
  if(is.null(survMod)) return(NULL)
  #Being passed a model that failed.
  suppressWarnings(if(is.na(survMod[[1]])) return(NULL))
  if(type == "Integrated" & is.null(integratedBrierTimes))
    integratedBrierTimes = c(0, max(c(survMod[[2]]$time, survMod[[3]]$time))) #max time of the training and testing set combined (the entire dataset).
  score = ifelse(type =="Single",
                 singleBrier(survMod, singleBrierTime),
                 integratedBrier(survMod, integratedBrierTimes, numPoints))
  return(score)
}

singleBrier = function(survMod, singleBrierTime){
  eventTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  trainingEventTimes = survMod[[3]]$time
  trainingCensorStatus = survMod[[3]]$delta
  #Default brier time will be the 50th quantile of training and testing event times combined.
  if(is.null(singleBrierTime)){
    singleBrierTime = quantile(c(eventTimes, trainingEventTimes), .5)
  }
  inverseCensorTrain = 1 - trainingCensorStatus
  invProbCensor = prodlim(Surv(trainingEventTimes,inverseCensorTrain)~1)
  #Here we are ordering event times and then using predict with level.chaos = 1 which returns predictions ordered by time.
  orderOfTimes = order(eventTimes)
  #Category one is individuals with event time lower than the time of interest and were NOT censored.
  weightCat1 = (eventTimes[orderOfTimes] <= singleBrierTime & censorStatus[orderOfTimes])*predict(invProbCensor,
                                                                                            eventTimes,
                                                                                            level.chaos = 1)
  #Catch if event times goes over max training event time, i.e. predict gives NA
  weightCat1[is.na(weightCat1)] = 0
  #Category 2 is individuals whose time was greater than the time of interest (singleBrierTime) - both censored and uncensored individuals.
  weightCat2 = (eventTimes[orderOfTimes] > singleBrierTime)*predict(invProbCensor,
                                                                    singleBrierTime,
                                                              level.chaos = 1)
  #predict returns NA if the passed in time is greater than any of the times used to build the inverse probability of censoring model.
  weightCat2[is.na(weightCat2)] = 0
  
  predictedTimes = survMod[[1]][,1]
  #Take the survival curves, remove the times column, and then order the curves by the order in which they died. We have to order them to 
  #line up with the weight vectors.
  survivalCurvesOrdered = survMod[[1]][,-1][orderOfTimes]
  predictions = apply(survivalCurvesOrdered,2, function(z) predictProbabilityFromCurve(z,predictedTimes,singleBrierTime))
  bScore = mean(predictions^2*weightCat1 + (1-predictions)^2*weightCat2)
  return(bScore)
}

#integrated Brier Score based on the single time point brier score, found above.
integratedBrier = function(survMod, integratedBrierTimes,numPoints){
  censorStatus = survMod[[2]]$delta
  eventTimes = survMod[[2]]$time[censorStatus==1]
  orderOfTimes = order(eventTimes)
  sortedEvents = sort(eventTimes)
  if(is.null(numPoints)){
    points = sortedEvents[sortedEvents >= integratedBrierTimes[1] & sortedEvents <=integratedBrierTimes[2]]
  } else{
    points = seq(integratedBrierTimes[1], integratedBrierTimes[2], length.out = numPoints)
  }
  bscores = singleBrierMultiplePoints(survMod,points)
  indexs = 2:length(points)
  trapezoidIntegralVal = diff(points) %*% ((bscores[indexs - 1] + bscores[indexs])/2)
  intBScore = trapezoidIntegralVal/diff(range(points))
  return(intBScore)
}

singleBrierMultiplePoints = function(survMod, BrierTimes){
  eventTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  trainingEventTimes = survMod[[3]]$time
  trainingCensorStatus = survMod[[3]]$delta
  inverseCensorTrain = 1 - trainingCensorStatus
  invProbCensor = prodlim(Surv(trainingEventTimes,inverseCensorTrain)~1)
  orderOfTimes = order(eventTimes)

  bsPointsMat = matrix(rep(BrierTimes, length(eventTimes)), nrow = length(eventTimes),byrow = T)
  
  #Each column represents the indicator for a single brier score 
  weightCat1Mat = (eventTimes[orderOfTimes] <= bsPointsMat & censorStatus[orderOfTimes])*predict(invProbCensor, eventTimes, level.chaos = 1)
  #Catch if event times goes over max training event time, i.e. predict gives NA
  weightCat1Mat[is.na(weightCat1Mat)] = 0

  weightCat2Mat = t(t((eventTimes[orderOfTimes] > bsPointsMat))*predict(invProbCensor, BrierTimes,level.chaos = 2))
  #Catch if BrierTimes goes over max event time, i.e. predict gives NA
  weightCat2Mat[is.na(weightCat2Mat)] = 0
  
  predictedTimes = survMod[[1]][,1]
  #Take the survival curves, remove the times column, and then order the curves by the order in which they died. We have to order them to
  #line up with the weight matricies.
  survivalCurvesOrdered = survMod[[1]][,-1][orderOfTimes]
  predictions = t(apply(survivalCurvesOrdered,2, function(curve) predictProbabilityFromCurve(curve,predictedTimes,BrierTimes)))
  bscores = apply(predictions^2*weightCat1Mat + (1-predictions)^2*weightCat2Mat, 2, mean)
  return(bscores)
}

