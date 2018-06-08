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

#Input 1: A list of (1) a matrix of survival curves, and (2) the true death times and event/censoring indicator (delta =1 implies death/event)
#of the TESTING data, and (3) the true death times and event/censoring indicator (delta =1 implies death/event) of the TRAINING data.
#Input 2: The time for the Brier score to be evaluated. Either a vector of length 2 e.g. c(0,100) or a single numeric value, e.g. 42.
#Input 3: A boolean indicating if integration should be done using the event times of the testing data or if we should integrate over all
#possible time points in the vector. In the first case an example is if the death times of the test cases were 5, 10, 42, and 100 and the time
#passed in was c(6, 75) then only 2 time points would be used, namely, 10 and 42. In the second case we numerically integrate over c(6,75).
#Output: The desired type of Brier Score (single or integrated).
###############################################################################################################################################
#Library dependencies:
#prodlim gives a faster KM implementation and also gives a predict function for KM
library(prodlim)
#Helper functions:
source("Evaluations/EvaluationHelperFunctions.R")

BrierScore = function(survMod, BrierTime = NULL, basedOnEvents=F){
  if(is.null(survMod)) return(NULL)
  if(is.null(BrierTime)) 
    BrierTime = c(0, max(c(survMod[[2]]$time, survMod[[3]]$time))) #max time of the training and testing set combined (the entire dataset).
  score = ifelse(length(BrierTime) ==1,
                 singleBrier(survMod, BrierTime),
                 ifelse(basedOnEvents,integratedBrier(survMod, BrierTime),
                        integrate(singleBrierMultiplePoints, lower= BrierTime[1],upper= BrierTime[2],survMod=survMod,
                                  subdivisions = 700)[[1]]/diff(BrierTime)))
  return(score)
}

singleBrier = function(survMod, BrierTime){
  eventTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  trainingEventTimes = survMod[[3]]$time
  trainingCensorStatus = survMod[[3]]$delta
  inverseCensorTrain = 1 - trainingCensorStatus
  invProbCensor = prodlim(Surv(trainingEventTimes,inverseCensorTrain)~1)
  #Here we are ordering event times and then using predict with level.chaos = 1 which returns predictions ordered by time.
  orderOfTimes = order(eventTimes)
  #Category one is individuals with event time lower than the time of interest and were NOT censored.
  weightCat1 = (eventTimes[orderOfTimes] <= BrierTime & censorStatus[orderOfTimes])*predict(invProbCensor,
                                                                                            eventTimes,
                                                                                            level.chaos = 1)
  #Catch if event times goes over max training event time, i.e. predict gives NA
  weightCat1[is.na(weightCat1)] = 0
  #Category 2 is individuals whose time was greater than the time of interest (BrierTime) - both censored and uncensored individuals.
  weightCat2 = (eventTimes[orderOfTimes] > BrierTime)*predict(invProbCensor,
                                                              BrierTime,
                                                              level.chaos = 1)
  #predict returns NA if the passed in time is greater than any of the times used to build the inverse probability of censoring model.
  weightCat2[is.na(weightCat2)] = 0
  
  predictedTimes = survMod[[1]][,1]
  #Take the survival curves, remove the times column, and then order the curves by the order in which they died. We have to order them to 
  #line up with the weight vectors.
  survivalCurvesOrdered = survMod[[1]][,-1][orderOfTimes]
  predictions = apply(survivalCurvesOrdered,2, function(z) predictProbabilityFromCurve(z,predictedTimes,BrierTime))
  bScore = mean(predictions^2*weightCat1 + (1-predictions)^2*weightCat2)
  return(bScore)
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

#integrated Brier Score based on the single time point brier score, found above.
integratedBrier = function(survMod, BrierTime){
  eventTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  trainingEventTimes = survMod[[3]]$time
  trainingCensorStatus = survMod[[3]]$delta
  inverseCensorTrain = 1 - trainingCensorStatus
  invProbCensor = prodlim(Surv(trainingEventTimes,inverseCensorTrain)~1)
  orderOfTimes = order(eventTimes)
  sortedEvents = sort(eventTimes)
  points = sortedEvents[sortedEvents >= BrierTime[1] & sortedEvents <=BrierTime[2]]
  bsPointsMat = matrix(rep(points, length(eventTimes)), nrow = length(eventTimes),byrow = T)
  
  #Each column represents the indicator for a single brier score 
  weightCat1Mat = (eventTimes[orderOfTimes] <= bsPointsMat & censorStatus[orderOfTimes])*predict(invProbCensor, eventTimes, level.chaos = 1)
  #Catch if event times goes over max training event time, i.e. predict gives NA
  weightCat1Mat[is.na(weightCat1Mat)] = 0
  
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
  intBScore = trapezoidIntegralVal/diff(range(points))
  return(intBScore)
}

















