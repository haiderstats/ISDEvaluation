#### File Information #####################################################################################################################
#File Name: BrierScore.R
#Date Created: May 28th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#This file was created to implement a single time and integrated time Brier score. The integrated Brier Score is evaluted numerically 
#two different ways. Either one can specifcy the number of points to evaluate the brier score and a uniform selection of these points will
#be chosen from 0 to the max event time. If numPoints is left to be null then the event times in the dataset will be evaluated.

### Functions #############################################################################################################################


## Function 1: BrierScore(survMod, type = "Integrated", singleBrierTime = NULL, integratedBrierTimes = NULL, numPoints=NULL)

#Inputs:
#   survMod: A list of 4 items:(1) TestCurves - The survival curves for the testing set.
#                            (2) TestData - The censor/death indicator and event time for the testing set. 
#                            (3) TrainData - The censor/death indicator and event time for the training set. 
#                            (4) TrainCurves - The survival curves for the training set.
#   type: The type of Brier score to compute. Either "Single" or "Integrated"
#   singleBrierTime: If "Single" is selected, indicate the time at which to evaluate the Brier Score. If left as NULL,
#                    the median event time is used.
#   integratedBrierTimes: A vector of length 2 indicating the interval to use for the Brier score (e.g. c(0,100)) would be evaluating
#                         from t= 0 to t = 100. If left NULL then the interval is t = 0 to the maximum event time (using both the
#                         training and testing sets).
#   numPoints: The number of points to use for the integrated Brier score. If left NULL then the points will correspond to the event times 
#               in the testing dataset.

# Output: The desired Brier score.

# Usage: Evaluate the Brier score on a given survival model, e.g. MTLR, CoxENKP, RSF.


## Function 2: singleBrier(survMod,singleBrierTime)

# Inputs: See BrierScore().

# Output: The Brier score at the indicated time.

# Usage: Evaluate the SINGLE TIME Brier score on a given survival model at a given time. (Helper Function for BrierScore).


## Function 3: integratedBrier(survMod,integratedBrierTimes,numPoints)

# Inputs: See BrierScore().

# Output: The Integrated Brier Score for the given interval.

# Usage: Evaluate the Integrated Brier score on a given survival model at a given time. (Helper Function for BrierScore).


## Function 4: singleBrierMultiplePoints(survMod,BrierTimes)

# Inputs:
#   survMod: See BrierScore().
#   BrierTimes: A vector of time to evaluate the Brier Score.

# Output: The single time Brier score evaluated at all BrierTimes.

# Usage: Evaluate the Single Brier score for an entire vector of times (does not work for a single time).
#        (Helper Function for integratedBrier).


### Code #############################################################################################################################
#Library Dependencies:
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

