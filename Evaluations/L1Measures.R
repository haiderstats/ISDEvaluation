#### File Information #####################################################################################################################
#File Name: L1Measures.R
#Date Created: May 29th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#This file was created to implement multiple L measures (L1, L1 hinge, L1 margin) as an evaluation metric for individual survival curves.

### Functions #############################################################################################################################

## Function 1: L1(survMod, type = "Margin", logScale = F, method = "Median")

#Inputs:
#   survMod: A list of 4 items:(1) TestCurves - The survival curves for the testing set.
#                              (2) TestData - The censor/death indicator and event time for the testing set. 
#                              (3) TrainData - The censor/death indicator and event time for the training set. 
#                              (4) TrainCurves - The survival curves for the training set.
#   type: A string specifying the type of L1, options being "Uncensored","Hinge", or "Margin". 
#   logScale: A boolean saying whether or not to use log scale, default is FALSE.
#   method: A string indicating whether the "Mean" or "Median" should be used to calculate a patient's survival time.

# Output: The desired L1-loss.

# Usage: Calculate the L1-loss given a survival model.

### Code #############################################################################################################################
#Library Dependencies:
#We need the Surv function from survival for the prodlim implementation of Kaplan Meier
library(survival)
#We use this for the prodlim function.
library(prodlim)
#Helper Functions: predictMeanSurvivalTimeSpline(survivalCurve,predictedTimes)
source("Evaluations/EvaluationHelperFunctions.R")

L1 = function(survMod, type = "Margin", logScale = F, method = "Median"){
  #Being passed an empty model.
  if(is.null(survMod)) return(NULL)
  #Being passed a model that failed.
  suppressWarnings(if(is.na(survMod[[1]])) return(NULL))
  predictedTimes = survMod[[1]]$time
  survivalCurves = survMod[[1]][-1]
  trueDeathTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  censorTimes = trueDeathTimes[as.logical(1-censorStatus)]
  trainingDeathTimes = survMod[[3]]$time
  trainingCensorStatus = survMod[[3]]$delta
  
  predictMethod = switch(method,
                         Mean = predictMeanSurvivalTimeSpline,
                         Median = predictMedianSurvivalTimeSpline)
  averageUncensored = unlist(lapply(which(as.logical(censorStatus)),
                                    function(index) predictMethod(survivalCurves[,index],
                                                                                  predictedTimes)))
  
  averageCensored = unlist(lapply(which(as.logical(1-censorStatus)),
                                  function(index) predictMethod(survivalCurves[,index],
                                                                                predictedTimes)))
  #Sometimes the mean is infinite or extrememly, unreasonable large. Since this  will largely sway the L1 measure, we enforce that 
  #a model cannot predict a value higher than the end of the KM curve with the linear extension since this is what we believe to be 
  #the "true" survival curve of the population. 
  KMCurve = prodlim(Surv(trainingDeathTimes, trainingCensorStatus)~1)
  KMLinearZero = -1/((1-min(KMCurve$surv))/(0 - max(KMCurve$time)))
  #If every patient is censored we choose the last time point to be the maximum time.
  if(is.infinite(KMLinearZero))
    KMLinearZero = max(KMCurve$time)
  averageUncensored = pmin(averageUncensored, KMLinearZero)
  averageCensored = pmin(averageCensored, KMLinearZero)
  
  
  uncensoredPiece = ifelse(!logScale,
                           sum(abs(trueDeathTimes[as.logical(censorStatus)] - averageUncensored)),
                           sum(abs(log(trueDeathTimes[as.logical(censorStatus)]) - log(averageUncensored))))
  L1Measure = switch(type,
                     Uncensored = {
                       L1Measure = (1/sum(censorStatus))*uncensoredPiece
                     },
                     Hinge = {

                       hingePiece = ifelse(!logScale,
                                           censorTimes - averageCensored,
                                           log(censorTimes) - log(averageCensored))
                       hingePieceCorrected = ifelse(hingePiece >=0, hingePiece,0)
                       L1Measure = (1/(length(censorStatus)))*(uncensoredPiece + hingePieceCorrected)                     
                    },
                     Margin = {
                       KMLinearPredict = function(time){
                         prediction = predict(KMCurve,time)
                         slope = (1-min(KMCurve$surv))/(0 - max(KMCurve$time))
                         predictedProbabiliteis = ifelse(is.na(prediction), pmax(1+time*slope,0), prediction)
                         return(predictedProbabiliteis)
                       }
                       bestGuess = unlist(lapply(censorTimes,
                              function(time) time + integrate(KMLinearPredict,
                                                              lower = time, 
                                                              upper = KMLinearZero,subdivisions = 2000,
                                                              rel.tol = .01)[[1]]/KMLinearPredict(time)))
                       bestGuess[censorTimes > KMLinearZero] = censorTimes[censorTimes > KMLinearZero]
                       weights = 1- KMLinearPredict(censorTimes)
                       marginPiece = ifelse(!logScale,
                                            #Use weights!=0 incase there is an infinitiy, we don't get NaN.
                                            sum(weights*(abs(bestGuess - averageCensored))),
                                            sum(weights*(abs(log(bestGuess) - log(averageCensored))))
                       )
                       L1Measure = (uncensoredPiece+marginPiece)/(sum(censorStatus)+sum(weights))
                     }
  )
  return(L1Measure)
}
