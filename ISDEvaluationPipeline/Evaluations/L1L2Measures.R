#File Name: DCalibration.R
#Date Created: May 29th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file was created to implement multiple L measures (L1, L2, L1 hinge etc.) as an evaluation metric for individual survival curves.
#Input 1: A list of (1) a matrix of survival curves, and (2) the true death times and event/censoring indicator (delta =1 implies death/event).
#Output: The desired L-measure value.
##############################################################################################################################################
#Library Dependencies
#We use this for the sindex function.
library(prodlim)
#We use this for ldply, a combiner of lists.
library(plyr)
#Helper Functions: predictMeanSurvivalTimeKM(survivalCurve,predictedTimes)
source("Evaluations/EvaluationHelperFunctions.R")

L1 = function(survMod, Lmeasure = "meanLinear", type = "Uncensored", logScale = F){
  predictedTimes = survMod[[1]][,1]
  survivalCurves = survMod[[1]][-1]
  trueDeathTimes = survMod[[2]][,1]
  censorStatus = survMod[[2]][,2]
  
  averageMeasure = switch(Lmeasure,
                          meanLinear = predictMeanSurvivalTimeLinear,
                          meanKM = predictMeanSurvivalTimeKM,
                          medianLinear = predictMedianSurvivalTimeLinear,
                          medianKM = predictMedianSurvivalTimeKM)
  averageUncensored = unlist(lapply(which(as.logical(censorStatus)),
                                    function(index) averageMeasure(survivalCurves[,index],
                                                                                  predictedTimes)))
  L1Measure = switch(type,
                     Uncensored = {
                       L1Measure = ifelse(!logScale,
                                          (1/(sum(censorStatus)))*sum(abs(trueDeathTimes[as.logical(censorStatus)] - averageUncensored)),
                                          (1/(sum(censorStatus)))*sum(abs(log(trueDeathTimes[as.logical(censorStatus)]) - 
                                                                            log(averageUncensored))))
                     },
                     Hinge = {
                       averageCensored = unlist(lapply(which(as.logical(1-censorStatus)),
                                                       function(index) averageMeasure(survivalCurves[,index],
                                                                                      predictedTimes)))
                       hingePiece = trueDeathTimes[as.logical(1-censorStatus)] - averageCensored
                       hingePieceCorrected = ifelse(hingePiece >=0, hingePiece,0)
                       L1Measure = ifelse(!logScale,
                                          (1/(length(censorStatus)))*(sum(abs(trueDeathTimes[as.logical(censorStatus)] - averageUncensored)) +
                                                                        sum(hingePieceCorrected)),
                                          (1/(length(censorStatus)))*(sum(abs(log(trueDeathTimes[as.logical(censorStatus)]) - 
                                                                                log(averageUncensored))) +
                                                                        sum(log(hingePieceCorrected[hingePieceCorrected!=0]))))
                     }
  )
}


L2 = function(survMod, Lmeasure = "meanLinear", type = "Uncensored", logScale = F){
  predictedTimes = survMod[[1]][,1]
  survivalCurves = survMod[[1]][-1]
  trueDeathTimes = survMod[[2]][,1]
  censorStatus = survMod[[2]][,2]
  
  averageMeasure = switch(Lmeasure,
                          meanLinear = predictMeanSurvivalTimeLinear,
                          meanKM = predictMeanSurvivalTimeKM,
                          medianLinear = predictMedianSurvivalTimeLinear,
                          medianKM = predictMedianSurvivalTimeKM)
  averageUncensored = unlist(lapply(which(as.logical(censorStatus)),
                                    function(index) averageMeasure(survivalCurves[,index],
                                                                   predictedTimes)))
  L2Measure = switch(type,
                     Uncensored = {
                       L2Measure = ifelse(!logScale,
                                          (1/(sum(censorStatus)))*sum((trueDeathTimes[as.logical(censorStatus)] - averageUncensored)^2),
                                          (1/(sum(censorStatus)))*sum((log(trueDeathTimes[as.logical(censorStatus)]) -
                                                                         log(averageUncensored))^2))
                     },
                     Hinge = {
                       
                       averageCensored = unlist(lapply(which(as.logical(1-censorStatus)),
                                                       function(index) averageMeasure(survivalCurves[,index],
                                                                                      predictedTimes)))
                       hingePiece = trueDeathTimes[as.logical(1-censorStatus)] - averageCensored
                       hingePieceCorrected = ifelse(hingePiece >=0, hingePiece,0)
                       L2Measure = ifelse(!logScale,
                                          (1/(length(censorStatus)))*(sum((trueDeathTimes[as.logical(censorStatus)] - averageUncensored)^2) +
                                                                        sum(hingePieceCorrected)),
                                          (1/(length(censorStatus)))*(sum((log(trueDeathTimes[as.logical(censorStatus)]) -
                                                                             log(averageUncensored))^2) +
                                                                        sum(log(hingePieceCorrected[hingePieceCorrected!=0])))
                       )
                     }
  )
    return(L2Measure)
}








