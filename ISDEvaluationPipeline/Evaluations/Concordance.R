#File Name: Concordance.R
#Date Created: May 28th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file was created to implement concordance as an evaluation measure for individual survival curves.
#Input 1: A list of (1) a matrix of survival curves, and (2) the true death times and event/censoring indicator (delta =1 implies death/event).
#Output: The C-index score.
##############################################################################################################################################
#Library Dependencies
#We use this for the sindex function.
library(prodlim)
#We use this for ldply, a combiner of lists.
library(plyr)
#Helper Functions: predictProbabilityFromCurve(survivalCurve,predictedTimes, timeToPredict)
source("Evaluations/EvaluationHelperFunctions.R")

Concordance = function(survMod){
  predictedTimes = survMod[[1]][,1]
  survivalCurves = survMod[[1]][-1]
  trueDeathTimes = survMod[[2]][,1]
  censorStatus = survMod[[2]][,2]
  #This retrieves the death probability of the survival curve at the true time of death.
  averageSurvivalTimes = unlist(lapply(seq_along(trueDeathTimes),
                                     function(index) predictMeanSurvivalTimeLinear(survivalCurves[,index],
                                                                                 predictedTimes)))
  orderDeathTimes = order(trueDeathTimes)
  sortedDeathTimes = sort(trueDeathTimes)
  #TODO:
  #Currently we use the mean time as a "risk score". How do survival curves produce concordance? We have survival probabilities for every time...
  #but what time should we use? When an individual dies.... and compare it to what - the same time for every other individual?
  #Another option is to test if mean/median survival times are concordant. If one individual dies before another we would expect
  #their average expected time should be lower.
  #TODO:
  #How should we handle ties?
  trueMatrix = ldply(lapply(seq_along(sortedDeathTimes), function(x) as.numeric(sortedDeathTimes[x] < sortedDeathTimes  &
                                                                                  censorStatus[orderDeathTimes][x])),rbind)
  estimatedMatrix = ldply(lapply(seq_along(sortedDeathTimes),
                                 function(x) as.numeric(sortedDeathTimes[x] < sortedDeathTimes  &
                                                        censorStatus[orderDeathTimes][x]&
                                                          averageSurvivalTimes[orderDeathTimes][x] < averageSurvivalTimes[orderDeathTimes])),rbind)
  sum(estimatedMatrix)/sum(trueMatrix)
   return()
}








