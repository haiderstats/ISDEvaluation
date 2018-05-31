#File Name: DCalibration.R
#Date Created: May 28th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file was created to implement D-Calibration as an evaluation measure for individual survival curves.
#Input 1: A list of (1) a matrix of survival curves, and (2) the true death times and event/censoring indicator (delta =1 implies death/event).
#Input 2: The number of bins to evaluate D-calibration.
#Output: The p-value associated with D-calibration.
##############################################################################################################################################
#Library Dependencies
#We use this for the sindex function.
library(prodlim)
#We use this for ldply, a combiner of lists.
library(plyr)
#Helper Functions: predictProbabilityFromCurve(survivalCurve,predictedTimes, timeToPredict)
source("Evaluations/EvaluationHelperFunctions.R")

DCalibration = function(survMod, numBins){
  if(is.null(survMod)) return(NULL)
  combinedBins = getBinned(survMod, numBins)
  pvalue = chisq.test(combinedBins)$p.value
  return(pvalue)
}

DCalibrationCumulative = function(listOfSurvivalModels, numBins){
  if(length(listOfSurvivalModels) == 0) return(NULL)
  #Here we apply getBinned to every survival model, stack the bins into a matrix, and sum the columns to get the total bin values.
  combinedBins =colSums(ldply(lapply(seq_along(listOfSurvivalModels), function(x) getBinned(listOfSurvivalModels[[x]], numBins)), rbind))
  pvalue = chisq.test(combinedBins)$p.value
  return(pvalue)
}


getBinned = function(survMod,numBins){
  predictedTimes = survMod[[1]][,1]
  survivalCurves = survMod[[1]][-1]
  trueDeathTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  quantiles = 1 - seq(0,1,length.out = numBins+1)
  #This retrieves the death probability the survival curve gave at the true time of death.
  deathProbabilities = unlist(lapply(seq_along(trueDeathTimes),
                                     function(index) predictProbabilityFromCurve(survivalCurves[,index],
                                                                                 predictedTimes,
                                                                                 trueDeathTimes[index])))
  uncensoredProbabilities = deathProbabilities[as.logical(censorStatus)]
  binIndex = seq_along(quantiles)[-1]
  uncensoredBinning = unlist(lapply(binIndex, function(x) length(which(uncensoredProbabilities >= quantiles[x] &
                                                                         uncensoredProbabilities < quantiles[x-1]))))
  censoredProbabilities = deathProbabilities[as.logical(1-censorStatus)]
  censoredBinPositions = sindex(quantiles, censoredProbabilities)
  percentToAdd = 1/((numBins - censoredBinPositions +1))
  #This is a list where each element is is a numBins length vector with 0s prior to the bin where the censored individual landed and 
  #a fraction that they contribute to each bin in all the following bins. This list is then composed into a matrix and columns are summed
  #to add the approprate contribution to each bin for the censored individuals.
  listOfContributions = lapply(seq_along(censoredBinPositions),function(x) c(rep(0, censoredBinPositions[x] -1),
                                                                             rep(percentToAdd[x], numBins - censoredBinPositions[x] +1)))
  censoredBinning = colSums(ldply(listOfContributions,rbind))
  combinedBins = uncensoredBinning + censoredBinning
  return(combinedBins)
}







