#### File Information #####################################################################################################################
#File Name: DCalibration.R
#Date Created: May 28th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#This file was created to implement D-Calibration as an evaluation measure for individual survival curves.
#Currently there are two options, DCalibration and DCalibrationCumulative. Since we are doing Cross validation it doesn't make much sense
#to be averaging p-values. For this reason we have DCalibrationCumulative which takes in a list of lists containing survival curves, one entry
#for each fold of the cross validation. 

### Functions #############################################################################################################################

## Function 1: DCalibration(survMod, numBins = 10)

#Inputs:
#   survMod: A list of 4 items:(1) TestCurves - The survival curves for the testing set.
#                              (2) TestData - The censor/death indicator and event time for the testing set. 
#                              (3) TrainData - The censor/death indicator and event time for the training set. 
#                              (4) TrainCurves - The survival curves for the training set.
#   numBins: The number of Buckets/Bins/Groups to use for D-Calibration.

# Output: The p-value for D-Calibration (for a given test set).

# Usage: Calculate D-Calibration on a test set given a survival model.


## Function 2: DCalibrationCumulative(listOfSurvivalModels, numBins = 10)

# Inputs:
#   listOfSurvivalModels: A list of models, each corresponding to survMod in DCalibration()
#   numBins:  The number of Buckets/Bins/Groups to use for D-Calibration.

# Output: The p-value for D-Calibration (across an entire dataset).

# Usage: Calculate D-Calibration on the entire dataset using all folds of data.


## Function 3: getBinned(survMod, numBins)

# Inputs: See DCalibration()

# Output: The counts for each bin.

# Usage: Given a survival model get death counts for each bin corresponding to the bins percentiles. (Helper function for DCalibration/
#        DCalibrationCumulative).

### Code ##################################################################################################################################
#Library Dependencies:
#We use this for the sindex function.
library(prodlim)
#We use this for ldply, a combiner of lists.
library(plyr)
#Helper Functions: predictProbabilityFromCurve(survivalCurve,predictedTimes, timeToPredict)
source("Evaluations/EvaluationHelperFunctions.R")

DCalibration = function(survMod, numBins = 10){
  #Being passed an empty model.
  if(is.null(survMod)) return(NULL)
  #Being passed a model that failed.
  suppressWarnings(if(is.na(survMod[[1]])) return(NULL))
  combinedBins = getBinned(survMod, numBins)
  pvalue = chisq.test(combinedBins)$p.value
  return(pvalue)
}

DCalibrationCumulative = function(listOfSurvivalModels, numBins = 10){
  if(length(listOfSurvivalModels) == 0) return(NULL)
  suppressWarnings(if(any(unlist(lapply(listOfSurvivalModels, is.na)))) return(NA))
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
  if(length(censoredProbabilities) >0){
    censoredBinPositions = sindex(quantiles,censoredProbabilities, comp = "greater", strict = T)
    #Sometimes the probability will be 1 in which case we just want to put them in the first bin. 
    censoredBinPositions = ifelse(censoredBinPositions ==0, 1,censoredBinPositions)
    quantileWidth = 1/numBins
    firstBin = ifelse(censoredBinPositions == numBins,1,(censoredProbabilities - quantiles[censoredBinPositions+1])/censoredProbabilities)
    restOfBins = ifelse(censoredProbabilities ==0,1,1/(numBins*censoredProbabilities))
    
    listOfContributions = lapply(seq_along(censoredBinPositions),function(x) c(rep(0, censoredBinPositions[x] -1),
                                                                               rep(firstBin[x],1),
                                                                               rep(restOfBins[x], numBins - censoredBinPositions[x])))
    censoredBinning = colSums(ldply(listOfContributions,rbind)) 
  } else censoredBinning = 0

  combinedBins = uncensoredBinning + censoredBinning
  names(combinedBins) = quantiles[-length(quantiles)]
  return(combinedBins)
}

