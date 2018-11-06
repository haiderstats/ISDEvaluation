#### File Information #####################################################################################################################
#File Name: Concordance.R
#Date Created: May 28th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#This file was created to implement concordance as an evaluation measure for individual survival curves. 

### Functions #############################################################################################################################

## Function 1: Concordance(survMod, ties = "None", method = "Median")

#Inputs:
#   survMod: A list of 4 items:(1) TestCurves - The survival curves for the testing set.
#                              (2) TestData - The censor/death indicator and event time for the testing set. 
#                              (3) TrainData - The censor/death indicator and event time for the training set. 
#                              (4) TrainCurves - The survival curves for the training set.
#   ties: A string indicating the way ties should be handled. Options: "None" will throw out all ties in survival time and all ties from
#          risk scores. "Time" includes ties in survival time but removes ties in risk scores. "Risk" includes ties in risk scores but 
#          not in survival time. "All" includes all ties (both in survival time and in risk scores). Note the concordance calculation is
#          given by (Concordant Pairs + (Number of Ties/2))/(Concordant Pairs + Discordant Pairs + Number of Ties).
#   method: A string indicating whether the "Mean" or "Median" should be used to calculate a patient's risk score.

# Output: The C-index.

# Usage: Calculate Concordance given a survival model.

### Code ##################################################################################################################################
#Library Dependencies:
#We use this for the survConcordance function.
library(survival)

#Helper Functions: predictMeanSurvivalTimeSpline(survivalCurve,predictedTimes) (Or median)
source("Evaluations/EvaluationHelperFunctions.R")

#The following function is split into 2 parts. Part 1 retrieves all the relevant pieces from the passed in survMod object, e.g. the survival
#curves and the true death times of test subjects. Part 2 uses survConcordance to calculate classical concordance measures. Additionally, this
#is were tied data is handled.
Concordance = function(survMod, ties = "None", method = "Median"){
  #Part 1:
  #Being passed an empty model.
  if(is.null(survMod)) return(NULL)
  #Being passed a model that failed.
  suppressWarnings(if(is.na(survMod[[1]])) return(NULL))
  if(!ties %in% c("None","Risk","Time","All"))
    stop("Please enter one of: 'None', 'Risk','Time', or 'All' as the ties argument.")
  predictedTimes = survMod[[1]][,1]
  survivalCurves = survMod[[1]][-1]
  trueDeathTimes = survMod[[2]][,1]
  censorStatus = survMod[[2]][,2]
  
  predictMethod = switch(method,
                         Mean = predictMeanSurvivalTimeSpline,
                         Median = predictMedianSurvivalTimeSpline)
  #This retrieves the mean death probability of the survival curve.
  averageSurvivalTimes = unlist(lapply(seq_along(trueDeathTimes),
                                       function(index) predictMethod(survivalCurves[,index],
                                                                                     predictedTimes)))
  #Part 2:
  #The risk score should be higher for subjects that live shorter (i.e. lower average survival time).
  risk = -1*averageSurvivalTimes
  concordanceInfo = survConcordance(Surv(trueDeathTimes, censorStatus)~ risk)
  concordantPairs= concordanceInfo$stats[1] 
  discordantPairs = concordanceInfo$stats[2] 
  riskTies = concordanceInfo$stats[3]
  timeTies = concordanceInfo$stats[4]
  
  CIndex = switch(ties,
                  None = concordantPairs/(concordantPairs + discordantPairs),
                  Time = (concordantPairs+timeTies/2)/(concordantPairs + discordantPairs + timeTies),
                  Risk = (concordantPairs+riskTies/2)/(concordantPairs + discordantPairs + riskTies),
                  All = (concordantPairs+(riskTies +timeTies)/2)/(concordantPairs + discordantPairs + timeTies + riskTies)
  )
  return(CIndex)
}
