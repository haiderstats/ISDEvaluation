#File Name: Concordance.R
#Date Created: May 28th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file was created to implement concordance as an evaluation measure for individual survival curves. 
#Input 1: A list of (1) a matrix of survival curves, (2) the true death times and event/censoring indicator (delta =1 implies death/event) of 
#the TESTING data, and (3) the true death times and event/censoring indicator (delta =1 implies death/event) of the TRAINING data
#Input 2: A string indicating the way ties should be handled. Options: "None" will throw out all ties in survival time and all ties from
#risk scores. "Time" includes ties in survival time but removes ties in risk scores. "Risk" includes ties in risk scores but not in survival
#time. "All" includes all ties (both in survival time and in risk scores). Note the concordance calculation is given by
#(Concordant Pairs + (Number of Ties/2))/(Concordant Pairs + Discordant Pairs + Number of Ties)
#Currently a "risk" score for a survival distribution model is considered to be the mean survival time using a spline fit to make a continuous
#survival distribution. 
#Input 3: A boolean on whether or not to include censored patients in the calculation of concordance. See the paper for details on how 
#concordance is calculated for censored individuals.
#Output: The C-index.
##############################################################################################################################################
#Library Dependencies
#We use this for the prodlim function to make a Kaplan-Meier curve. Specifically, we use prodlim due to its predict function.
library(prodlim)
#We use this for the survConcordance function.
library(survival)
#We use this for ldply, a combiner of lists.
library(plyr)
#Helper Functions: predictMeanSurvivalTimeLinear(survivalCurve,predictedTimes)
source("Evaluations/EvaluationHelperFunctions.R")

#The following function is split into 3 parts. Part 1 retrieves all the relevant pieces from the passed in survMod object, e.g. the survival
#curves and the true death times of test subjects. Part 2 is only used if we are to include concordance calculations for two censored patients,
#see the paper for the details on this calculation. Part 3 uses survConcordance to calculate classical concordance measures. We then 
#add in the censored piece from part 2 if includeCensored == True. Additionally, this is were tied data is handled.
Concordance = function(survMod, ties = "None",includeCensored =F){
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
  trainingDeathTimes = survMod[[3]]$time
  trainingCensorStatus = survMod[[3]]$delta
  
  #This retrieves the mean death probability of the survival curve.
  averageSurvivalTimes = unlist(lapply(seq_along(trueDeathTimes),
                                       function(index) predictMeanSurvivalTimeSpline(survivalCurves[,index],
                                                                                     predictedTimes)))
  #Part 2:
  if(includeCensored){
    #Here we make the Kaplan-Meier curve and alter the predict function to be linear outside of the last time point of the Kaplan-Meier 
    #curve.
    KMCurve = prodlim(Surv(trainingDeathTimes, trainingCensorStatus)~1)
    KMLinearPredict = function(time){
      prediction = predict(KMCurve,time)
      slope = (1-min(KMCurve$surv))/(0 - max(KMCurve$time))
      predictedProbabiliteis = ifelse(is.na(prediction), pmax(1+time*slope,0), prediction)
      return(predictedProbabiliteis)
    }
    
    orderDeathTimes = order(trueDeathTimes)
    sortedDeathTimes = sort(trueDeathTimes)
    #Each row of the following matrix represents a single individual. Each column is the (ith +1) individual in the testing data.
    #A value of 1 in the ith row and jth column means that the ith individual (sorted by event time/censored time) is censored and
    #has a smaller event time than the jth +1 individual. Note that we ignore ties by using sortedDeathTimes[x] < sortedDeathTimes[-x]
    #as opposed to sortedDeathTimes[x] <= sortedDeathTimes[-x].
    
    trueCensMatrix = ldply(lapply(seq_along(sortedDeathTimes), function(x) as.numeric(sortedDeathTimes[x] < sortedDeathTimes[-x]  &
                                                                                        !censorStatus[orderDeathTimes][x])),rbind)
    
    rowIndex = 0
    estimatedCensMatrix = apply(trueCensMatrix,1,function(x){
      rowIndex <<- rowIndex+1
      #Need to add 1 since there are number of individuals -1 columns.
      toCompare = which(x == 1) +1
      uncensoredInd = as.logical(censorStatus[orderDeathTimes][toCompare])
      if(length(toCompare)>0){
        s2 = KMLinearPredict(sortedDeathTimes[toCompare])
        s1 = rep(KMLinearPredict(sortedDeathTimes[rowIndex]), length(toCompare))
        greaterProb = ifelse(uncensoredInd,1-(s2/s1),1 - 0.5*(s2/s1))
        #If s1 evaluates to 0 then the greater probability should be 1 - they should both already be dead anyway.
        greaterProb = ifelse(is.nan(greaterProb),1,greaterProb)
        predictedProb = ifelse(averageSurvivalTimes[orderDeathTimes][rowIndex] <= averageSurvivalTimes[orderDeathTimes][toCompare],
                               greaterProb,1-greaterProb)
        x[toCompare-1] = predictedProb
      }
      return(x)
    })
  }
  #Part 3:
  #The risk score should be higher for subjects that live shorter (i.e. lower average survival time).
  risk = -1*averageSurvivalTimes
  concordanceInfo = survConcordance(Surv(trueDeathTimes, censorStatus)~ risk)
  concordantPairs= concordanceInfo$stats[1] + ifelse(includeCensored,sum(estimatedCensMatrix),0)
  discordantPairs = concordanceInfo$stats[2] + ifelse(includeCensored,sum(trueCensMatrix) - sum(estimatedCensMatrix),0)
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








