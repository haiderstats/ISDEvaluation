#File Name: validateAndClean.R
#Date Created: May 25, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#The purpose of this file is to take in a survival dataset and validate that it has all parameters required for analysis. Specifically, 
#we need to validate that the dataset has two columns names time and delta which are our time and event indicator variables respectively.
#We will also check that these variables do not have things such as nonnegative times, and validate that delta only takes on values of
#0 and 1, i.e. right censor and event indicators - we do not plan to handle any other type of censoring. From there we will do basic
#cleaning procedures such as removing heavily missing columns (greater than 25% of observations).

#Input: Survival Dataset
#Output: Survival Dataset which has been validated and had minor cleaning.
############################################################################################################################################
#We call caret for the nearZeroVar function. This function will let us remove variables with all the same value.
library(caret)
#Functions:

validateAndClean = function(survivalDataset, imputeZero=T){
  validatedData = validate(survivalDataset, imputeZero=T)
  cleanData = clean(validatedData)
  return(cleanData)
}

validate = function(survivalDataset,imputeZero){
  if(is.null(survivalDataset) || nrow(survivalDataset) == 0)
    stop("You must enter a nonempty survival dataset.")
  if(!"time" %in% names(survivalDataset))
    stop("The variable 'time' in not included in the given dataset.")
  if(!"delta" %in% names(survivalDataset))
    stop("The variable 'delta' in not included in the given dataset.")
  survivalDataset$time = as.numeric(survivalDataset$time)
  survivalDataset$delta = as.numeric(survivalDataset$delta)
  #We want to change the ordering so that the time column is first and the delta column is second.
  timeIndex = which(names(survivalDataset) == "time")
  deltaIndex = which(names(survivalDataset) == "delta")
  survivalDataset = survivalDataset[,c(timeIndex,deltaIndex, (1:ncol(survivalDataset))[-c(timeIndex,deltaIndex)])]
  if(any(is.na(survivalDataset$time))){
    warning("'time' includes NA values. These rows have been removed.")
    NAtimeIndex = which(is.na(survivalDataset$time))
    survivalDataset = survivalDataset[-NAtimeIndex,]
  }
  if(any(is.infinite(survivalDataset$time))){
    warning("'time' includes Inf values. These rows have been removed.")
    InftimeIndex = which(is.infinite(survivalDataset$time))
    survivalDataset = survivalDataset[-InftimeIndex,]
  }
  if(any(survivalDataset$time <0)){
    warning("'time' includes negative values. These rows have been removed.")
    nonPosTimeIndex = which(survivalDataset$time <0)
    survivalDataset = survivalDataset[-nonPosTimeIndex,]
  }
  if(any(survivalDataset$time ==0) & imputeZero){
    imputeVal  =min(survivalDataset$time[survivalDataset$time > 0])/2
    warning(paste("'time' includes 0 valued times. These have been imputed to half the minimum positive time point:",
                  imputeVal))
    zeroTimeIndex = which(survivalDataset$time ==0)
    survivalDataset[zeroTimeIndex,]$time = imputeVal
  }
  if(any(is.na(survivalDataset$delta))){
    warning("'delta' includes NA values. These rows have been removed.")
    NAdeltaIndex = which(is.na(survivalDataset$delta))
    survivalDataset = survivalDataset[-NAdeltaIndex,]
  }
  if(any(is.infinite(survivalDataset$delta))){
    warning("'delta' includes Inf values. These rows have been removed.")
    InfdeltaIndex = which(is.infinite(survivalDataset$delta))
    survivalDataset = survivalDataset[-InfdeltaIndex,]
  }
  if(length(unique(survivalDataset$delta)) >2)
    stop("'delta' contains more than two unique values. 'delta' should only contain 0 and 1.")
  if(!all(sort(unique(survivalDataset$delta)) == c(0,1)))
    stop("'delta' contains values other than 0 and 1. 0 should indicate RIGHT censoring and 1 should indicate the event occuring.")
  return(survivalDataset)
}

clean = function(survivalDataset){
  #Remove columns with no variance.
  allSame = nearZeroVar(survivalDataset,freqCut = 100/0)
  if(length(allSame)>0){
    namesToRemoveNoVar = names(survivalDataset)[allSame]
    survivalDataset = survivalDataset[,-allSame]
    warning(paste("The variable(s)",paste(namesToRemoveNoVar, collapse = ", "),"contained only 1 unique value. They have been removed.",
                  sep = " "))
  }
  #Remove columns with over 25% of data missing.
  naProp = apply(survivalDataset, 2, function(x) sum(is.na(x))/nrow(survivalDataset))
  toRemoveNA = which(naProp > 0.25)
  if(length(toRemoveNA)>0){
    namesToRemoveNA = names(survivalDataset)[toRemoveNA]
    survivalDataset = survivalDataset[,-toRemoveNA]
    warning(paste("The variable(s)",paste(namesToRemoveNA, collapse = ", "),"contained over 25% NA values. They have been removed.",
                  sep = " "))
  }
  return(survivalDataset)
}

