#File Name: validateAndClean.R
#Date: May 25, 2018
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
#No library calls needed.

#Functions:

validateAndClean = function(survivalDataset){
  validatedData = validate(survivalDataset)
  cleanData = clean(validatedData)
  return(cleanData)
}

validate = function(survivalDataset){
  if(is.null(survivalDataset) || nrow(survivalDataset) == 0)
    stop("You must enter a nonempty survival dataset.")
  if(!"time" %in% names(survivalDataset))
    stop("The variable 'time' in not included in the given dataset.")
  if(!"delta" %in% names(survivalDataset))
    stop("The variable 'delta' in not included in the given dataset.")
  survivalDataset$time = as.numeric(survivalDataset$time)
  survivalDataset$delta = as.numeric(survivalDataset$delta)
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
  #Remove columns with over 25% of data missing.
  naProp = apply(survivalDataset, 2, function(x) sum(is.na(x))/nrow(survivalDataset))
  toRemove = which(naProp > 0.25) 
  if(length(toRemove)>0){
    namesToRemove = names(survivalDataset)[toRemove]
    survivalDataset = survivalDataset[,-toRemove]
    warning(paste("The variable(s)",paste(namesToRemove, collapse = ", "),"contained over 25% NA values. They have been removed.", sep = " "))
  }
  return(survivalDataset)
}

