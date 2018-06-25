#File Name: createFoldsAndNormalize.R
#Date Created: May 25, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#The purpose of this file is to take in a dataset which has already been validated and had some (initial) cleaning. Here we expand upon
#that by creating cross validation folds and doing basic mean imputation and normalizing the data. Specifically, we will do a one hot 
#encoding for factor variables and subtract the mean and divide by the standard deviation for numeric variables. We make sure not to do
#this for the 'time' and 'delta' variables.

#Input #1: Survival Dataset which has been validated and had initial cleaning
#Input #2: Desired number of folds.
#
#Output: List containing two lists, a training list and a testing list. Each of these inner lists will contain K datasets where K is the
#number of desired folds.
#Diagram of Output:
#                             
#                                  |.--> Training Dataset #1
#                                  |.--> Training Dataset #2
#              |---> Training List |
#              |                   |.--> Training Dataset #(K-1)
#              |                   |.--> Training Dataset #K
#Starting List |
#              |                   |.--> Testing Dataset #1
#              |                   |.--> Testing Dataset #2
#              |--->  Testing List |
#                                  |.--> Testing Dataset #(K-1)
#                                  |.--> Testing Dataset #K
############################################################################################################################################
#Library dependencies:
#caret is a common R machine learning package. Specifically, I call it here to avoid implementing my own cross validation function. 
#Specifically, I want to harness the stratefied cross validation as we would like to have the same amount of censoring in each fold.
#Additionally, we use the one hot encoding function built into caret (dummyVars).
library(caret)

#Functions:
createFoldsAndNormalize = function(survivalDataset, numberOfFolds){
  folds = createFoldsOfData(survivalDataset, numberOfFolds)
  originalIndexing = folds[[1]]
  listOfDatasets = folds[[2]]
  listOfImputedDatasets = meanImputation(listOfDatasets, numberOfFolds)
  listOfNormalizedDatasets = normalizeVariables(listOfImputedDatasets, numberOfFolds)
  return(list(originalIndexing, listOfNormalizedDatasets))
}

createFoldsOfData = function(survivalDataset, numberOfFolds){
  #Create folds with the equal amounts of censoring and quasi equal ranges. Here we order by censor status and time to event
  #and then match these off to folds one by one. Doing it this way makes cross validation deterministic which is not ideal
  #but lets us have equally matched folds.
  time = survivalDataset$time
  delta = survivalDataset$delta
  Order= order(delta,time)
  foldIndex = lapply(1:numberOfFolds, function(x) Order[seq(x,nrow(survivalDataset), by = numberOfFolds)])
  listOfTestingSets = lapply(foldIndex, function(indexs) survivalDataset[indexs,])
  listOfTrainingSets = lapply(foldIndex, function(indexs) survivalDataset[-indexs,])
  listOfDatasets = list(Training = listOfTrainingSets, Testing = listOfTestingSets)
  return(list(foldIndex,listOfDatasets))
}

meanImputation = function(listOfDatasets, numberOfFolds){
  #Here we take the means (numeric variables) and modes (factor variables) of the training and data and impute the training AND the 
  #testing data with the same mean/mode of the respective training data.
  for(i in 1:numberOfFolds){
    train = listOfDatasets$Training[[i]]
    test = listOfDatasets$Testing[[i]]
    #Note that we are also imputing time and delta, but validation removed all NA instances of time and delta so we are really imputing 
    #nothing.
    
    #For numeric variables:
    trainNumeric = train[,sapply(train,is.numeric), drop = FALSE]
    #drop = FALSE makes it stay as a dataframe instead of turning to a vector in the event there is only 1 variable with the condition.
    testNumeric = test[,sapply(test,is.numeric), drop = FALSE]
    varMeans = apply(trainNumeric,2,function(x) mean(x, na.rm = T))
    trainNumericImputed = apply(trainNumeric, 2, function(x) ifelse(is.na(x),mean(x, na.rm = T),x))
    #The function is more complex for test because we have to specify the means and sd of the TRAINING set as we go through 
    #the columns of the test set.
    testNumericImputed = as.data.frame(sapply(1:length(varMeans), function(x) ifelse(is.na(testNumeric[,x]),varMeans[x],testNumeric[,x])))
    #Names are getting lost in sapply
    names(testNumericImputed) = names(testNumeric)
    
    
    #For factor variables:
    trainFactor = train[,sapply(train,is.factor),drop=FALSE]
    testFactor = test[,sapply(test,is.factor),drop=FALSE]    
    if(ncol(trainFactor) > 0){
      #We need to introduce a mode function to find the mode of the factor variables.
      Mode <- function(x) {
        #Modified from https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
        x = x[!is.na(x)]
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
      }
      varModes = apply(trainFactor,2,function(x) Mode(x))
      trainFactorListed = lapply(trainFactor, as.character)
      trainFactorImputed = as.data.frame(lapply(trainFactorListed, function(x) ifelse(is.na(x), Mode(x), x)))
      #We run into problems if the training dataset didn't have a factor level in it which was included in the 
      #original dataset. Here we are just adding back the original factor levels.
      for(j in 1:ncol(trainFactorImputed)){
        missedLevels = levels(trainFactor[,j])[which(!levels(trainFactor[,j]) %in% levels(trainFactorImputed[,j]))]
        levels(trainFactorImputed[,j]) = c(levels(trainFactorImputed[,j]), missedLevels)
      }
      testFactorImputed = as.data.frame(sapply(1:length(varModes),
                                               function(x) factor(ifelse(is.na(testFactor[,x]),
                                                                         varModes[x],
                                                                         paste(testFactor[,x])),
                                                                  levels = levels(trainFactorImputed[,x]))))
      names(testFactorImputed) = names(testFactor)
      #Combine numeric and factor variables and save over the previous datasets.
      listOfDatasets$Training[[i]] = cbind.data.frame(trainNumericImputed, trainFactorImputed)
      listOfDatasets$Testing[[i]] = cbind.data.frame(testNumericImputed, testFactorImputed)
    }
   else{
     listOfDatasets$Training[[i]] = as.data.frame(trainNumericImputed)
     listOfDatasets$Testing[[i]] = as.data.frame(testNumericImputed)
   }
  }
  return(listOfDatasets)
}

normalizeVariables = function(listOfImputedDatasets, numberOfFolds){
  listOfNormalizedDatasets = list()
  for(i in 1:numberOfFolds){
    train = listOfImputedDatasets$Training[[i]]
    test = listOfImputedDatasets$Testing[[i]]
    #We need to remove time an delta so they don't get normalized. Here we set them apart and make an indicator for their location.
    timeDeltaInd = which(names(train) %in% c("time","delta"))
    timeDeltaTrain = train[,timeDeltaInd]
    timeDeltaTest = test[,timeDeltaInd]
    
    #fullRank =T drops one of the extra columns for the one-hot encoding.
    oneHotEncoder = dummyVars("~.",data = train[-timeDeltaInd], fullRank = T)
    trainEncoded = predict(oneHotEncoder, train[-timeDeltaInd])
    testEncoded = predict(oneHotEncoder, test[-timeDeltaInd])
    
    trainCentered = cbind.data.frame(timeDeltaTrain,scale(trainEncoded))
    testCentered = cbind.data.frame(timeDeltaTest,scale(testEncoded))
    
    trainCentered = apply(trainCentered,2, function(x){
      if(all(is.nan(x))) 
        return(rep(0,length(x)))
      return(x)
    })
    testCentered = apply(testCentered,2, function(x){
      if(all(is.nan(x))) 
        return(rep(0,length(x)))
      return(x)
    })

    listOfNormalizedDatasets$Training[[i]] = as.data.frame(trainCentered)
    listOfNormalizedDatasets$Testing[[i]] = as.data.frame(testCentered)
    #We need to make sure that there is no comma in any of the file names or this can wreck functions using csvs.
    #We will search for these and remove them. Ideally we could do this once for the entire dataset but this would require
    #checking all variables and then all the levels of all the factor variables. Since this is linear in the number of 
    #features this check shouldn't be computationally difficult so we will do the "lazy" way and check for commas here 
    #and remove them for every fold.
    names(listOfNormalizedDatasets$Training[[i]]) = make.names(names(listOfNormalizedDatasets$Training[[i]]),unique=T)
    names(listOfNormalizedDatasets$Testing[[i]]) = make.names(names(listOfNormalizedDatasets$Testing[[i]]),unique=T)
    }

  return(listOfNormalizedDatasets)
}





















