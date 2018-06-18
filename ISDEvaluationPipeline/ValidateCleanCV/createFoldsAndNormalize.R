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
  listOfDatasets = createFoldsOfData(survivalDataset, numberOfFolds)
  listOfImputedDatasets = meanImputation(listOfDatasets, numberOfFolds)
  listOfNormalizedDatasets = normalizeVariables(listOfImputedDatasets, numberOfFolds)
  return(listOfNormalizedDatasets)
}

createFoldsOfData = function(survivalDataset, numberOfFolds){
  #Create folds with the equal amounts of censoring.
  foldIndex = createFolds(as.factor(survivalDataset$delta), k = numberOfFolds, list = TRUE)
  listOfTestingSets = lapply(foldIndex, function(indexs) survivalDataset[indexs,])
  listOfTrainingSets = lapply(foldIndex, function(indexs) survivalDataset[-indexs,])
  listOfDatasets = list(Training = listOfTrainingSets, Testing = listOfTestingSets)
  return(listOfDatasets)
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
  return(listOfDatasets)
}

normalizeVariables = function(listOfImputedDatasets, numberOfFolds){
  for(i in 1:numberOfFolds){
    train = listOfImputedDatasets$Training[[i]]
    test = listOfImputedDatasets$Testing[[i]]
    #We need to remove time an delta so they don't get normalized. Here we set them apart and make an indicator for their location.
    timeDeltaInd = which(names(train) %in% c("time","delta"))
    timeDeltaTrain = train[,timeDeltaInd]
    timeDeltaTest = test[,timeDeltaInd]
    
    #For numeric variables (need to normalize):
    trainNumeric = train[,-timeDeltaInd][,sapply(train[,-timeDeltaInd],is.numeric), drop=FALSE]
    testNumeric = test[,-timeDeltaInd][,sapply(test[,-timeDeltaInd],is.numeric), drop = FALSE]
    varMeans = apply(trainNumeric,2,mean)
    varSD = apply(trainNumeric, 2, sd)
    trainNumericNormalized = apply(trainNumeric, 2, function(x) (x-mean(x))/sd(x))
    #The function is more complex for test because we have to specify the means and sd of the TRAINING set as we go through 
    #the columns of the test set.
    testNumericNormalized = as.data.frame(sapply(1:length(varMeans), function(x) (testNumeric[x] - varMeans[x])/varSD[x]))
    names(testNumericNormalized) = names(testNumeric)
    
    #For the factor variables (need to don one hot encoding):
    trainFactor = train[,-timeDeltaInd][,sapply(train[,-timeDeltaInd],is.factor), drop=FALSE]
    testFactor = test[,-timeDeltaInd][,sapply(test[,-timeDeltaInd],is.factor), drop = FALSE]
    
    #fullRank =T drops one of the extra columns for the one-hot encoding.
    oneHotEncoder = dummyVars("~.",data = trainFactor, fullRank = T)
    trainFactorEncoded = predict(oneHotEncoder, trainFactor)
    testFactorEncoded = predict(oneHotEncoder, testFactor)
    
    #Combine one-hot encoded factor variables and normalized numeric variables and save them as one data frame.
    listOfImputedDatasets$Training[[i]] = cbind.data.frame(timeDeltaTrain,trainNumericNormalized, trainFactorEncoded)
    listOfImputedDatasets$Testing[[i]] = cbind.data.frame(timeDeltaTest,testNumericNormalized, testFactorEncoded)
    
    #We need to make sure that there is no comma in any of the file names or this can wreck functions using csvs.
    #We will search for these and remove them. Ideally we could do this once for the entire dataset but this would require
    #checking all variables and then all the levels of all the factor variables. Since this is linear in the number of 
    #features this check shouldn't be computationally difficult so we will do the "lazy" way and check for commas here 
    #and remove them for every fold.
    namesWithCommas = which(grepl(",",names(listOfImputedDatasets$Training[[i]])))
    if(length(namesWithCommas) > 0){
      names(listOfImputedDatasets$Training[[i]])[namesWithCommas] = sub(",","COMMA",
                                                                        names(listOfImputedDatasets$Training[[i]])[namesWithCommas])
      names(listOfImputedDatasets$Testing[[i]])[namesWithCommas] = sub(",","COMMA",
                                                                        names(listOfImputedDatasets$Testing[[i]])[namesWithCommas])
    }
  }
  return(listOfImputedDatasets)
}





















