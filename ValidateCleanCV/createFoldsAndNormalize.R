#### File Information #####################################################################################################################
#File Name: createFoldsAndNormalize.R
#Date Created: May 25, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#The purpose of this file is to take in a dataset which has already been validated and had some (initial) cleaning. Here we expand upon
#that by creating cross validation folds and doing basic mean imputation and normalizing the data. Specifically, we will do a one hot 
#encoding for factor variables and subtract the mean and divide by the standard deviation for numeric variables. We make sure not to do
#this for the 'time' and 'delta' variables.


### Functions #############################################################################################################################

## Function 1: createFoldsAndNormalize(survivalDataset, numberOfFolds)

# Inputs:
#   survivalDataset: A survival dataset which has been validated and cleaned by validateAndClean().
#   numberOfFolds:   The number of folds to split the data into (must be greater than 1).

# Output: A list of 2 items. The first being the original indexing of the data so that test data is placed back in order at the end.
# The second item is a list containing two lists, a training list and a testing list. Each of these inner lists will contain K datasets
# where K is the number of desired folds.
# Diagram of Output (2nd item):
#                             
#                                  |.--> Training Dataset #1
#                                  |.--> Training Dataset #2
#              |---> Training List |
#              |                   |.--> Training Dataset #(K-1)
#              |                   |.--> Training Dataset #K
# Starting List|
#              |                   |.--> Testing Dataset #1
#              |                   |.--> Testing Dataset #2
#              |--->  Testing List |
#                                  |.--> Testing Dataset #(K-1)
#                                  |.--> Testing Dataset #K

# Usage: Use this function to create folds of data and normalize the features. Then pass those folds to survival models.


## Function 2: createFoldsOfData(survivalDataset, numberOfFolds)

# Inputs: See createFoldsAndNormalize().

# Output: A list containing (1) the fold indicies and (2) a list containing the training/testing folds. See above diagram.

# Usage: Helper function for createFoldsAndNormalize. 
#        This function splits the data into folds. (Helper function of createFoldsAndNormalize).


## Function 3: meanImputation(listOfDatasets)

# Inputs:
#   listOfDatasets: The second output of createFoldOfData (see diagram).

# Output: The same list of datasets passed in but with missing values imputed to the mean of each feature.

# Usage: Helper function for createFoldsAndNormalize. 
# This function imputes all the missing values.


## Function 4: normalizeVariables(listOfImputedDatasets)

# Inputs:
#   listOfImputedDatasetsThe output of meanImputation.

# Output: The same list of datasets passed in but with features normalized.

# Usage: Helper function for createFoldsAndNormalize. 
# This function normalizes all the features.

### Code ##################################################################################################################################
#Library Dependencies:
#We use the one hot encoding function built into caret (dummyVars).
library(caret)
#We use this for the build_scales function to apply normalization to variables.
library(dataPreparation)

createFoldsAndNormalize = function(survivalDataset, numberOfFolds){
  folds = createFoldsOfData(survivalDataset, numberOfFolds)
  originalIndexing = folds[[1]]
  listOfDatasets = folds[[2]]
  listOfImputedDatasets = meanImputation(listOfDatasets)
  listOfNormalizedDatasets = normalizeVariables(listOfImputedDatasets)
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

meanImputation = function(listOfDatasets){
  #Here we take the means of the training and data and impute the training AND the 
  #testing data with the same mean/mode of the respective training data. Not that factor variables are one-hot encoded and then
  #the binary vector is where missing values are imputed as the mean value (of the binary vector).
  for(i in 1:length(listOfDatasets$Training)){
    train = listOfDatasets$Training[[i]]
    test = listOfDatasets$Testing[[i]]
    #Note that we are also imputing time and delta, but validation removed all NA instances of time and delta so we are really imputing 
    #nothing.
    
    #We need to remove time an delta so they don't get one hot encoded (the names would change).
    #Here we set them apart and make an indicator for their location.
    timeDeltaInd = which(names(train) %in% c("time","delta"))
    #fullRank =T drops one of the extra columns for the one-hot encoding.
    oneHotEncoder = dummyVars("~.",data = rbind.data.frame(test,train)[-timeDeltaInd], fullRank = T)
    trainEncoded = predict(oneHotEncoder, train[-timeDeltaInd])
    testEncoded = predict(oneHotEncoder, test[-timeDeltaInd])
    
    varMeans = apply(trainEncoded,2,function(x) mean(x, na.rm = T))
    trainImputed = as.data.frame(apply(trainEncoded, 2, function(x) ifelse(is.na(x),mean(x, na.rm = T),x)))
    
    #The function is more complex for test because we have to specify the means and sd of the TRAINING set as we go through 
    #the columns of the test set.
    testImputed = as.data.frame(sapply(1:length(varMeans), function(x) ifelse(is.na(testEncoded[,x]),varMeans[x],testEncoded[,x])))
    
    #We need to make sure that there is no comma in any of the file names or this can wreck functions using csvs.
    #We will search for these and remove them. Ideally we could do this once for the entire dataset but this would require
    #checking all variables and then all the levels of all the factor variables. Since this is linear in the number of 
    #features this check shouldn't be computationally difficult so we will do the "lazy" way and check for commas here 
    #and remove them for every fold.
    names(trainImputed) = make.names(names(trainImputed), unique = T)
    names(testImputed) = make.names(names(trainImputed), unique = T)
    
    listOfDatasets$Training[[i]] = cbind.data.frame(train[,timeDeltaInd], trainImputed)
    listOfDatasets$Testing[[i]] = cbind.data.frame(test[,timeDeltaInd], testImputed)
  }
  return(listOfDatasets)
}

normalizeVariables = function(listOfImputedDatasets){
  listOfNormalizedDatasets = list()
  for(i in 1:length(listOfImputedDatasets$Training)){
    train = listOfImputedDatasets$Training[[i]]
    test = listOfImputedDatasets$Testing[[i]]
    #We need to remove time an delta so they don't get normalized. Here we set them apart and make an indicator for their location.
    timeDeltaInd = which(names(train) %in% c("time","delta"))
    timeDeltaTrain = train[,timeDeltaInd]
    timeDeltaTest = test[,timeDeltaInd]
    
    scales = build_scales(train[-timeDeltaInd], verbose = F)
    zeroVarianceFeatures = which(unlist(scales)[seq(2,length(unlist(scales)),by=2)]== 0)
    normalizedTrain = fast_scale(train[-timeDeltaInd], scales=scales, verbose=F)
    normalizedTest = fast_scale(test[-timeDeltaInd], scales=scales,verbose=F)
    #Note that if there is zero variance in a feature of the training set the models *should* be assigning it a weight of zero, thus
    #the variable shouldn't matter. Regardless, we turn these values back into an "unormalized" version -- note the normalized version is 
    #NaN since the sd would be zero, i.e. we are dividing by 0. Dividing by 0 is bad for your health and safety.
    if(length(zeroVarianceFeatures) > 0){
      normalizedTrain[,zeroVarianceFeatures] = train[-timeDeltaInd][,zeroVarianceFeatures]
      normalizedTest[,zeroVarianceFeatures] = test[-timeDeltaInd][,zeroVarianceFeatures]
    }
    trainCentered = cbind.data.frame(timeDeltaTrain, normalizedTrain)
    testCentered = cbind.data.frame(timeDeltaTest,normalizedTest)
    
    #If a variable had zero variance in the test set we will end up dividing by 0 and getting NaN values.
    #Here we simply make everything 0, the mean value, if this occurs.
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
    }

  return(listOfNormalizedDatasets)
}





















