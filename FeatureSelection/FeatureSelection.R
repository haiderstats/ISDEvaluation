#### File Information #####################################################################################################################
#File Name: FeatureSelection.R
#Date Created: June 21st, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#This file allows for feature selection in the pre processing steps for ISD evaluations. Future feature selection methods should be written
#and added to the switch statement.

### Functions #############################################################################################################################

## Function 1: FeatureSelection(dataset, type ="UniCox", obs_thresh = 0, pThresh = 0.1)

# Inputs:
#   dataset:    The survival dataset post validation.
#   type:       The type of feature selection. Currently only "UniCox" is supported.
#   obs_thresh: The number of observations that must be included in that feature (e.g. for categorical featrures) to be considered.
#   pThresh:    The p-value threshold required to keep the feature (univariate cox model p-value must be lower than threshold).

# Output: The reduced survival dataset (also oneHot encoded).

# Usage: Use this function to perform feature selection on the dataset.


## Function 2: uniCox(dataset, obs_thresh, pThresh)

# Inputs: See FeatureSelection().

# Output: The reduced survival dataset.

# Usage: This performs univariate cox feature selection.

### Code ##################################################################################################################################
#Library Dependencies:
#We use caret for the dummyVars function.
library(caret)
FeatureSelection = function(dataset, type ="UniCox", obs_thresh = 0, pThresh = 0.1){
  selectedData = switch(type,
                        UniCox = {
                         oneHotEncoder = dummyVars("~.",data = dataset, fullRank = T)
                         toSelect = as.data.frame(predict(oneHotEncoder, dataset))
                         names(toSelect) = make.names(names(toSelect), unique = T)
                         uniCox(toSelect, obs_thresh, pThresh)
                        }
  )
}

#The following function was developed by Bret Hoehn (bhoehn@ualberta.ca) at some date and given to Humza Haider for doing feature selection.
#Humza then made some modifications from the original and this modified version is given below.
#Here we are simply seeing if variables are significant in a univariate cox model. If the p-value is below a specified value (default is 0.1),
#Then it is returned in the dataset, otherwise it is removed.
uniCox <- function(dataset, obs_thresh, pThresh){
  # select variables with p-val <pThresh via univariate Cox

  # determine the threshold for whether or not a variable is acceptable (each allowed value should correspond to a
  #minimum number of observations)
  num_obs = sum(dataset$delta)
  num_obs_thresh = num_obs*obs_thresh
  
  min_val <- pThresh
  min_var <- ""
  chosen_vars <- c()   
  var_names = names(dataset)[-which(names(dataset) %in% c("time","delta"))]
  for (var_name in var_names){
    formula = as.formula(paste("Surv(time, delta) ~", var_name))
    unique_vals = unique(dataset[[var_name]])
    cox.fit <- coxph(formula, data = dataset, singular.ok = TRUE)
    
    min_obs_restriction_met=TRUE
    if (length(unique_vals)<=3){
      # check if each setting of this variable satisfies the required number of observed events
      for (val in unique_vals){
        dset_valIdx = (dataset[[var_name]]==val)
        num_obs_val = sum(dataset$delta[dset_valIdx],na.rm=T)
        if (num_obs_val < num_obs_thresh)
          min_obs_restrictionMet=FALSE
      }
    }
    
    #get the p-value of this feature
    if (min_obs_restriction_met==FALSE){
        print(paste(var_name," does not have the required number of observations for each value setting"))
        pValue <- 1
      }else {
        pValue <- summary(cox.fit)$coefficients[,5]
      }
    if (is.nan(pValue) | is.na(pValue)) {
      pValue <- 1
    }
    if (pValue < min_val) {
      min_val <- pValue
      min_var <- var_name
    }
    
    if (!is.nan(pValue) && pValue < 0.1 & !is.na(pValue)){
      chosen_vars <- c(chosen_vars, var_name)
    }
  }
  if (length(chosen_vars) == 0) {
    print(paste("All features eliminated by univariate selection, including ", min_var, " by default"))
    chosen_vars <- c(min_var)
  }
  chosenVarIndex = which(names(dataset) %in% c("time","delta",chosen_vars))
  return (dataset[,chosenVarIndex])
}
