#### File Information #####################################################################################################################
#File Name: AcceleratedFailureTime.R
#Date Created: May 28, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################

#This file is used to run the Accelerated Failure Time (AFT) model for some specified distribution to create individual survival curves.
#This implimentaiton only supports the following distributions: weibull, exponential, lognormal, gaussian, loglogistic, and logistic.

### Functions #############################################################################################################################

## Function 1: AFT(training, testing, AFTDistribution)

#Inputs:
#   training: The training dataset (after normalization and imputation).
#   testing: The testing dataset (after normalization and imputation).
#   AFTDistribution: The distribution to use for the AFT model.

# Output: A list of 4 items:(1) TestCurves - The survival curves for the testing set.
#                           (2) TestData - The censor/death indicator and event time for the testing set. 
#                           (3) TrainData - The censor/death indicator and event time for the training set. 
#                           (4) TrainCurves - The survival curves for the training set.

# Usage: Train and evaluate the AFT model.


## Function 2:  survfunc(AFTMod, newdata, t)

# Inputs:
#   AFTMod: The AFT model returned by survreg().
#   newdata: The data on which to evaluate the AFT model.
#   t: The times at which to evaluate the AFT model.

# Output: The survival curves for each patient.

# Usage: Evaluate the AFT model to get survival probabilites. (Helper function for AFT).


### Code #############################################################################################################################
#Library Dependencies:
#survival is needed to get survreg for the AFT model.
library(survival)

AFT = function(training, testing, AFTDistribution){
  #Sometimes the AFT model will fail to converge (though rare) we want to catch this.
  tryCatch({
    AFTMod = survreg(Surv(time,delta)~., data = training, dist = AFTDistribution)
    
    trainingTimes = sort(unique(training$time))
    if(0 %in% trainingTimes){
      timesToPredict = trainingTimes
    } else {
      timesToPredict = c(0,trainingTimes)
    }
    survivalCurves = survfunc(AFTMod, newdata = testing, t = timesToPredict)
    survivalCurvesTrain = survfunc(AFTMod, newdata = training, t = timesToPredict)
  },
  error = function(e) {
    message(e)
    warning("AFT failed to converge.")
  })
  if(!exists("AFTMod") | !exists("survivalCurves")){
    return(NA)
  }

  probabilities = survivalCurves$sur
  probabilitiesTrain = survivalCurvesTrain$sur
  #Since survfunc returns survival probabilities with the first time point for every individual (ordered by how the testing individuals)
  #were passed in, we can simply fill a matrix by row to have each individual curve be a column. This can be verified by checking
  #the survival probabilties (sur) for any ID_SurvivalCurves against any column, e.g. the survival probabilities for ID_SurvivalCurves == 2,
  #correspond to the probabilites found in the second column of the matrix below.
  probabilityMatrix = matrix(probabilities, ncol = nrow(testing),byrow = T)
  probabilityTrainMatrix = matrix(probabilitiesTrain, ncol = nrow(training),byrow = T)
  
  curvesToReturn = cbind.data.frame(time = timesToPredict, probabilityMatrix)
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  trainingCurvesToReturn = cbind.data.frame(time = timesToPredict, probabilityTrainMatrix)
  return(list(TestCurves = curvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

#The following was taken and altered from http://rstudio-pubs-static.s3.amazonaws.com/161203_6ee743eb28df4cd68089a519aa705123.html.
#This code is used to pass in a time point and get the predicted probability. The predict function for survreg objects only return 
#times from probabilities so we need to reverse enginner this using the code below.
survfunc = function (object, t, newdata, name = "t") {
  #Altered from origina: I am going to add an ID to every row so we can retrieve the individuals easily from the output.
  #I gave a weird ID variable name so if the original data came in with a variable ("ID") it won't break our system.
  newdata$ID_SurvivalCurves = 1:nrow(newdata)
  newdata <- do.call(rbind, rep(list(newdata), length(t)))
  t <- rep(t, each = nrow(newdata)/length(t))
  if (class(object) != "survreg") 
    stop("not a survreg object")
  lp <- predict(object, newdata = newdata, type = "lp")
  if (object$dist %in% c("weibull", "exponential")) {
    newdata$pdf <- dweibull(t, 1/object$scale, exp(lp))
    newdata$cdf <- ifelse(t == 0,0,
                          ifelse(is.nan(pweibull(t, 1/object$scale, exp(lp))),1,pweibull(t, 1/object$scale, exp(lp))))
    newdata$haz <- exp(dweibull(t, 1/object$scale, exp(lp), 
                                log = TRUE) - pweibull(t, 1/object$scale, exp(lp), 
                                                       lower.tail = FALSE, log.p = TRUE))
  }
  else if (object$dist == "lognormal") {
    newdata$pdf <- dlnorm(t, lp, object$scale)
    newdata$cdf <- plnorm(t, lp, object$scale)
    newdata$haz <- exp(dlnorm(t, lp, object$scale, log = TRUE) - 
                         plnorm(t, lp, object$scale, lower.tail = FALSE, log.p = TRUE))
  }
  else if (object$dist == "gaussian") {
    newdata$pdf <- dnorm(t, lp, object$scale)
    newdata$cdf <- pnorm(t, lp, object$scale)
    newdata$haz <- exp(dnorm(t, lp, object$scale, log = TRUE) - 
                         pnorm(t, lp, object$scale, lower.tail = FALSE, log.p = TRUE))
  }
  else if (object$dist == "loglogistic") {
    newdata$pdf <- dlogis(log(t), lp, object$scale)/t
    newdata$cdf <- plogis(log(t), lp, object$scale)
    newdata$haz <- exp(dlogis(log(t), lp, object$scale, log = TRUE) - 
                         log(t) - plogis(log(t), lp, object$scale, lower.tail = FALSE, 
                                         log.p = TRUE))
  }
  else if (object$dist == "logistic") {
    newdata$pdf <- dlogis(t, lp, object$scale)
    newdata$cdf <- plogis(t, lp, object$scale)
    newdata$haz <- exp(dlogis(t, lp, object$scale, log = TRUE) - 
                         dlogis(t, lp, object$scale, lower.tail = FALSE, log.p = TRUE))
  }
  else {
    stop("unknown distribution")
  }
  newdata$sur <- 1 - newdata$cdf
  newdata[name] <- t
  return(newdata)
}

