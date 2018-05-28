#File Name: AcceleratedFailureTime.R
#Date Created: May 28, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file is used to run the Accelerated Failure Time (AFT) model for some specified distribution. The most common distributions for
#AFT are the Weibull, Log-logistic, and Log-Normal distributions. This implimentaiton only supports the following distributions:
#weibull, exponential, lognormal, gaussian, loglogistic, and logistic.
#Here we will take in a training and a testing set. The function will then train on the training set and return a list containing (1) a 
#matrix of survival curves, where the first column is time values and all following columns are survival probabilities of test subjects and 
#(2) the true death times and event indicator (i.e. time and delta) of the test subjects.

#Input #1: Survival Dataset post normalization and imputation.
#Input #2: The distribution for AFT to use.
#Output: A list of (1) matrix of survival curves, and (2) the true death times.
############################################################################################################################################
#Library Dependencies
#survival is needed to get survreg for the AFT model.
library(survival)

AFT = function(training, testing, AFTDistribution){
  aftMod = survreg(Surv(time,delta)~., data = training, dist = AFTDistribution)
  timesToPredict = c(0,sort(unique(training$time)))
  survivalCurves = survfunc(aftMod, newdata = testing, t = timesToPredict)
  probabilities = survivalCurves$sur
  #Since survfunc returns survival probabilities with the first time point for every individual (ordered by how the testing individuals)
  #were passed in, we can simply fill a matrix by row to have each individual curve be a column. This can be verified by checking
  #the survival probabilties (sur) for any ID_SurvivalCurves against any column, e.g. the survival probabilities for ID_SurvivalCurves == 2,
  #correspond to the probabilites found in the second column of the matrix below.
  probabilityMatrix = matrix(probabilities, ncol = nrow(testing),byrow = T)
  curvesToReturn = cbind.data.frame(time = timesToPredict, probabilityMatrix)
  timesAndCens = cbind.data.frame(time = testing$time, delta = testing$delta)
  return(list(curvesToReturn, timesAndCens))  
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
    newdata$cdf <- pweibull(t, 1/object$scale, exp(lp))
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

