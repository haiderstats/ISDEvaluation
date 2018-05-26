#File Name: coxPH_KP.R
#Date Created: May 26, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file is used to run the Cox proportional hazards model and then use the  Kalbfleisch-Prentice method to attain survival probabilities.
#Here we will take in a training and a testing set. The function will then train on the training set and return a list containing (1) a 
#matrix of survival curves, where the first column is time values and all following columns are survival probabilities of test subjects and 
#(2) the true death times and event indicator (i.e. time and delta) of the test subjects.

#Input: Survival Dataset post normalization and imputation.
#Output: A list of (1) matrix of survival curves, and (2) the true death times.
############################################################################################################################################
#Library Dependencies
#survival is needed to get survfit and the implimentation of the KP estimator.
library(survival)

CoxPH_KP = function(training, testing){
  coxMod = coxph(Surv(time,delta)~., data = training)
  survivalCurves = survfit(coxMod, testing, type = "kalbfleisch-prentice")
  timePoints = c(0,survivalCurves$time)
  probabilities = rbind(1,survivalCurves$surv)
  curvesToReturn = cbind.data.frame(time = timePoints, probabilities)
  timesAndCens = cbind.data.frame(time = testing$time, delta = testing$delta)
  return(list(curvesToReturn, timesAndCens))  
}




