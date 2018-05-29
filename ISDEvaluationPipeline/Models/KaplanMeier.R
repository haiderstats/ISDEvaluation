#File Name: KaplanMeier.R
#Date Created: May 28, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file is used to run the Kaplan Meier to attain individual survival probabilities. While typically people use the implimentation given
#in the survival package, we instead use prodlim as prodlim has a predict function built in whereas survival does not.
#Here we will take in a training and a testing set. The function will then train on the training set and return a list containing (1) a 
#matrix of survival curves, where the first column is time values and all following columns are survival probabilities of test subjects and 
#(2) the true death times and event indicator (i.e. time and delta) of the test subjects.

#Input: Survival Dataset post normalization and imputation.
#Output: A list of (1) matrix of survival curves, and (2) the true death times.
############################################################################################################################################
#Library Dependencies
#prodlim has a Kaplan Meier implementation that includes a predict function.
library(prodlim)

KM = function(training, testing){
  kmMod = prodlim(Surv(time,delta)~1, data = training)
  timesToPredict = c(0,sort(unique(training$time)))
  probabilities = predict(kmMod,timesToPredict)
  curvesToReturn = cbind.data.frame(time = timesToPredict,matrix(rep(probabilities,nrow(testing)),ncol = nrow(testing)))
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  return(list(curvesToReturn, timesAndCensTest,timesAndCensTrain))   
}



