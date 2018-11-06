#File Name: MTLR2.R
#Date Created: October 23, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:


#Input 1: training set
#Output:
############################################################################################################################################

library(Rcpp)

sourceCpp("Models/AdditionalMTLRFiles/C/MTLRRCpp.cpp")
source("ValidateCleanCV/createFoldsAndNormalize.R")

MTLR = function(training, testing, C1 = NULL, numFolds = 5){
  if(is.null(C1)){
    C1 = internalCV(training, numFolds)
  }
  
  ordTrain = order(training$delta)
  ordTest = order(testing$delta)
  
  training = training[ordTrain,]
  testing = testing[ordTest,]
  
  m = floor(sqrt(nrow(training))+1)
  quantileVals = seq(0,1,length.out = m+2)[-c(1,m+2)]
  timePoints = unname(quantile(training$time, quantileVals))
  timePoints = timePoints[!duplicated(timePoints)]
  d = as.matrix(training[,-c(1,2)])
  dAsZero = matrix(0,ncol = ncol(d), nrow = nrow(d))
  yval = matrix(1 -Reduce(c,Map(function(x) training$time > timePoints[x],1:length(timePoints))), ncol = nrow(training), byrow = T)
  biasPar = optim(par = rep(0,length(timePoints)*(ncol(d) +1)),fn = objValC_Cens_LogTrick,gr = gradC_Cens_LogTrick, yval = yval, featureVal = dAsZero,C1=C1, delta = sort(training$delta), 
                  method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45041e11))
  
  allParamsUnc = optim(par = biasPar$par,fn = objValC_Cens_LogTrick,gr = gradC_Cens_LogTrick, yval = yval, featureVal = d,C1=C1,delta = rep(1,nrow(training)),
                    method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45041e11))
  
  allParams = optim(par = allParamsUnc$par,fn = objValC_Cens_LogTrick,gr = gradC_Cens_LogTrick, yval = yval, featureVal = d,C1=C1,delta = sort(training$delta),
                    method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45041e11))
  
  survivalProbabilitiesTrain =  predictC(allParams$par, d)
  survivalProbabilitiesTest =  predictC(allParams$par, as.matrix(testing[,-c(1,2)]))
  #Issue due to machine precision, we get survival probabilities of 1+e-16. So here we adjust for that.
  survivalProbabilitiesTrain[survivalProbabilitiesTrain > 1] = 1
  survivalProbabilitiesTest[survivalProbabilitiesTest > 1] = 1
  
  survivalProbabilitiesTrain = survivalProbabilitiesTrain[,order(ordTrain)]
  survivalProbabilitiesTest = survivalProbabilitiesTest[,order(ordTest)]
  training = training[order(ordTrain),]
  testing = testing[order(ordTest),]
  
  if(!(0 %in% timePoints)){
    timePoints = c(0,timePoints)
    survivalProbabilitiesTest = rbind(1,survivalProbabilitiesTest)
    survivalProbabilitiesTrain = rbind(1,survivalProbabilitiesTrain)
  }
  
  testCurvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilitiesTest) 
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  trainingCurvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilitiesTrain) 
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}


internalCV = function(training, numFolds =5){
 foldedData = createFoldsOfData(training, numFolds)[[2]] 
 
  resultsMatrix = matrix(rep(0,numFolds*7), ncol = 7,nrow =5) #7 vals of C1
  for(i in 1:numFolds){
    trainingFold = foldedData[[1]][[i]]
    testingFold = foldedData[[2]][[i]]
    
    trainingFold = trainingFold[order(trainingFold$delta),]
    testingFold = testingFold[order(testingFold$delta),]
    
    m = floor(sqrt(nrow(training))+1)
    quantileVals = seq(0,1,length.out = m+2)[-c(1,m+2)]
    timePoints = unname(quantile(trainingFold$time, quantileVals))
    timePoints = timePoints[!duplicated(timePoints)]
    d = as.matrix(trainingFold[,-c(1,2)])
    dAsZero = matrix(0,ncol = ncol(d), nrow = nrow(d))
    yval = matrix(1 -Reduce(c,Map(function(x) trainingFold$time > timePoints[x],1:length(timePoints))), ncol = nrow(trainingFold), byrow = T)
    
    #Train biases first and then all parameters.
    resultVec = c()
    for(C1 in c(0.001,0.01,0.1, 1, 10,100,1000)){
      biasPar = optim(par = rep(0,length(timePoints)*(ncol(d) +1)),fn = objValC_Cens_LogTrick,gr = gradC_Cens_LogTrick, yval = yval, featureVal = dAsZero,C1=C1, delta = sort(trainingFold$delta), 
                      method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45041e11))
      
      allParamsUnc = optim(par = biasPar$par,fn = objValC_Cens_LogTrick,gr = gradC_Cens_LogTrick, yval = yval, featureVal = d,C1=C1,delta = rep(1,nrow(trainingFold)),
                           method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45041e11))
      
      allParams = optim(par = allParamsUnc$par,fn = objValC_Cens_LogTrick,gr = gradC_Cens_LogTrick, yval = yval, featureVal = d,C1=C1,delta = sort(trainingFold$delta),
                        method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45041e11))
      
      resultVec = c(resultVec,avgLogLikLoss(allParams$par,testingFold,timePoints))
    }
    resultsMatrix[i,] = resultVec
  }
  meanResults = apply(resultsMatrix, 2, mean)
  bestC1 = c(0.001,0.01,0.1, 1, 10, 100,1000)[which.min(meanResults)]
  return(bestC1)
}
  
avgLogLikLoss = function(params, dat, timePoints){
  d = as.matrix(dat[,-c(1,2)])
  survivalCurves = predictC(params,d)
  NCens = sum(1-  dat$delta)
  #Cens
  logloss = 0
  censorTimes = dat$time[1:NCens]
  x = sapply(seq_along(censorTimes),
             function(index) predictProbabilityFromCurve(survivalCurves[,index],
                                                         timePoints,
                                                         censorTimes[index]))
  logloss = logloss - sum(log(x + 1e-5))
  
  deathTimes = dat$time[(NCens+1):nrow(dat)]
  uncenSurvival = survivalCurves[,(NCens+1):nrow(dat)]
  uncenSurvival = rbind(1,uncenSurvival,0)
  pmfProbs = -diff(uncenSurvival)
  indexRow = sapply(deathTimes, function(x) findInterval(x, timePoints)) + 1
  indexCol = 1:(nrow(dat) - NCens)
  indexMat = matrix(c(indexRow,indexCol),ncol = 2)
  probs = pmfProbs[indexMat]
  logloss = logloss  - sum(log(probs))
  
  return(logloss/nrow(dat))
}
  
  
  