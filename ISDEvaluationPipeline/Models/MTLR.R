#File Name: MTLR.R
#Date Created: June 14, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file is used to run MTLR to generate individual survival curves. All the MTLR code has been turned into an executable generated from
#the C++ code it was written in. For this reason we use system commands to executa the MTLR Code.
#The function will then train on the training set and return a list containing (1) a matrix of survival curves, where the first column 
#is time values and all following columns are survival probabilities of test subjects and  (2) the true death times and event indicator
#(i.e. time and delta) of the test subjects.
#Input 1: Survival Dataset post normalization and imputation.
#Input 2: Number of Trees
#Output: A list of (1) matrix of survival curves, (2) the true death times and censor status of the TESTING SET, and (3) the true death
#times and censor status of the TRAINING set.
############################################################################################################################################

MTLR = function(training, testing,addLastTimePoint=F){
  #The idea here is to have the working directory sitting in ISDEvaluationPipeline. We move the working directory into the folder
  #with executables for ease of execution and move back to the original working directory before exiting the function.
  executablesPath = "Models/AdditionalMTLRFiles/"
  originalWd = getwd()
  setwd(paste(originalWd,"/",executablesPath,sep=""))
  #Write csv files to be called by CovertDataFiles to make the correct format of input for MTLR.
  write.csv(training, paste("training.csv",sep=""), row.names = F)
  write.csv(testing, paste("testing.csv",sep=""), row.names = F)
  system2("java", args=c('-cp ./',"ConvertDataFiles", "convert2MTLR", "training.csv", "training.mtlr", "FlipCensoredBit"))
  system2("java", args=c('-cp ./',"ConvertDataFiles", "convert2MTLR", "testing.csv", "testing.mtlr", "FlipCensoredBit"))
  system2("./mtlr_opt",args=c("-i", "training.mtlr"),stdout = FALSE)
  system2("./mtlr_test", args=c("-i", "testing.mtlr", "-s","training.mtlr","-p", "-o", "./fold1_modelfile > MTLR_output.txt"))
  times = unlist((unname(read.table("fold1_modelfile",skip = 1,sep = ",",nrows = 1))))
  if(addLastTimePoint){
    #It appears the last survival probability is always appearing to be zero. We need to add some additional time point to the end to account
    #for this so we choose the value of the last time point  squared divided by the second to last time point. 
    lastTimePoint = round(times[length(times)]^2/times[length(times) -1])
    #Add the last time point and a 0 time point if needed.
    if(0 %in% times){
      timePoints = c(times,lastTimePoint)
    } else{
      timePoints = c(0,times, lastTimePoint)
    } 
  }
  else{
    if(0 %in% times){
      timePoints = times
    } else{
      timePoints = c(0,times)
    } 
  }

  testingPoints = read.table("MTLR_output.txt")
  
  #######Get training data
  system2("./mtlr_test", args=c("-i", "training.mtlr", "-s","training.mtlr", "-p","-o", "./fold1_modelfile > MTLR_output.txt"))
  timesTrain = unlist((unname(read.table("fold1_modelfile",skip = 1,sep = ",",nrows = 1))))
  if(addLastTimePoint){
    #It appears the last survival probability is always appearing to be zero. We need to add some additional time point to the end to account
    #for this so we choose the value of the last time point  squared divided by the second to last time point. 
    lastTimePoint = round(timesTrain[length(timesTrain)]^2/timesTrain[length(timesTrain) -1])
    #Add the last time point and a 0 time point if needed.
    if(0 %in% timesTrain){
      timePointsTrain = c(timesTrain,lastTimePoint)
    } else{
      timePointsTrain = c(0,timesTrain, lastTimePoint)
    } 
  }
  else{
    if(0 %in% timesTrain){
      timePointsTrain = timesTrain
    } else{
      timePointsTrain = c(0,timesTrain)
    } 
  }
  trainingPoints = read.table("MTLR_output.txt")
  #######
  
  #Clean up directory:
  system2("rm",args=c("fold*_modelfile", "tmp_*","CI_log","*.csv","*.mtlr","*.txt"))
  #Replace original working directory.
  setwd(originalWd)
  #the first 4 columns are the true time of death, 1- censoring status, and 2 different averaged survival times.
  trueDeathTimes = testingPoints[,1]
  censorStatus = 1 -testingPoints[,2]
  
  trueDeathTimesTrain = trainingPoints[,1]
  censorStatusTrain = 1 -trainingPoints[,2]
  #We don't want the averages (first 3 columns) and the last column is some evaluation of survival probability at the true time of death.
  #We will discard this since later on we have a method used for every survival curve. Further we minus 1 because we dont want to include the 
  #0th time point, i.e. since we use length(timePoints) we need to subtract an extra value it since we included a 0.
  #So we have ignore first 4 columns through the number of time points, subtracting 1 for the 0th time point (if we included a zero).
  if(0 %in% times){
      survivalProbabilities = testingPoints[5:(4+length(timePoints))]
      survivalProbabilitiesTrain = trainingPoints[5:(4+length(timePointsTrain))]
      
  } else{
      survivalProbabilities = testingPoints[5:(4+length(timePoints)-1)]
      survivalProbabilitiesTrain = trainingPoints[5:(4+length(timePointsTrain)-1)]
    }
  #Survival probabilities were read in a factors and contain commas. We will clean this out by turning them into character vectors and 
  #trimming commas.
  survivalProbabilities = apply(survivalProbabilities,c(1,2),as.character)
  survivalProbabilities = apply(survivalProbabilities,c(1,2), function(x) as.numeric(gsub(",","",x)))
  
  survivalProbabilitiesTrain = apply(survivalProbabilitiesTrain,c(1,2),as.character)
  survivalProbabilitiesTrain = apply(survivalProbabilitiesTrain,c(1,2), function(x) as.numeric(gsub(",","",x)))
  #Some of the last survival probabilities were negative so we turned these to zero. There were values like -1.2239e-17 so effectively 
  #zero anyways. Additionally we need to transpose the survival probabilities to match up with the survival time estimates and then
  #add a survival probability of 1 to the 0th time point if there was no original 0 time point.
  if(0 %in% times){
    survivalProbabilities = t(apply(survivalProbabilities, c(1,2), function(x) ifelse(x < 0,0,x)))
    survivalProbabilitiesTrain = t(apply(survivalProbabilitiesTrain, c(1,2), function(x) ifelse(x < 0,0,x)))
    
  } else {
    survivalProbabilities = rbind(1,t(apply(survivalProbabilities, c(1,2), function(x) ifelse(x < 0,0,x))))
    survivalProbabilitiesTrain = rbind(1,t(apply(survivalProbabilitiesTrain, c(1,2), function(x) ifelse(x < 0,0,x))))
    
  }
  curvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilities) 
  timesAndCensTest = cbind.data.frame(time = trueDeathTimes, delta = censorStatus)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  trainingCurvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilities) 
  return(list(TestCurves = curvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
  
}
