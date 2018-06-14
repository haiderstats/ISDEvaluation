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

MTLR = function(training, testing){
  executablesPath = "Models/AdditionalMTLRFiles/"
  originalWd = getwd()
  setwd(paste(originalWd,"/",executablesPath,sep=""))
  #Write csv files to be called by CovertDataFiles to make the correct format of input for MTLR.
  write.csv(training, paste("Data/training.csv",sep=""), row.names = F)
  write.csv(testing, paste("Data/testing.csv",sep=""), row.names = F)
  #We have hardcoded the paths to the data sets. Future use of this code may need to change these 
  system("java -cp ./ ConvertDataFiles convert2MTLR Data/training.csv Data/training.mtlr FlipCensoredBit")
  system("java -cp ./ ConvertDataFiles convert2MTLR Data/testing.csv Data/testing.mtlr FlipCensoredBit")
  system("./mtlr_opt -i Data/training.mtlr",ignore.stdout = T)
  system("rm Ptrain1")
  system("rm Pmodel1")
  system("./mtlr_test -i Data/testing.mtlr -s Data/training.mtlr -o ./fold1_modelfile > Output/MTLR_output.txt")
  system("rm CI_log")
  timePoints = unlist((unname(read.table("fold1_modelfile",skip = 1,sep = ",",nrows = 1))))
  #There are more survival probabilities than survival time points. Further, the previous writer of code (Fatima) used the last time point
  #squared divided by the second to last time point. I'm guessing that that this was done using some information not known to me at 
  #this time so until further notice we will do the same.
  lastTimePoint = round(timePoints[length(timePoints)]^2/timePoints[length(timePoints) -1])
  #Add the last time point and a 0 time point.
  timePoints = c(0,timePoints, lastTimePoint)
  testingPoints = read.table("Output/MTLR_output.txt")
  #Replace original working directory.
  setwd(originalWd)
  #the first 4 columns are the true time of death, 1- censoring status, and 2 different averaged survival times.
  trueDeathTimes = testingPoints[,1]
  censorStatus = 1 -testingPoints[,2]
  #We don't want the averages and the last column is some evaluation of survival probability at the true time of death. We will 
  #discard this since later on we have a method used for every survival curve. Further we minus 1 because we dont want to include the 
  #0th time point. So we have ignore first 4 columns through the number of time points, subtracting 1 for the 0th time point.
  survivalProbabilities = testingPoints[5:(4+length(timePoints)-1)]
  #Survival probabilities were read in a factors and contain commas. We will clean this out by turning them into character vectors and 
  #trimming commas.
  survivalProbabilities = apply(survivalProbabilities,c(1,2),as.character)
  survivalProbabilities = apply(survivalProbabilities,c(1,2), function(x) as.numeric(gsub(",","",x)))
  #Some of the last survival probabilities were negative so we turned these to zero. There were values like -1.2239e-17 so effectively 
  #zero anyways. Additionally we need to transpose the survival probabilities to match up with the survival time estimates and then
  #add a survival probability of 1 to the 0th time point.
  survivalProbabilities = rbind(1,t(apply(survivalProbabilities, c(1,2), function(x) ifelse(x < 0,0,x))))
  curvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilities) 
  timesAndCensTest = cbind.data.frame(time = trueDeathTimes, delta = censorStatus)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  return(list(curvesToReturn, timesAndCensTest,timesAndCensTrain))  
}
  
  
  
  
  
  
  
  
  
  
  


