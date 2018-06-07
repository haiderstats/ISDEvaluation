#File Name: analysisMaster.R
#Date Created: May 26, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file should act as the master file to analyze a given dataset with all modeling techniques. As it currently stands we will fill in
#this file as more files get completed. We will need to source all the files and make calls to their respective major functions.

#Input: Survival Dataset, Number of Folds, Time of Interest for One Calibration (possibly a vector), a time point for brier score or 
#a vector of two time points (e.g. c(0,50)) for integrated Brier Score, and the number of bins for D-Calibration. Further, for the 
#L1 and L2 measures we need to specify if we want to use the mean or median.
#Also indicate evaluation metrics and models wanted, default is for all models and evaluation metrics to be computed.
#Additionally, we allow to take in additional model information, e.g. number of survival trees and distribution for AFT.
#Output: A CSV containing averaged evaluation results for each model across the K folds.
############################################################################################################################################
#Dependencies
#Change the path and file names as needed.
#Pre-modelling files:
source("ValidateCleanCV/validateAndClean.R")
source("ValidateCleanCV/createFoldsAndNormalize.R")

#Modeling files:
source("Models/CoxPH_KP.R")
source("Models/KaplanMeier.R")
source("Models/RandomSurvivalForests.R")
source("Models/AcceleratedFailureTime.R")

#As of May 26, I'm not sure if PSSP has an R implementation or if it is in some other language.
#We will leave a commented out R file for now.
#source("Models/PSSP.R")
#Evaluation files:
source("Evaluations/DCalibration.R")
source("Evaluations/OneCalibration.R")
source("Evaluations/Concordance.R")
source("Evaluations/L1Measures.R")
source("Evaluations/BrierScore.R")

#For now we will use an example dataset.
survivalDataset = lung[,-1]
names(survivalDataset)[2] = "delta"
survivalDataset$delta = ifelse(survivalDataset$delta == 2, 1,0)
survivalDataset$ph.ecog = as.factor(survivalDataset$ph.ecog)
survivalDataset$sex = factor(survivalDataset$sex,levels = c(1,2), labels= c("Male","Female"))
#add some NA values to sex for fun
survivalDataset[1:5, 4] = NA


analysisMaster = function(survivalDataset, numberOfFolds,
                          CoxKP = T, KaplanMeier = T, RSFModel = T, AFTModel = T, #Models
                          DCal = T, OneCal = T, Concor = T, L1Measure = T, Brier = T, #Evaluations
                          DCalBins = 10, OneCalTime = NULL,  concordanceTies = "None", #Evaluation args
                          BrierTime = NULL, BrierBasedOnEvents = F, Ltype = "Margin", Llog = F, #Evaluation args
                          typeOneCal = "BucketKM", oneCalBuckets = 10, #Evaluation args
                          AFTDistribution = "weibull", ntree = 1000 #Model args
                          ){
  set.seed(42)
  validatedData = validateAndClean(survivalDataset)
  normalizedData = createFoldsAndNormalize(validatedData, numberOfFolds)
  evaluationResults = data.frame()
  combinedTestResults = list(Cox = list(), KM = list(), AFT = list(), RSF = list())
  for(i in 1:numberOfFolds){
    #Models - We evaluate values to NULL so we can pass them to evaluations, regardless if the models were ran or not.
    coxMod = NULL; kmMod = NULL; rsfMod = NULL; aftMod = NULL;
    training = normalizedData[[1]][[i]]
    testing = normalizedData[[2]][[i]]
    if(CoxKP){
      coxMod = CoxPH_KP(training, testing)
      combinedTestResults$Cox[[i]] = coxMod
    }
    if(KaplanMeier){
      kmMod = KM(training, testing)
      combinedTestResults$KM[[i]] = kmMod
    }
    if(RSFModel){
      rsfMod = RSF(training, testing,ntree = ntree)
      combinedTestResults$RSF[[i]] = rsfMod
    }
    if(AFTModel){
      aftMod = AFT(training, testing, AFTDistribution)
      combinedTestResults$AFT[[i]] = aftMod
    }
    #Evaluations - Note that if evaluations are passed a NULL value they return a NULL.
    DCalResults = NULL;OneCalResults = NULL;ConcResults = NULL;BrierResults = NULL;L1Results = NULL; L2Results = NULL; 
    if(DCal){
      coxDcal = DCalibration(coxMod, DCalBins)
      kmDcal = DCalibration(kmMod, DCalBins)
      rsfDcal = DCalibration(rsfMod, DCalBins)
      aftDcal = DCalibration(aftMod, DCalBins)
      DCalResults = rbind(coxDcal, kmDcal, rsfDcal,aftDcal)
    }
    if(OneCal){
      coxOneCal = OneCalibration(coxMod, OneCalTime, typeOneCal, oneCalBuckets)
      kmOneCal = OneCalibration(kmMod, OneCalTime, typeOneCal, oneCalBuckets)
      rsfOneCal = OneCalibration(rsfMod, OneCalTime, typeOneCal, oneCalBuckets)
      aftOneCal = OneCalibration(aftMod, OneCalTime, typeOneCal, oneCalBuckets)
      OneCalResults = rbind(coxOneCal, kmOneCal, rsfOneCal,aftOneCal)
    }
    if(Concor){
      coxConcCens = Concordance(coxMod, concordanceTies, T)
      kmConcCens = Concordance(kmMod, concordanceTies, T)
      rsfConcCens = Concordance(rsfMod, concordanceTies, T)
      aftConcCens = Concordance(aftMod, concordanceTies, T)
      coxConcUncens = Concordance(coxMod, concordanceTies, F)
      kmConcUncens = Concordance(kmMod, concordanceTies, F)
      rsfConcUncens = Concordance(rsfMod, concordanceTies, F)
      aftConcUncens = Concordance(aftMod, concordanceTies, F)
      ConcCensResults = rbind(coxConcCens, kmConcCens, rsfConcCens, aftConcCens)
      ConcUncensResults = rbind(coxConcUncens, kmConcUncens, rsfConcUncens, aftConcUncens)
    }
    if(Brier){
      print(i)
      coxBrier = BrierScore(coxMod, BrierTime, BrierBasedOnEvents)
      print("coxDone")
      kmBrier = BrierScore(kmMod, BrierTime, BrierBasedOnEvents)
      print("KMDone")
      rsfBrier = BrierScore(rsfMod, BrierTime, BrierBasedOnEvents)
      print("rsfDone")
      aftBrier = BrierScore(aftMod, BrierTime, BrierBasedOnEvents)
      print("aftDone")
      BrierResults = rbind(coxBrier, kmBrier, rsfBrier, aftBrier)
    }
    if(L1Measure){
      coxL1 = L1(coxMod, Ltype, Llog)
      kmL1 = L1(kmMod, Ltype, Llog)
      rsfL1 = L1(rsfMod, Ltype, Llog)
      aftL1 = L1(aftMod, Ltype, Llog)
      L1Results = rbind(coxL1,kmL1,rsfL1,aftL1)
    }
    toAdd = as.data.frame(cbind(DCalResults, OneCalResults, ConcCensResults,ConcUncensResults, BrierResults, L1Results))
    metricsRan = c(DCal, OneCal, Concor,Concor, L1Measure, Brier)
    names(toAdd) = c("DCalibration","OneCalibration","ConcordanceCensored","ConcordanceUncensensore","BrierResults", "L1Results")[metricsRan]
    modelsRan = c(CoxKP, KaplanMeier, RSFModel, AFTModel)
    models = c("CoxKP","Kaplan-Meier","RSF","AFT")[modelsRan]
    toAdd = cbind.data.frame(Model = models,FoldNumer = i, toAdd)
    evaluationResults = rbind.data.frame(evaluationResults, toAdd)
  }
  if(DCal){
    coxCumDcal = DCalibrationCumulative(combinedTestResults$Cox,numBins = numBins)
    kmCumDcal = DCalibrationCumulative(combinedTestResults$KM,numBins = numBins)
    rsfCumDcal = DCalibrationCumulative(combinedTestResults$RSF,numBins = numBins)
    aftCumDcal = DCalibrationCumulative(combinedTestResults$AFT,numBins = numBins)
    DCalCumResults = c(coxCumDcal, kmCumDcal, rsfCumDcal, aftCumDcal)
  }
  if(OneCal){
    coxCum1cal = OneCalibrationCumulative(combinedTestResults$Cox, OneCalTime, typeOneCal, oneCalBuckets)
    kmCum1cal = OneCalibrationCumulative(combinedTestResults$KM, OneCalTime, typeOneCal, oneCalBuckets)
    rsfCum1cal = OneCalibrationCumulative(combinedTestResults$RSF, OneCalTime, typeOneCal, oneCalBuckets)
    aftCum1cal = OneCalibrationCumulative(combinedTestResults$AFT, OneCalTime, typeOneCal, oneCalBuckets)
    OneCalCumResults = c(coxCum1cal, kmCum1cal, rsfCum1cal, aftCum1cal)
  }
  return(list(evaluationResults, DCalCumResults,OneCalCumResults) )
}










