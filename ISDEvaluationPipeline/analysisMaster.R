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
source("Evaluations/L1L2Measures.R")
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
                          CoxKP = T, KaplanMeier = T, RSF = T, AFT = T,
                          DCal = T, OneCal = T, Concordance = T, L1 = T, L2 = T, Brier = T,
                          OneCalTime = c(), BrierTime = c(), DCalBins = 10, LMeasure = "mean"){
  validatedData = validateAndClean(survivalDataset)
  normalizedData = createFoldsAndNormalize(validatedData, numberOfFolds)
  evaluationResults = data.frame()
  for(i in 1:numberOfFolds){
    #Models - We evaluate values to NULL so we can pass them to evaluations, regardless if the models were ran or not.
    coxMod = NULL; kmMod = NULL; rsfMod = NULL; aftMod = NULL;
    training = normalizedData[[1]][[i]]
    testing = normalizedData[[2]][[i]]
    if(CoxKP)
      coxMod = CoxPH_KP(training, testing)
    if(KaplanMeier)
      kmMod = KaplanMeier(training, testing)
    if(RSF)
      rsfMod = RSF(training, testing)
    if(AFT)
      aftMod = AFT(training, testing)
    #Evaluations - Note that if evaluations are passed a NULL value they return a NULL.
    DCalResults = NULL;OneCalResults = NULL;ConcResults = NULL;BrierResults = NULL;L1Results = NULL; L2Results = NULL; 
    if(DCal){
      coxDcal = DCalibration(coxMod, DCalBinds)
      kmDcal = DCalibration(kmMod, DCalBinds)
      rsfDcal = DCalibration(rsfMod, DCalBinds)
      aftDcal = DCalibration(aftMod, DCalBinds)
      DCalResults = rbind(coxDcal, kmDcal, rsfDcal,aftDcal)
    }
    if(OneCal){
      if(length(OneCalTime)==0)
        stop("Please enter a time for One Calibration or change OneCal to FALSE")
      coxOneCal = OneCalibration(coxMod, OneCalTime)
      kmOneCal = OneCalibration(kmMod, OneCalTime)
      rsfOneCal = OneCalibration(rsfMod, OneCalTime)
      aftOneCal = OneCalibration(aftMod, OneCalTime)
      OneCalResults = rbind(coxOneCal, kmOneCal, rsfOneCal,aftOneCal)
    }
    if(Concordance){
      coxConc = ConcordanceEval(coxMod)
      kmConc = ConcordanceEval(kmMod)
      rsfConc = ConcordanceEval(rsfMod)
      aftConc = ConcordanceEval(aftMod)
      ConcResults = rbind(coxConc, kmConc, rsfConc, aftConc)
    }
    if(Brier){
      if(length(BrierTime)==0)
        stop("Please enter a time for the Brier score or change Brier to FALSE")
      coxBrier = BrierScore(coxMod, BrierTime)
      kmBrier = BrierScore(kmMod, BrierTime)
      rsfBrier = BrierScore(rsfMod, BrierTime)
      aftBrier = BrierScore(aftMod, BrierTime)
      BrierResults = rbind(coxBrier, kmBrier, rsfBrier, aftBrier)
    }
    if(L1){
      coxL1 = L1Eval(coxMod, LMeasure)
      kmL1 = L1Eval(kmMod, LMeasure)
      rsfL1 = L1Eval(rsfMod, LMeasure)
      aftL1 = L1Eval(aftMod, LMeasure)
      L1Results = rbind(coxL1,kmL1,rsfL1,aftL1)
    }
    if(L2){
      coxL2 = L2Eval(coxMod, LMeasure)
      kmL2 = L2Eval(kmMod, LMeasure)
      rsfL2 = L2Eval(rsfMod, LMeasure)
      aftL2 = L2Eval(aftMod, LMeasure)
      L2Results = rbind(coxL2,kmL2,rsfL2,aftL2)
    }
    toAdd = as.data.frame(cbind(DCalResults, OneCalResults, ConcResults, BrierResults, L1Results, L2Results))
    metricsRan = c(DCal,OneCal,Concordance, Brier,L1,L2)
    names(toAdd) = c("DCalibration","OneCalibration","Concordance","BrierResults", "L1Results","L2Results")[metricsRan]
    toAdd = cbind.data.frame(FoldNumer = i, toAdd)
    evaluationResults = rbind.data.frame(evaluationResults, toAdd)
  }
  return(evaluationResults)
}










