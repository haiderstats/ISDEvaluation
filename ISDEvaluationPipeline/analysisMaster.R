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
source("Models/MTLR.R")

#Evaluation files:
source("Evaluations/DCalibration.R")
source("Evaluations/OneCalibration.R")
source("Evaluations/Concordance.R")
source("Evaluations/L1Measures.R")
source("Evaluations/BrierScore.R")

#Misc files:
source("FeatureSelection/FeatureSelection.R")

analysisMaster = function(survivalDataset, numberOfFolds,
                          CoxKP = T,CoxKPEN = T, KaplanMeier = T, RSFModel = T, AFTModel = T, MTLRModel =T, #Models
                          DCal = T, OneCal = T, Concor = T, L1Measure = T, Brier = T, #Evaluations
                          DCalBins = 10, OneCalTime = NULL,  concordanceTies = "Risk", #Evaluation args
                          BrierTime = NULL, numBrierPoints = 1000, Ltype = "Margin", Llog = F, #Evaluation args
                          typeOneCal = "BucketKM", oneCalBuckets = 10, survivalPredictionMethod = "Mean", #Evaluation args
                          AFTDistribution = "weibull", ntree = 1000, #Model args,
                          FeatureSelection = T, imputeZero=T # Misc args
                          ){
  validatedData = validateAndClean(survivalDataset, imputeZero)
  if(FeatureSelection)
    validatedData = FeatureSelection(validatedData, type = "UniCox")
  foldsAndNormalizedData = createFoldsAndNormalize(validatedData, numberOfFolds)
  originalIndexing = foldsAndNormalizedData[[1]]
  normalizedData = foldsAndNormalizedData[[2]]
  evaluationResults = data.frame()
  combinedTestResults = list(Cox = list(),CoxEN = list(), KM = list(), AFT = list(), RSF = list(), MTLR = list())
  coxTimes = NULL;coxENTimes = NULL; kmTimes = NULL; rsfTimes = NULL; aftTimes = NULL; mtlrTimes = NULL;
  for(i in 1:numberOfFolds){
    print(Sys.time())
    print(paste("Starting fold",i,"of", numberOfFolds, "total folds."))
    #Models - We evaluate values to NULL so we can pass them to evaluations, regardless if the models were ran or not.
    coxMod = NULL;coxENMod =NULL; kmMod = NULL; rsfMod = NULL; aftMod = NULL; mtlrMod = NULL;
    training = normalizedData[[1]][[i]]
    testing = normalizedData[[2]][[i]]
    print(paste("Beginning model training."))
    if(CoxKP){
      print("Starting Cox Proportional Hazards.")
      coxMod = CoxPH_KP(training, testing)
      if(length(coxMod) ==1){
        combinedTestResults$Cox = list()
        coxTimes = NULL
        CoxKP = F
        if(i > 1)
          evaluationResults = with(evaluationResults,evaluationResults[-which(Model == "CoxKP"),])
      }
      else{
        combinedTestResults$Cox[[i]] = coxMod
        coxTimes = c(coxTimes,coxMod[[1]]$time)
      }
    }
    if(CoxKPEN){
      print("Starting Cox Proportional Hazards - Elastic Net.")
      coxENMod = CoxPH_KP(training, testing,ElasticNet = T)
      combinedTestResults$CoxEN[[i]] = coxENMod
      coxENTimes = c(coxTimes,coxENMod[[1]]$time)
    }
    if(KaplanMeier){
      print("Starting Kaplan Meier.")
      kmMod = KM(training, testing)
      combinedTestResults$KM[[i]] = kmMod
      kmTimes = c(kmTimes,kmMod[[1]]$time)
    }
    if(RSFModel){
      print("Starting Random Survival Forests.")
      rsfMod = RSF(training, testing,ntree = ntree)
      combinedTestResults$RSF[[i]] = rsfMod
      rsfTimes = c(rsfTimes,rsfMod[[1]]$time)
    }
    if(AFTModel){
      print("Starting Accelerated Failure Time.")
      aftMod = AFT(training, testing, AFTDistribution)
      if(length(aftMod)==1){
          combinedTestResults$AFT = list()
          aftTimes = NULL
          AFTModel = F
          if(i >1)
            evaluationResults = with(evaluationResults,evaluationResults[-which(Model == "AFT"),])
        }
      else{
          combinedTestResults$AFT[[i]] = aftMod
          aftTimes = c(aftTimes,aftMod[[1]]$time)
      }
    }
    if(MTLRModel){
      print("Starting Multi-task Logistic Regression (PSSP).")
      mtlrMod = MTLR(training, testing)
      combinedTestResults$MTLR[[i]] = mtlrMod
      mtlrTimes = c(mtlrTimes,mtlrMod[[1]]$time)
    }
    #Evaluations - Note that if evaluations are passed a NULL value they return a NULL.
    DCalResults = NULL;OneCalResults = NULL;ConcCensResults = NULL;ConcUncensResults = NULL;
    BrierResultsInt = NULL;BrierResultsSingle = NULL;L1Results = NULL; L2Results = NULL; 
    if(DCal){
      print("Staring Evaluation: D-Calibration")
      coxDcal = DCalibration(coxMod, DCalBins)
      coxENDcal = DCalibration(coxENMod, DCalBins)
      kmDcal = DCalibration(kmMod, DCalBins)
      rsfDcal = DCalibration(rsfMod, DCalBins)
      aftDcal = DCalibration(aftMod, DCalBins)
      mtlrDcal = DCalibration(mtlrMod, DCalBins)
      DCalResults = rbind(coxDcal,coxENDcal, kmDcal, rsfDcal, aftDcal,mtlrDcal)
    }
    if(OneCal){
      print("Staring Evaluation: One-Calibration")
      coxOneCal = OneCalibration(coxMod, OneCalTime, typeOneCal, oneCalBuckets)
      coxENOneCal = OneCalibration(coxENMod, OneCalTime, typeOneCal, oneCalBuckets)
      kmOneCal = OneCalibration(kmMod, OneCalTime, typeOneCal, oneCalBuckets)
      rsfOneCal = OneCalibration(rsfMod, OneCalTime, typeOneCal, oneCalBuckets)
      aftOneCal = OneCalibration(aftMod, OneCalTime, typeOneCal, oneCalBuckets)
      mtlrOneCal = OneCalibration(mtlrMod, OneCalTime, typeOneCal, oneCalBuckets)
      OneCalResults = rbind(coxOneCal,coxENOneCal, kmOneCal, rsfOneCal,aftOneCal, mtlrOneCal)
    }
    if(Concor){
      print("Staring Evaluation: Concordance")
      coxConcCens = Concordance(coxMod, concordanceTies, T,survivalPredictionMethod)
      coxENConcCens = Concordance(coxENMod, concordanceTies, T,survivalPredictionMethod)
      kmConcCens = Concordance(kmMod, concordanceTies, T,survivalPredictionMethod)
      rsfConcCens = Concordance(rsfMod, concordanceTies, T,survivalPredictionMethod)
      aftConcCens = Concordance(aftMod, concordanceTies, T,survivalPredictionMethod)
      mtlrConcCens = Concordance(mtlrMod, concordanceTies, T,survivalPredictionMethod)
      
      coxConcUncens = Concordance(coxMod, concordanceTies, F,survivalPredictionMethod)
      coxENConcUncens = Concordance(coxENMod, concordanceTies, F,survivalPredictionMethod)
      kmConcUncens = Concordance(kmMod, concordanceTies, F,survivalPredictionMethod)
      rsfConcUncens = Concordance(rsfMod, concordanceTies, F,survivalPredictionMethod)
      aftConcUncens = Concordance(aftMod, concordanceTies, F,survivalPredictionMethod)
      mtlrConcUncens = Concordance(mtlrMod, concordanceTies, F,survivalPredictionMethod)
      
      ConcCensResults = rbind(coxConcCens, coxENConcCens,kmConcCens, rsfConcCens, aftConcCens, mtlrConcCens)
      ConcUncensResults = rbind(coxConcUncens,coxENConcUncens, kmConcUncens, rsfConcUncens, aftConcUncens, mtlrConcUncens)
    }
    if(Brier){
      print("Staring Evaluation: Brier Score")
      coxBrierInt = BrierScore(coxMod, type = "Integrated", numPoints = numBrierPoints)
      coxENBrierInt = BrierScore(coxENMod, type = "Integrated", numPoints = numBrierPoints)
      kmBrierInt = BrierScore(kmMod, type = "Integrated", numPoints = numBrierPoints)
      rsfBrierInt = BrierScore(rsfMod, type = "Integrated",numPoints =  numBrierPoints)
      aftBrierInt = BrierScore(aftMod, type = "Integrated", numPoints = numBrierPoints)
      mtlrBrierInt = BrierScore(mtlrMod, type = "Integrated", numPoints =  numBrierPoints)
      
      coxBrierS = BrierScore(coxMod, type = "Single",singleTime = BrierTime)
      coxENBrierS = BrierScore(coxENMod, type = "Single",singleTime = BrierTime)
      kmBrierS = BrierScore(kmMod, type = "Single",singleTime = BrierTime)
      rsfBrierS = BrierScore(rsfMod, type = "Single",singleTime = BrierTime)
      aftBrierS = BrierScore(aftMod, type = "Single",singleTime = BrierTime)
      mtlrBrierS = BrierScore(mtlrMod, type = "Single",singleTime = BrierTime)
      
      BrierResultsInt = rbind(coxBrierInt,coxENBrierInt, kmBrierInt, rsfBrierInt, aftBrierInt, mtlrBrierInt)
      BrierResultsSingle = rbind(coxBrierS,coxENBrierS, kmBrierS, rsfBrierS, aftBrierS, mtlrBrierS)
    }
    if(L1Measure){
      print("Staring Evaluation: L1 Loss")
      coxL1 = L1(coxMod, Ltype, Llog,survivalPredictionMethod)
      coxENL1 = L1(coxENMod, Ltype, Llog,survivalPredictionMethod)
      kmL1 = L1(kmMod, Ltype, Llog,survivalPredictionMethod)
      rsfL1 = L1(rsfMod, Ltype, Llog,survivalPredictionMethod)
      aftL1 = L1(aftMod, Ltype, Llog,survivalPredictionMethod)
      mtlrL1 = L1(mtlrMod, Ltype, Llog,survivalPredictionMethod)
      
      L1Results = rbind(coxL1,coxENL1,kmL1,rsfL1,aftL1,mtlrL1)
    }
    toAdd = as.data.frame(cbind(DCalResults, OneCalResults, ConcCensResults,ConcUncensResults,
                                BrierResultsInt,BrierResultsSingle, L1Results))
    metricsRan = c(DCal, OneCal, Concor,Concor,Brier,Brier, L1Measure)
    names(toAdd) = c("DCalibration","OneCalibration","ConcordanceCensored","ConcordanceUncensensored",
                     "BrierResultsInt","BrierResultsSingle", "L1Results")[metricsRan]
    modelsRan = c(CoxKP,CoxKPEN, KaplanMeier, RSFModel, AFTModel, MTLRModel)
    models = c("CoxKP","CoxKPEN","Kaplan-Meier","RSF","AFT", "MTLR")[modelsRan]
    toAdd = cbind.data.frame(Model = models,FoldNumer = i, toAdd)
    evaluationResults = rbind.data.frame(evaluationResults, toAdd)
    print(evaluationResults)
  }
  if(DCal){
    print("Staring Evaluation: Cumulative D-Calibration")
    coxCumDcal = DCalibrationCumulative(combinedTestResults$Cox,DCalBins)
    coxENCumDcal = DCalibrationCumulative(combinedTestResults$CoxEN,DCalBins)
    kmCumDcal = DCalibrationCumulative(combinedTestResults$KM,DCalBins)
    rsfCumDcal = DCalibrationCumulative(combinedTestResults$RSF,DCalBins)
    aftCumDcal = DCalibrationCumulative(combinedTestResults$AFT,DCalBins)
    mtlrCumDcal = DCalibrationCumulative(combinedTestResults$MTLR,DCalBins)
    
    DCalCumResults = c(coxCumDcal,coxENCumDcal, kmCumDcal, rsfCumDcal, aftCumDcal, mtlrCumDcal)
    evaluationResults$DCalCumResults = rep(DCalCumResults, numberOfFolds)
  }
  if(OneCal){
    print("Staring Evaluation: Cumulative One-Calibration")
    coxCum1cal = OneCalibrationCumulative(combinedTestResults$Cox, OneCalTime, typeOneCal, oneCalBuckets)
    coxENCum1cal = OneCalibrationCumulative(combinedTestResults$CoxEN, OneCalTime, typeOneCal, oneCalBuckets)
    kmCum1cal = OneCalibrationCumulative(combinedTestResults$KM, OneCalTime, typeOneCal, oneCalBuckets)
    rsfCum1cal = OneCalibrationCumulative(combinedTestResults$RSF, OneCalTime, typeOneCal, oneCalBuckets)
    aftCum1cal = OneCalibrationCumulative(combinedTestResults$AFT, OneCalTime, typeOneCal, oneCalBuckets)
    mtlrCum1cal = OneCalibrationCumulative(combinedTestResults$MTLR, OneCalTime, typeOneCal, oneCalBuckets)
    
    OneCalCumResults = c(coxCum1cal,coxENCum1cal, kmCum1cal, rsfCum1cal,aftCum1cal, mtlrCum1cal)
    evaluationResults$OneCalCumResults = rep(OneCalCumResults, numberOfFolds)
  }
  #We will add some basic information about the dataset.
  evaluationResults$N = nrow(validatedData)
  #Note we subtract 2 to not count `time` and `delta`.
  evaluationResults$NumFeatures = ncol(validatedData) - 2
  evaluationResults$NumFeaturesOneHot = ncol(training) - 2
  evaluationResults$PercentCensored = sum(!validatedData$delta)/nrow(validatedData)
  survivalCurves = getSurvivalCurves(coxTimes,coxENTimes, kmTimes, aftTimes, rsfTimes, mtlrTimes,
                                     CoxKP,CoxKPEN, KaplanMeier, RSFModel, AFTModel, MTLRModel,
                                     combinedTestResults, numberOfFolds,originalIndexing)
  names(survivalCurves) = models
  return(list(datasetUsed = validatedData, survivalCurves = survivalCurves, results = evaluationResults))
}




getSurvivalCurves = function(coxTimes,coxENTimes, kmTimes, aftTimes, rsfTimes, mtlrTimes,
                             CoxKP = T,CoxKPEN=T, KaplanMeier = T, RSFModel = T, AFTModel = T, MTLRModel =T,
                             combinedTestResults, numberOfFolds, originalIndexing){
  originalIndexOrder = order(unname(unlist(originalIndexing)))
  if(!is.null(coxTimes))
    coxTimes = sort(unique(coxTimes))
  if(!is.null(coxENTimes))
    coxENTimes = sort(unique(coxENTimes))
  if(!is.null(kmTimes))
    kmTimes = sort(unique(kmTimes))
  if(!is.null(rsfTimes))
    rsfTimes = sort(unique(rsfTimes))
  if(!is.null(aftTimes))
    aftTimes = sort(unique(aftTimes))
  if(!is.null(mtlrTimes))
    mtlrTimes = sort(unique(mtlrTimes))
  models = c(CoxKP,CoxKPEN, KaplanMeier, AFTModel,RSFModel,MTLRModel)
  allTimes = list(coxTimes,coxENTimes,kmTimes,aftTimes,rsfTimes,mtlrTimes)
  survivalCurves = list()
  count = 0
  for(j in which(models)){
    count =count+1
    fullCurves = data.frame(row.names = 1:length(allTimes[[j]]))
    for(i in 1:numberOfFolds){
      #Index method -> fold -> survival curves
      times = combinedTestResults[[j]][[i]][[1]]$time
      maxTime = max(times)
      curves  = combinedTestResults[[j]][[i]][[1]][,-1]
      timesToEvaluate = setdiff(allTimes[[j]],times)
      fullCurves = cbind.data.frame(fullCurves,sapply(curves,
                                                      function(x){
                                                        curveSpline = splinefun(times,x,method='hyman')
                                                        maxSpline = curveSpline(maxTime)
                                                        curveSplineConstant = function(time){
                                                          timeToEval = ifelse(time > maxTime, maxTime,time)
                                                          toReturn = rep(NA,length(time))
                                                          toReturn[timeToEval== maxTime] = 0
                                                          toReturn[timeToEval !=maxTime] = curveSpline(timeToEval[timeToEval!=maxTime])
                                                          return(toReturn)
                                                        }
                                                        extraPoints =curveSplineConstant(timesToEvaluate)
                                                        toReturn = rep(NA, length(allTimes[[j]]))
                                                        originalIndex = which(!allTimes[[j]] %in% timesToEvaluate)
                                                        newIndex = which(allTimes[[j]] %in% timesToEvaluate)
                                                        toReturn[originalIndex] = x
                                                        toReturn[newIndex] = extraPoints
                                                        return(toReturn)
                                                      }
      ))
    }
    fullCurves = cbind.data.frame(allTimes[j], fullCurves)
    fullCurves[originalIndexOrder]
    colnames(fullCurves) = c("time",1:(ncol(fullCurves)-1))
    survivalCurves[[count]] = fullCurves
  }
  return(survivalCurves)
}






