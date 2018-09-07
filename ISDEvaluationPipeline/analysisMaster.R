#File Name: analysisMaster.R
#Date Created: May 26, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#This file can act as a master file to analyze a given dataset with all modeling techniques and evaluation metrics.

#Funtion 1: analysisMaster()

#Input: 
#survivalDataset - This is the dataset one wishes to analyze. This must include 'time', 'delta', and at least 1 more feature. No default.
#numberOfFolds - The number of desired cross-validation folds. No default.
#CoxKP, CoxKPEN, KaplanMeier, RSFModel, AFTModel, MTLRModel: Booleans specifying whether or not to run that model. Default is TRUE.
#DCal, OneCal, Concor, L1Measure, BrierSingle, BrierInt: Booleans specifying whether or not to run that evaluation metric. Default is TRUE.
#DCalBins: Number of bins for D-Calibration. Default is 10.
#OneCalTime: An int specifying the time to evaluate 1-Calibration.
#If left as NULL but OneCal = TRUE, then the 10th, 25th, 50th, 75th, and 90th percentiless of all event times are used. Default is NULL.
#concordanceTies: A string ("None", "Time", "Risk","All") indicating how to handle ties in concordance. Default is "Risk".
#SingleBrierTime: The time to evaluate the Brier Score. If left as null, the 50th percentile of all event times is used. Default is NULL.
#IntegratedBrierTimes: A 2 length vector (e.g. c(0,100)) specifying the lower and upper bounds on the integrated Brier score. If NULL then
#the default is 0 as a lower bound and the max event time of the entire dataset is used as an upper bound. Default is NULL.
#numBrierPoints: The number of points to evaluate the integrated Brier score. A simple trapezoidal numerical approximation is used. Default
#is 1000 points.
#Ltype: The type of L1-loss. Must be one of "Uncensored","Hinge", or "Margin". Default is "Margin".
#Llog: A boolean specifying whether or not to use log-L1 metric. Default is FALSE.
#typeOneCal: A string indicating the type of 1-Calibrtion to use. Must be one of "DN" or "Uncensored". Default is "DN".
#oneCalBuckets: An int specifying number of bins for 1-Calibration. Default is 10.
#survivalPredictionMethod: The way in which to estimate average surival times. Must be one of "Mean" or "Median". Default is "Mean".
#AFTDistribution: The distribution to use for AFT, default is "weibull". Must be one of "weibull","exponential","lognormal","gaussian",
#"loglogistic","logistic".
#ntree: The number of trees for RSF. Default is 1000.
#FS: A boolean specifying whether or not to use feature selection. Default is TRUE.
#imputeZero: A boolean specifying whether 0 valued times should be imputed (AFT breaks for 0 valued times). If TRUE then 0 valued times are
#imputed to half the minimum non-zero time. Default is TRUE.

#Output: A list of (3) items:
#(1) datasetUsed: This is the dataset that is actually used post feature selection but pre normalization and imputaiton. datasetUsed
#will have all the patients who had acceptable time and delta values and the features that were selected.
#(2) survivalCurves: This is a list containing the survival curves for all patients for each model that was tested. 
#(3) results: This is a dataframe containing all the evaluation results with specified model and fold number. Additionally the sample size
#feature size, and censoring percnetage are returned. Notice that the feature sizes before and after one hot encoding are returned. 
#If none of the features were factors then NumFeatures should equal NumFeaturesOneHot.

#Note that survivalCurves can be plotted by plotSurvivalCurves().

#Function 2: getSurvivalCurves()
#This is a helper function for analysisMaster(). It simply retrieves all the survival curves from all test folds and fits all test points
#from the test folds so all curves are evaluated on the same time points.
############################################################################################################################################
#Data processing files:
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
source("Plotting/plotSurvivalCurves.R")

analysisMaster = function(survivalDataset, numberOfFolds,
                          CoxKP = T,CoxKPEN = T, KaplanMeier = T, RSFModel = T, AFTModel = T, MTLRModel =T, #Models
                          DCal = T, OneCal = T, Concor = T, L1Measure = T, BrierInt = T, BrierSingle = T, #Evaluations
                          DCalBins = 10, OneCalTime = NULL,  concordanceTies = "Risk", #Evaluation args
                          SingleBrierTime = NULL, IntegratedBrierTimes = NULL, numBrierPoints = 1000, Ltype = "Margin", #Evaluation args
                          Llog = F, typeOneCal = "DN", oneCalBuckets = 10, survivalPredictionMethod = "Mean", #Evaluation args
                          AFTDistribution = "weibull", ntree = 1000, #Model args,
                          FS = T, imputeZero=T # Misc args
                          ){
  validatedData = validateAndClean(survivalDataset, imputeZero)
  if(FS)
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
      coxENTimes = c(coxENTimes,coxENMod[[1]]$time)
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
    if(Concor){
      print("Staring Evaluation: Concordance")
      coxConcUncens = Concordance(coxMod, concordanceTies,survivalPredictionMethod)
      coxENConcUncens = Concordance(coxENMod, concordanceTies,survivalPredictionMethod)
      kmConcUncens = Concordance(kmMod, concordanceTies,survivalPredictionMethod)
      rsfConcUncens = Concordance(rsfMod, concordanceTies,survivalPredictionMethod)
      aftConcUncens = Concordance(aftMod, concordanceTies,survivalPredictionMethod)
      mtlrConcUncens = Concordance(mtlrMod, concordanceTies,survivalPredictionMethod)
      
      ConcUncensResults = rbind(coxConcUncens,coxENConcUncens, kmConcUncens, rsfConcUncens, aftConcUncens, mtlrConcUncens)
    }
    if(BrierInt){
      print("Staring Evaluation: Brier Score- Integrated")
      coxBrierInt = BrierScore(coxMod, type = "Integrated", numPoints = numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      coxENBrierInt = BrierScore(coxENMod, type = "Integrated", numPoints = numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      kmBrierInt = BrierScore(kmMod, type = "Integrated", numPoints = numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      rsfBrierInt = BrierScore(rsfMod, type = "Integrated",numPoints =  numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      aftBrierInt = BrierScore(aftMod, type = "Integrated", numPoints = numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      mtlrBrierInt = BrierScore(mtlrMod, type = "Integrated", numPoints =  numBrierPoints, integratedBrierTimes = IntegratedBrierTimes)
      
      BrierResultsInt = rbind(coxBrierInt,coxENBrierInt, kmBrierInt, rsfBrierInt, aftBrierInt, mtlrBrierInt)
      
    }
    if(BrierSingle){
      print("Staring Evaluation: Brier Score - Single")
      coxBrierSingle = BrierScore(coxMod, type = "Single", singleBrierTime =SingleBrierTime )
      coxENBrierSingle = BrierScore(coxENMod, type = "Single", singleBrierTime =SingleBrierTime )
      kmBrierSingle = BrierScore(kmMod, type = "Single", singleBrierTime =SingleBrierTime )
      rsfBrierSingle = BrierScore(rsfMod, type = "Single", singleBrierTime =SingleBrierTime )
      aftBrierSingle = BrierScore(aftMod, type = "Single", singleBrierTime =SingleBrierTime )
      mtlrBrierSingle = BrierScore(mtlrMod, type = "Single", singleBrierTime =SingleBrierTime )
      
      BrierResultsSingle = rbind(coxBrierSingle,coxENBrierSingle, kmBrierSingle, rsfBrierSingle, aftBrierSingle, mtlrBrierSingle)
      
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
    toAdd = as.data.frame(cbind(ConcUncensResults,
                                BrierResultsInt, BrierResultsSingle,L1Results))
    metricsRan = c(Concor,BrierInt,BrierSingle, L1Measure)
    names(toAdd) = c("ConcordanceUncensensored",
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
    
    numTimes = max(sapply(list(coxCum1cal,coxENCum1cal, kmCum1cal, rsfCum1cal,aftCum1cal, mtlrCum1cal),length))
    
    for(times in 1:numTimes){
      varName = paste("OneCalCumResults",times, sep="")
      assign(varName,c(coxCum1cal[times],coxENCum1cal[times], kmCum1cal[times], rsfCum1cal[times],aftCum1cal[times], mtlrCum1cal[times]))
      evaluationResults[varName] = rep(eval(parse(text=varName)), numberOfFolds)
    }
    print(evaluationResults)
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
  names(survivalCurves) = c("Cox","CoxEN","KM","AFT","RSF","MTLR")[c(CoxKP,CoxKPEN, KaplanMeier, AFTModel,RSFModel, MTLRModel)]
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
                                                          toReturn[timeToEval== maxTime] = max(maxSpline,0)
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
    fullCurves =  fullCurves[originalIndexOrder]
    fullCurves = cbind.data.frame(allTimes[j], fullCurves)
    colnames(fullCurves) = c("time",1:(ncol(fullCurves)-1))
    survivalCurves[[count]] = fullCurves
  }
  return(survivalCurves)
}






