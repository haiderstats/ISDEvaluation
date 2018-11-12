# ISDEvaluation
This is a repository of the code used for the paper entitled "Effective Ways to Build and Evaluate Individual
Survival Distributions". This code can be used to run 6 different survival prediction models (Kaplan-Meier, Cox proportional-hazards, Accelerated Failure Time, Random Survival Forests, Cox Elastic-Net, and Multi-task Logistic Regression) and 5 evaluate those models across 5 different metrics (Concordance, Single/Integrated Brier score, L1-loss, 1-Calibration, and D-Calibration).


If you are interested in using the datasets from the paper you can access the NACD (and subsequently the NACD-Col) data from the [Patient Specific Survival Prediction (PSSP) website](http://pssp.srv.ualberta.ca) under "Public Predictors" or use this [direct download link](http://pssp.srv.ualberta.ca/system/predictors/datasets/000/000/032/original/All_Data_updated_may2011_CLEANED.csv?1350302245). Note that the here the censored bit is flipped from the notation in the paper (CENSORED = 1 implies a censored patient). The data from TCGA can be found by accessing [TCGA's firebrowse website](http://firebrowse.org/), selecting a cancer cohort, and downloading the clinical features. All details for the High-Dimensional datasets can be found in the paper ["A Multi-Task Learning Formulation for Survival Analysis"](http://dmkd.cs.vt.edu/papers/KDD16.pdf) (left column of page 6) by Li et al. Supplemental dataset details and results can be found on [RPubs](http://rpubs.com/haiderstats/ISDEvaluationSupplement).


A full tutorial of the code in this repository can be found on Humza Haider's [RPubs site](http://rpubs.com/haiderstats/ISDEvaluation), however, we also give a brief example below. Before attempting to run any of the code please make sure you have installed all the required packages: `caret`, `dataPreparation`, `ggplot2`, `reshape2`,`randomForestSRC`, `Rcpp`,`RcppArmadillo`, `prodlim`, `survival`, `fastcox`, `plyr`,and `dplyr`. Once your working directory is in `ISDEvaluation`, you can run the following

```
#First load all the code.

source('analysisMaster.R')

#We will use some example data from the "survival" library. For details on this dataset see 'help(lung)'.

survivalDataset = survival::lung

#If you run 'head(survivalDataset' you will see that the time to event feature is already named 'time'
#but the censor bit is named 'status' and has values of 2 for dead and 1 for censored. We change this
#name to 'delta' and set these values to 0 and 1 below.

names(survivalDataset)[c(3)] = c("delta")
survivalDataset$delta = survivalDataset$delta - 1

#Then we can run analysisMaster() to get the evaluation results for all the models.

ISD = analysisMaster(survivalDataset, numberOfFolds = 5)
```

`analysisMaster()` should print it's progress as it continues to run (set verbose = FALSE to ignore this output). Once the model has finished running (a few minutes to allow the different models models to select their hyper parameters) `ISD` will be a list contains three items: (1) `datasetUsed` -- the survival dataset we passed in post validation and feature selection, (2) `survivalCurves` -- a list of matrices where each column represents the survival probabilities for each patient (one matrix for every model), and (3) `results` -- these are the evaluation results for each model. 



You can plot these survival curves using `plotSurvivalCurves()`; for example to see the first 10 survival curves (for the first 10 rows os ISD$datasetUsed) of the Multi-task Logistic Regression (MTLR) model you can use the following:
```
plotSurvivalCurves(ISD$survivalCurves$MTLR, 1:10)
```

The `results` are a dataframe returning each models performance across the number of folds specified (here 5). Under the default settings these evaluations will be Concordance, Integrated Brier score, single time Brier score evaluated at the median event time, Margin-L1-loss, D-Calibration, and 1-Calibration evaluated at the 10th, 25th, 50th, 75th, and 90th percentiles of event times (for 1-Calibration these correspond to the column names OneCalibration_1, OneCalibration_2,...,OneCalibration_5, respectively). Note that while Kaplan-Meier reports values for 1-Calibration, these values are meaningless since Kaplan-Meier assigns equal probabilities for all patients.

Please email Humza Haider at hshaider@ualberta.ca with any comments/questions you may have.
