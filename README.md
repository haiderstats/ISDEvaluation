# ISDEvaluation
This is a repository of the code used for the paper entitled "Effective Ways to Build and Evaluate Individual
Survival Distributions". This code can be used to run 6 different survival prediction models (Kaplan-Meier, Cox proportional-hazards, Accelerated Failure Time, Random Survival Forests, Cox Elastic-Net, and Multi-task Logistic Regression) and 5 evaluate those models across 5 different metrics (Concordance, Single/Integrated Brier score, L1-loss, 1-Calibration, and D-Calibration). There is a full tutorial of this code given on  [RPubs](http://rpubs.com/haiderstats/ISDEvaluation), however, we give a miniature example below as well.

Before attempting to run any of the code please make sure you have installed all the required packages: `caret`, `dataPreparation`, `ggplot2`, `reshape2`,`randomForestSRC`, `Rcpp`, `prodlim`, `survival`, `fastcox`, `plyr`,and `dplyr`. Once your working directory is in `ISDEvaluation`, you can run the following

```
source('analysisMaster.R')
#We will use the example data supplied here.
survivalDataset = read.csv(file = "https://raw.githubusercontent.com/haiderstats/ISDEvaluation/master/exampleData.csv")

#Set the event time feature and censor indicator names to  'time' and 'delta'.
names(survivalDataset)[c(1,2)] = c("time", "delta")


ISD = analysisMaster(survivalDataset, numberOfFolds = 5)
```
`analysisMaster()` should print it's progress as it continues to run (set verbose = FALSE to ignore this output). Once the model has finished running (a few minutes to allow RSF to select it's hyper parameters) `ISD` will be a list contains three items: (1) `datasetUsed` -- the survival dataset we passed in post validation and feature selection, (2) `survivalCurves` -- a list of matricies where each column represents the survival probabilities for each patient (one matrix for every model), and (3) `results` -- these are the evaluation results for each model. 



You can plot these survival curves using `plotSurvivalCurves()`; for example to see the first 10 survival curves (for the first 10 rows os ISD$datasetUsed) of the Multi-task Logistic Regression (MTLR) model you can use the following:
```
plotSurvivalCurves(ISD$survivalCurves$MTLR, 1:10)
```

The `results` are a dataframe returning each models performance across the number of folds specified (here 5). Under the default settings these evaluations will be Concordance, Integrated Brier score, single time Brier score evaluated at the median event time, Margin-L1-loss, D-Calibration, and 1-Calibration evaluated at the 10th, 25th, 50th, 75th, and 90th percentiles of event times (for 1-Calibration these correspond to the column names OneCalibration_1, OneCalibration_2,...,OneCalibration_5, respectively).

Please email Humza Haider at hshaider@ualberta.ca with any comments/questions you may have.
