# ISDEvaluation
A repository of the code used for the paper entitled "Effective Ways to Build and Evaluate Individual
Survival Distributions".  

Code was written by Humza Haider (hshaider@ualberta.ca). Note that all code here is written in R (with some Rcpp). Before attempting to run code please make sure you have installed are the required packages: `caret`, `dataPreparation`, `ggplot2`, `reshape2`,`randomForestSRC`, `Rcpp`, `prodlim`, `survival`, `fastcox`, `plyr`,and `dplyr`. The large majority of code is within the Models and Evaluations folders where model and evaluation metric implementations are given respecitvely. There is a master file, `analysisMaster.R` which can run many different options for all models and evaluations. Below we have given an example of how to use this master file and it's function, altough a more detailed tutorial is avaliable on [RPubs](http://rpubs.com/haiderstats/4140630).

Once your working directory is in `ISDEvaluation`, run the following

```
source('analysisMaster.R')
#We will use the lung dataset from the survival library.
survivalDataset = lung

#Set the censor indicator name to delta.
names(survivalDataset)[3] = "delta"

#Change the censor indicator so that delta =1 means the patient died
survivalDataset$delta = survivalDataset$delta - 1

ISD = analysisMaster(survivalDataset, numberOfFolds = 5)
```
`analysisMaster()` should print it's progress as it continues to run. Once the model has finished running (a few minutes to allow RSF to select it's hyper parameters) `ISD` will be a list contains three items: (1) `datasetUsed` - the survival dataset we passed in post validation and feature selection, (2) `survivalCurves` - a list of matricies where each column represents the survival probabilities for each patient (one matrix for every model), and (3) `results` - These are the evaluation results for each model. 

You can plot these survival curves using `plotSurvivalCurves()`; for exampe try:
```
plotSurvivalCurves(ISD$survivalCurves$MTLR, 1:10)
```
for the first 10 patients (see ISD$datasetUsed\[1:10,\]) of the MTLR model.

Please email Humza Haider at hshaider@ualberta.ca with any comments/questions you may have.
