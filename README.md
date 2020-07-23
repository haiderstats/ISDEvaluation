# ISDEvaluation
This is a repository of the code used for the paper entitled [Effective Ways to Build and Evaluate Individual
Survival Distributions](http://www.jmlr.org/papers/v21/18-772.html), published in the Journal of Machine Learning Research (JMLR). This code can be used to run 6 different survival prediction models (Kaplan-Meier, Cox proportional-hazards, Accelerated Failure Time, Random Survival Forests, Cox Elastic-Net, and Multi-task Logistic Regression) and 5 evaluate those models across 5 different metrics (Concordance, Single/Integrated Brier score, L1-loss, 1-Calibration, and D-Calibration).


If you are interested in using the datasets from the paper you can access the NACD (and subsequently the NACD-Col) data from the [Patient Specific Survival Prediction (PSSP) website](http://pssp.srv.ualberta.ca) under "Public Predictors" or use this [direct download link](http://pssp.srv.ualberta.ca/system/predictors/datasets/000/000/032/original/All_Data_updated_may2011_CLEANED.csv?1350302245). Note that the here the censored bit is flipped from the notation in the paper (CENSORED = 1 implies a censored patient). The data from TCGA can be found by accessing [TCGA's firebrowse website](http://firebrowse.org/), selecting a cancer cohort via dropdown, and downloading the clinical features by selecting the blue bar and choosing the `Clicnical_Pick_Tier1.md5` download link. All details for the High-Dimensional datasets can be found in the paper ["A Multi-Task Learning Formulation for Survival Analysis"](http://dmkd.cs.vt.edu/papers/KDD16.pdf) (left column of page 6) by Li et al. and the dataset downloadable links can be found on their [authors website](http://user.it.uu.se/~liuya610/download.html). Supplemental dataset details and results can be found on [RPubs](http://rpubs.com/haiderstats/ISDEvaluationSupplement) and their basic features and download links are given in the table below.


<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ACC/20160128/gdac.broadinstitute.org_ACC.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz">ACC</a> </th>
    <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BLCA/20160128/gdac.broadinstitute.org_BLCA.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> BLCA</a> </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/CESC/20160128/gdac.broadinstitute.org_CESC.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> CESC</a> </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/CHOL/20160128/gdac.broadinstitute.org_CHOL.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> CHOL </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/COAD/20160128/gdac.broadinstitute.org_COAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> COAD </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/COADREAD/20160128/gdac.broadinstitute.org_COADREAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> COADREAD</a> </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ESCA/20160128/gdac.broadinstitute.org_ESCA.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> ESCA</a> </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/FPPP/20160128/gdac.broadinstitute.org_FPPP.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> FPPP</a> </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/HNSC/20160128/gdac.broadinstitute.org_HNSC.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> HNSC </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/KICH/20160128/gdac.broadinstitute.org_KICH.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> KICH</a> </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/KIPAN/20160128/gdac.broadinstitute.org_KIPAN.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> KIPAN </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/KIRC/20160128/gdac.broadinstitute.org_KIRC.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> KIRC</a> </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/KIRP/20160128/gdac.broadinstitute.org_KIR{.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> KIRP</a> </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LGG/20160128/gdac.broadinstitute.org_LGG.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> LGG </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LIHC/20160128/gdac.broadinstitute.org_LIHC.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> LIHC </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUAD/20160128/gdac.broadinstitute.org_LUAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> LUAD </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUSC/20160128/gdac.broadinstitute.org_LUSC.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> LUSC</a> </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/OV/20160128/gdac.broadinstitute.org_OV.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> OV </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/PAAD/20160128/gdac.broadinstitute.org_PAAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> PAAD </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/PRAD/20160128/gdac.broadinstitute.org_PRAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> PRAD</a> </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/SARC/20160128/gdac.broadinstitute.org_SARC.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> SARC </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/SKCM/20160128/gdac.broadinstitute.org_SKCM.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> SKCM</a> </th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/STAD/20160128/gdac.broadinstitute.org_STAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> STAD </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/STES/20160128/gdac.broadinstitute.org_STES.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> STES </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/THCA/20160128/gdac.broadinstitute.org_THCA.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> THCA </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/UCEC/20160128/gdac.broadinstitute.org_UCEC.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> UCEC </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/UCS/20160128/gdac.broadinstitute.org_UCS.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> UCS </a></th>
   <th style="text-align:center;"> <a href = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/UVM/20160128/gdac.broadinstitute.org_UVM.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"> UVM </a></th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;min-width: 4cm; font-weight: bold;"> N </td>
   <td style="text-align:center;"> 92 </td>
   <td style="text-align:center;"> 409 </td>
   <td style="text-align:center;"> 307 </td>
   <td style="text-align:center;"> 45 </td>
   <td style="text-align:center;"> 456 </td>
   <td style="text-align:center;"> 626 </td>
   <td style="text-align:center;"> 185 </td>
   <td style="text-align:center;"> 38 </td>
   <td style="text-align:center;"> 526 </td>
   <td style="text-align:center;"> 112 </td>
   <td style="text-align:center;"> 939 </td>
   <td style="text-align:center;"> 537 </td>
   <td style="text-align:center;"> 290 </td>
   <td style="text-align:center;"> 513 </td>
   <td style="text-align:center;"> 376 </td>
   <td style="text-align:center;"> 513 </td>
   <td style="text-align:center;"> 498 </td>
   <td style="text-align:center;"> 576 </td>
   <td style="text-align:center;"> 185 </td>
   <td style="text-align:center;"> 499 </td>
   <td style="text-align:center;"> 261 </td>
   <td style="text-align:center;"> 460 </td>
   <td style="text-align:center;"> 436 </td>
   <td style="text-align:center;"> 621 </td>
   <td style="text-align:center;"> 503 </td>
   <td style="text-align:center;"> 546 </td>
   <td style="text-align:center;"> 57 </td>
   <td style="text-align:center;"> 80 </td>
  </tr>
  <tr>
   <td style="text-align:left;min-width: 4cm; font-weight: bold;"> % Censored </td>
   <td style="text-align:center;"> 63.04 </td>
   <td style="text-align:center;"> 55.99 </td>
   <td style="text-align:center;"> 76.55 </td>
   <td style="text-align:center;"> 51.11 </td>
   <td style="text-align:center;"> 77.63 </td>
   <td style="text-align:center;"> 79.39 </td>
   <td style="text-align:center;"> 58.38 </td>
   <td style="text-align:center;"> 81.58 </td>
   <td style="text-align:center;"> 57.6 </td>
   <td style="text-align:center;"> 89.29 </td>
   <td style="text-align:center;"> 75.19 </td>
   <td style="text-align:center;"> 67.04 </td>
   <td style="text-align:center;"> 84.83 </td>
   <td style="text-align:center;"> 75.63 </td>
   <td style="text-align:center;"> 64.89 </td>
   <td style="text-align:center;"> 64.13 </td>
   <td style="text-align:center;"> 56.83 </td>
   <td style="text-align:center;"> 40.45 </td>
   <td style="text-align:center;"> 45.95 </td>
   <td style="text-align:center;"> 98 </td>
   <td style="text-align:center;"> 62.07 </td>
   <td style="text-align:center;"> 51.96 </td>
   <td style="text-align:center;"> 61.01 </td>
   <td style="text-align:center;"> 60.23 </td>
   <td style="text-align:center;"> 96.82 </td>
   <td style="text-align:center;"> 83.33 </td>
   <td style="text-align:center;"> 38.6 </td>
   <td style="text-align:center;"> 71.25 </td>
  </tr>
  <tr>
   <td style="text-align:left;min-width: 4cm; font-weight: bold;"> Number of Features </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 6 </td>
  </tr>
</tbody>
</table>



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
