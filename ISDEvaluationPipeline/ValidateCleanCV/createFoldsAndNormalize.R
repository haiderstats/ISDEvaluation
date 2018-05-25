#File Name: createFoldsAndNormalize.R
#Date: May 25, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#The purpose of this file is to take in a dataset which has already been validated and had some (initial) cleaning. Here we expand upon
#that by creating cross validation folds and doing basic mean imputation and normalizing the data. Specifically, we will do a one hot 
#encoding for factor variables and subtract the mean and divide by the standard deviation for numeric variables. We make sure not to do
#this for the 'time' and 'delta' variables.

#Input #1: Survival Dataset which has been validated and had initial cleaning
#Input #2: Desired number of folds.
#
#Output: List containing two lists, a training list and a testing list. Each of these inner lists will contain K datasets where K is the
#number of desired folds.
#Diagram of Output:
#                             
#                                  |.--> Training Dataset #1
#                                  |.--> Training Dataset #2
#              |---> Training List |
#              |                   |.--> Training Dataset #(K-1)
#              |                   |.--> Training Dataset #K
#Starting List |
#              |                   |.--> Testing Dataset #1
#              |                   |.--> Testing Dataset #2
#              |--->  Testing List |
#                                  |.--> Testing Dataset #(K-1)
#                                  |.--> Testing Dataset #K
############################################################################################################################################























