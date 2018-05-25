#File Name: validateAndClean.R
#Date: May 25, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
#The purpose of this file is to take in a survival dataset and validate that it has all parameters required for analysis. Specifically, 
#we need to validate that the dataset has two columns names time and delta which are our time and event indicator variables respectively.
#We will also check that these variables do not have things such as nonnegative times, and validate that delta only takes on values of
#0 and 1, i.e. right censor and event indicators - we do not plan to handle any other type of censoring. From there we will do basic
#cleaning procedures such as removing heavily missing columns (greater than 25% of observations).

#Input: Survival Dataset
#Output: Survival Dataset which has been validated and had minor cleaning.
############################################################################################################################################













