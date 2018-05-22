#Date: May 22nd, 2018
#Author: Humza Haider

#Purpose: This file was originally created to impliment an integrated Brier score for the use in the PSSP paper.
#The reason we need to impliment our own Brier score is that other implimentations have assumed inputs of a single 
#survival curve whereas with PSSP and other ISD models we have specific survivor curves for each patient. 

#Given N patients to evaluate, the input to this evaluation metric is N survival curves, the event time for each patient
#and the censoring status of the patient (1 = event occured, 0 = censored).
