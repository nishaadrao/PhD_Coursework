## ECON675: ASSIGNMENT 5
## Q2: WEAK INSTRUMENTS SIMULATIONS
## Anirudh Yadav 
## 11/19/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(foreach)            #for looping
library(data.table)         #for data manipulation
library(Matrix)             #fast matrix calcs
library(ggplot2)            #for pretty plots
library(sandwich)           #for variance-covariance estimation 
library(xtable)             #for latex tables
library(rdrobust)           #for RD stuff
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
# Input data
######################################################################

data <- as.data.table(read.csv('PhD_Coursework/ECON675/HW6/HeadStart.csv'))

######################################################################
# [2.1] RD Plots and falsification tests
######################################################################





