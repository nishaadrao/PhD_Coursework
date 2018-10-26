## ECON675: ASSIGNMENT 3
## Q3: WHEN BOOTSTRAP FAILS
## Anirudh Yadav 
## 10/26/2018

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
library(boot)               #for bootstrapping
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation


######################################################################
# Q1: Simulate data, run bootsrap
######################################################################
set.seed(123)

N = 1000

# Simulate runif data
X = runif(N,0,1)

# Get max
x.max.obs = max(X)

# Write function for bootrap statistic
boot.stat = function(data, i){
  N*(x.max.obs-max(data[i]))
}

# Run bootsrap with 599 replications
boot.results = boot(data = X, R = 599, statistic = boot.stat)

# NEXT UP PLOT THE BOOT HISTOGRAM WITH EXP(1) OVERLAYED!
