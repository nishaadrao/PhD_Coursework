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
library(boot)               #for bootstrapping
library(mvtnorm)            #for MVN stuff
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
# Generate random data and simulate
######################################################################