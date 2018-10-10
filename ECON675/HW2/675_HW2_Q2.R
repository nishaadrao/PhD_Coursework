## ECON675: ASSIGNMENT 2
## Q2: LINEAR SMOOTHERS, CROSS-VALIDATION AND SERIES
## Anirudh Yadav 
## 10/10/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(foreach)            #for looping
library(dplyr)              #for data manipulation
library(data.table)         #for data manipulation
library(ggplot2)            #for pretty plots
library(boot)               #for bootstrapping
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation